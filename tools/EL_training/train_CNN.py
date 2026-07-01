import os
import torch
import torch.onnx
import matplotlib
import matplotlib.pyplot as plt
from train_common_EL import get_device, setup_problem, train_and_eval

# Is the active matplotlib backend interactive? When the pipeline runs headless
# it sets MPLBACKEND=Agg (non-interactive) so training never blocks on a GUI;
# under such a backend plt.ion()/plt.pause()/plt.show() only emit
# "FigureCanvasAgg is non-interactive" warnings. Detect the backend once and
# skip the interactive calls when it can't show a window — the figure is saved
# to a PNG instead so the loss history is still captured.
_INTERACTIVE_BACKEND = matplotlib.get_backend().lower() not in {
    'agg', 'pdf', 'ps', 'svg', 'cairo', 'template'
}

# ─────────────────────────────────────────────
# Configuration  (single definition of each)
# ─────────────────────────────────────────────
# The tensor filenames and the output model basename can be tagged per test via
# the environment (set by tools/EL_training/run_element_learning.sh) so that the
# data / models from different tests and grids never overwrite one another.
# Unset ⇒ the historical defaults.
inputfile      = os.environ.get('EL_INPUT_TENSOR',  'input_tensor.csv')
outputfile     = os.environ.get('EL_OUTPUT_TENSOR', 'output_tensor.csv')
dataname       = os.environ.get('EL_DATANAME',      'JX_NN')

num_epochs          = 1000
b_size              = 256
ln_rate             = 1e-4
grad_clip_norm      = 0.5
early_stop_patience = 200       # ← only defined once

modelname      = dataname + '_model.pth'
best_modelname = dataname + '_best_model.pth'
training_error = dataname + '_error.mat'
onnx_name      = dataname + '_model.onnx'

# ─────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────
device = get_device()
dataloader_train, dataloader_test, model, criterion = \
    setup_problem(inputfile, outputfile, b_size, device)

# ── Weight initialisation (Kaiming for ReLU networks) ───────────────────────
def init_weights(m):
    if isinstance(m, torch.nn.Linear):
        torch.nn.init.kaiming_normal_(m.weight, nonlinearity='relu')
        torch.nn.init.zeros_(m.bias)

model.apply(init_weights)

# ── Optimiser + LR scheduler ─────────────────────────────────────────────────
optimizer = torch.optim.Adam(model.parameters(), lr=ln_rate)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer,
    mode      = 'min',
    factor    = 0.5,
    patience  = 50,
    min_lr    = 1e-8,
    threshold = 1e-6,
)

# ─────────────────────────────────────────────
# Real-time plot setup
# ─────────────────────────────────────────────
if _INTERACTIVE_BACKEND:
    plt.ion()
fig, ax = plt.subplots(figsize=(10, 6))
line_train, = ax.plot([], [], 'b-', label='Train Loss', linewidth=1.5)
line_test,  = ax.plot([], [], 'r-', label='Test Loss',  linewidth=1.5)
ax.set_yscale('log')
ax.set_xlabel('Epochs')
ax.set_ylabel('Loss (Log Scale)')
ax.set_title(f'Real-time Training Monitor: {dataname}')
ax.legend()
ax.grid(True, which="both", ls="-", alpha=0.3)
plt.tight_layout()

# ── Live-plot callback — called from inside train_and_eval every epoch ───────
def update_plot(epoch_idx, train_hist, test_hist):
    # Live updates only make sense with an interactive backend; skip them
    # (and the warning-emitting plt.pause) when running headless.
    if not _INTERACTIVE_BACKEND:
        return
    if epoch_idx % 10 == 0:
        epochs_range = list(range(len(train_hist)))
        line_train.set_data(epochs_range, train_hist)
        line_test.set_data(epochs_range,  test_hist)
        ax.relim()
        ax.autoscale_view()
        fig.canvas.draw()
        fig.canvas.flush_events()
        plt.pause(0.001)

# ─────────────────────────────────────────────
# Training
# ─────────────────────────────────────────────
print(f"Starting training on {device} for up to {num_epochs} epochs...")

model, optimizer, stopped_epoch, train_loss_history, test_loss_history = \
    train_and_eval(
        model                = model,
        optimizer            = optimizer,
        dataloader_train     = dataloader_train,
        dataloader_test      = dataloader_test,
        criterion            = criterion,
        num_epochs           = num_epochs,
        modelname            = modelname,
        best_modelname       = best_modelname,
        training_error       = training_error,
        scheduler            = scheduler,
        early_stop_patience  = early_stop_patience,
        grad_clip_norm       = grad_clip_norm,
        plot_callback        = update_plot,   # ← drives the live plot
    )

# ─────────────────────────────────────────────
# Final plot — force one last redraw + info box
# ─────────────────────────────────────────────
epochs_range = list(range(len(train_loss_history)))
line_train.set_data(epochs_range, train_loss_history)
line_test.set_data(epochs_range,  test_loss_history)
ax.relim()
ax.autoscale_view()

hparam_text = (
    f"num_epochs          = {num_epochs}\n"
    f"b_size              = {b_size}\n"
    f"ln_rate             = {ln_rate:.2e}\n"
    f"grad_clip_norm      = {grad_clip_norm}\n"
    f"early_stop_patience = {early_stop_patience}"
)
ax.text(
    0.98, 0.97,
    hparam_text,
    transform           = ax.transAxes,
    fontsize            = 9,
    verticalalignment   = 'top',
    horizontalalignment = 'right',
    bbox                = dict(boxstyle='round,pad=0.4',
                               facecolor='white',
                               alpha=0.8,
                               edgecolor='gray'),
    fontfamily          = 'monospace',
)
fig.canvas.draw()
if _INTERACTIVE_BACKEND:
    fig.canvas.flush_events()

# Always persist the loss-history figure so it is captured even headless.
plot_name = dataname + '_training.png'
fig.savefig(plot_name, dpi=150)
print(f"Saved training plot : {plot_name}")

# ─────────────────────────────────────────────
# ONNX export — always from CPU + best weights
# ─────────────────────────────────────────────
cpu_device = torch.device("cpu")
best_ckpt  = torch.load(best_modelname, map_location=cpu_device)
model_cpu  = model.to(cpu_device)
model_cpu.load_state_dict(best_ckpt['model_state_dict'])
model_cpu.eval()

dummy_input = next(iter(dataloader_train))[0][:1].to(cpu_device)
torch.onnx.export(
    model_cpu,
    dummy_input,
    onnx_name,
    export_params       = True,
    opset_version       = 17,
    do_constant_folding = True,
    input_names         = ['input'],
    output_names        = ['output'],
    dynamic_axes        = {
        'input':  {0: 'batch_size'},
        'output': {0: 'batch_size'},
    }
)
print(f"Saved ONNX model : {onnx_name}  (exported from best checkpoint)")

# ─────────────────────────────────────────────
# Finalise plot
# ─────────────────────────────────────────────
if _INTERACTIVE_BACKEND:
    plt.ioff()
    plt.show()
print("Done.")
