"""
train_nn.py  —  MLP/CNN trainer for the Element-Learning surrogate
                WITH input AND output standardisation (the key fix).

Why this exists
---------------
The plain trainer fed the raw â feature and raw T^{ie} target straight into the
network. Two things then make the loss plateau at "predict the mean":

  1. The â feature is heterogeneous (positive diagonals O(1..5), signed
     off-diagonals) and, with the old sampler, carried a nuisance amplitude.
  2. The 600-dim T^{ie} target has entries of very different magnitudes, so an
     un-normalised L1/L2 loss is dominated by the few large entries and ignores
     the rest.

This script standardises BOTH the inputs and the outputs (zero mean / unit
variance, computed from the TRAIN split) before training, and BAKES the
standardisation (and the output de-standardisation) into the exported ONNX. The
Julia side therefore keeps feeding the RAW â feature and reading the RAW T^{ie}
— no change needed in Jexpresso.

It reuses your train_common_EL building blocks (setup_problem / train_and_eval)
so the model architecture and training loop are unchanged; only the data scaling
and the ONNX export wrapper are added. Set `standardize = False` to recover the
original behaviour.
"""

import os
import sys
# Make the helper modules (train_common_EL, IO_EL, NN_RFRC, …) importable when
# this script is launched by a path from elsewhere, e.g.
#   cd EL_Jexpresso && python ../tools/EL_training/train_nn.py
# Python only puts the SCRIPT's folder on sys.path, not the run directory, so we
# prepend the current working directory (where train_common_EL.py lives).
sys.path.insert(0, os.getcwd())

import torch
import torch.nn as nn
import torch.onnx
import matplotlib
import matplotlib.pyplot as plt
from torch.utils.data import TensorDataset, DataLoader

from train_common_EL import get_device, setup_problem, train_and_eval

# Skip interactive plotting calls under a non-interactive backend (e.g. the
# MPLBACKEND=Agg the headless pipeline sets), which would otherwise warn
# "FigureCanvasAgg is non-interactive". The figure is still saved to a PNG.
_INTERACTIVE_BACKEND = matplotlib.get_backend().lower() not in {
    'agg', 'pdf', 'ps', 'svg', 'cairo', 'template'
}

# ─────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────
# Tensor filenames / model basename are env-overridable so the pipeline can tag
# them per test/grid (see tools/EL_training/run_element_learning.sh).
inputfile      = os.environ.get('EL_INPUT_TENSOR',  'input_tensor.csv')
outputfile     = os.environ.get('EL_OUTPUT_TENSOR', 'output_tensor.csv')
dataname       = os.environ.get('EL_DATANAME',      'JX_NN')

num_epochs          = 2000
b_size              = 256
ln_rate             = 1e-3          # bumped from 1e-4: standardised data tolerates a larger LR
grad_clip_norm      = 0.5
early_stop_patience = 300

standardize         = True          # ← the fix; set False for the old behaviour

modelname      = dataname + '_model.pth'
best_modelname = dataname + '_best_model.pth'
training_error = dataname + '_error.mat'
onnx_name      = dataname + '_model.onnx'

# ─────────────────────────────────────────────
# Setup (data + model + loss), exactly as before
# ─────────────────────────────────────────────
device = get_device()
dataloader_train, dataloader_test, model, criterion = \
    setup_problem(inputfile, outputfile, b_size, device)

# ─────────────────────────────────────────────
# Standardise inputs AND outputs (train-split statistics)
# ─────────────────────────────────────────────
def _collect(loader):
    xs, ys = [], []
    for xb, yb in loader:
        xs.append(xb); ys.append(yb)
    return torch.cat(xs, 0), torch.cat(ys, 0)

Xtr_raw, Ytr_raw = _collect(dataloader_train)
Xte_raw, Yte_raw = _collect(dataloader_test)

if standardize:
    x_mean = Xtr_raw.mean(0, keepdim=True)
    x_std  = Xtr_raw.std(0, unbiased=False, keepdim=True).clamp_min(1e-8)
    y_mean = Ytr_raw.mean(0, keepdim=True)
    y_std  = Ytr_raw.std(0, unbiased=False, keepdim=True).clamp_min(1e-8)
    print(f"Standardise: inputs {tuple(Xtr_raw.shape)} and outputs {tuple(Ytr_raw.shape)} "
          f"(zero-mean/unit-var, baked into ONNX)")
else:
    x_mean = torch.zeros(1, Xtr_raw.shape[1]); x_std = torch.ones(1, Xtr_raw.shape[1])
    y_mean = torch.zeros(1, Ytr_raw.shape[1]); y_std = torch.ones(1, Ytr_raw.shape[1])
    print("Standardise: disabled (raw data)")

dataloader_train = DataLoader(TensorDataset((Xtr_raw - x_mean) / x_std,
                                            (Ytr_raw - y_mean) / y_std),
                              batch_size=b_size, shuffle=True)
dataloader_test  = DataLoader(TensorDataset((Xte_raw - x_mean) / x_std,
                                            (Yte_raw - y_mean) / y_std),
                              batch_size=b_size, shuffle=False)

# ── Weight initialisation (Kaiming for ReLU networks) ───────────────────────
def init_weights(m):
    if isinstance(m, torch.nn.Linear):
        torch.nn.init.kaiming_normal_(m.weight, nonlinearity='relu')
        torch.nn.init.zeros_(m.bias)
model.apply(init_weights)

# ── Optimiser + LR scheduler ─────────────────────────────────────────────────
optimizer = torch.optim.Adam(model.parameters(), lr=ln_rate)
scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, mode='min', factor=0.5, patience=50, min_lr=1e-8, threshold=1e-6)

# ─────────────────────────────────────────────
# Real-time plot
# ─────────────────────────────────────────────
if _INTERACTIVE_BACKEND:
    plt.ion()
fig, ax = plt.subplots(figsize=(10, 6))
line_train, = ax.plot([], [], 'b-', label='Train Loss', linewidth=1.5)
line_test,  = ax.plot([], [], 'r-', label='Test Loss',  linewidth=1.5)
ax.set_yscale('log'); ax.set_xlabel('Epochs'); ax.set_ylabel('Loss (Log Scale)')
ax.set_title(f'Real-time Training Monitor: {dataname}')
ax.legend(); ax.grid(True, which="both", ls="-", alpha=0.3)
plt.tight_layout()

def update_plot(epoch_idx, train_hist, test_hist):
    if not _INTERACTIVE_BACKEND:
        return
    if epoch_idx % 10 == 0:
        r = list(range(len(train_hist)))
        line_train.set_data(r, train_hist); line_test.set_data(r, test_hist)
        ax.relim(); ax.autoscale_view()
        fig.canvas.draw(); fig.canvas.flush_events(); plt.pause(0.001)

# ─────────────────────────────────────────────
# Training (unchanged loop)
# ─────────────────────────────────────────────
print(f"Starting training on {device} for up to {num_epochs} epochs...")
model, optimizer, stopped_epoch, train_loss_history, test_loss_history = \
    train_and_eval(
        model               = model,
        optimizer           = optimizer,
        dataloader_train    = dataloader_train,
        dataloader_test     = dataloader_test,
        criterion           = criterion,
        num_epochs          = num_epochs,
        modelname           = modelname,
        best_modelname      = best_modelname,
        training_error      = training_error,
        scheduler           = scheduler,
        early_stop_patience = early_stop_patience,
        grad_clip_norm      = grad_clip_norm,
        plot_callback       = update_plot,
    )

# ─────────────────────────────────────────────
# Final plot
# ─────────────────────────────────────────────
r = list(range(len(train_loss_history)))
line_train.set_data(r, train_loss_history); line_test.set_data(r, test_loss_history)
ax.relim(); ax.autoscale_view()
ax.text(0.98, 0.97,
        (f"num_epochs   = {num_epochs}\nb_size       = {b_size}\n"
         f"ln_rate      = {ln_rate:.1e}\nstandardize  = {standardize}"),
        transform=ax.transAxes, fontsize=9, va='top', ha='right',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.8, edgecolor='gray'),
        fontfamily='monospace')
fig.canvas.draw()
if _INTERACTIVE_BACKEND:
    fig.canvas.flush_events()

# Always persist the loss-history figure so it is captured even headless.
plot_name = dataname + '_training.png'
fig.savefig(plot_name, dpi=150)
print(f"Saved training plot : {plot_name}")

# ─────────────────────────────────────────────
# ONNX export — bake standardisation IN/OUT so Julia feeds raw â, reads raw T^ie
# ─────────────────────────────────────────────
class _ExportWrapper(nn.Module):
    """raw x → (x-x_mean)/x_std → base → *y_std + y_mean → raw output."""
    def __init__(self, base, x_mean, x_std, y_mean, y_std):
        super().__init__()
        self.base = base
        self.register_buffer('x_mean', x_mean); self.register_buffer('x_std', x_std)
        self.register_buffer('y_mean', y_mean); self.register_buffer('y_std', y_std)
    def forward(self, x):
        return self.base((x - self.x_mean) / self.x_std) * self.y_std + self.y_mean

cpu = torch.device("cpu")
best_ckpt = torch.load(best_modelname, map_location=cpu)
model_cpu = model.to(cpu)
model_cpu.load_state_dict(best_ckpt['model_state_dict'])
model_cpu.eval()

export_model = _ExportWrapper(model_cpu,
                              x_mean.to(cpu), x_std.to(cpu),
                              y_mean.to(cpu), y_std.to(cpu)).to(cpu)
export_model.eval()

dummy_input = Xtr_raw[:1].to(cpu)          # RAW feature row (wrapper standardises internally)
torch.onnx.export(
    export_model, dummy_input, onnx_name,
    export_params=True, opset_version=17, do_constant_folding=True,
    input_names=['input'], output_names=['output'],
    dynamic_axes={'input': {0: 'batch_size'}, 'output': {0: 'batch_size'}})
print(f"Saved ONNX model : {onnx_name}  "
      f"(standardisation {'baked in' if standardize else 'disabled'})")

if _INTERACTIVE_BACKEND:
    plt.ioff(); plt.show()
print("Done.")
