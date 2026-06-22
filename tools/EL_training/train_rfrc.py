"""
train_rfrc.py  —  Recurrence-Free Reservoir Computing (RF-RC) trainer
======================================================================

RF-RC replaces the iterative backpropagation loop of a classical NN with:

  Phase 1 – ONE-SHOT  : closed-form ridge regression on the frozen reservoir
                        states.  Solves  W = (HᵀH + αI)⁻¹ HᵀY  analytically.
                        No epochs, no gradient computation — extremely fast.

  Phase 2 – FINE-TUNE : optional short gradient-descent pass on the readout
                        layer only, using the Phase-1 solution as a warm start.

The ONNX export at the end is identical to train_gpu.py so the saved model is
a drop-in replacement.

──────────────────────────────────────────────────────────────────────────────
NON-CONSTANT DIFFUSIVITY UPDATE
──────────────────────────────────────────────────────────────────────────────
The element-learning sampler now supports a non-constant diffusivity feature
(Option 1, synthesized Jacobians). Instead of the old per-node SCALAR coefficient
(input size (k+1)^2), each element is described by the per-node 2×2 SPD reference
diffusivity tensor â(ξ) = (|K|/|ref|)·a·J_K⁻¹J_K⁻ᵀ, stored as its 3 unique
entries (â11, â12, â22) at every one of the (k+1)^2 nodes:

        input size  =  3·(k+1)^2          (was (k+1)^2)
        layout      =  [â11_1, â12_1, â22_1,  â11_2, â12_2, â22_2,  ...]
        output size =  len(T^{ie}(:))      (unchanged)

Two changes make this work cleanly:

  1. SIZE  — `N_in` is already read from input_tensor.csv, so the model is sized
     automatically. Nothing to hard-code; the script below just reports the
     detected layout and (optionally) the inferred polynomial order.

  2. SCALE — unlike the old homogeneous scalar feature, the â feature is
     heterogeneous: â11, â22 are positive and O(1), while â12 is signed and on a
     different scale. We standardise each input feature (zero mean, unit variance
     from the TRAIN split) so the reservoir is well conditioned. The
     standardisation statistics are BAKED INTO THE EXPORTED ONNX GRAPH (as a fixed
     affine pre-layer), so the Julia inference side keeps feeding the RAW â vector
     exactly as before — no change needed in Jexpresso.

Set `standardize_inputs = False` to recover the original (un-normalised) behaviour
and a plain ONNX export.
"""

import time
import torch
import torch.nn as nn
import torch.onnx
import matplotlib.pyplot as plt
from scipy.io import savemat
from torch.utils.data import TensorDataset, DataLoader

from IO_EL      import csv2pyt_fc
from NN_RFRC    import RFRC
from SLmodel_EL import save_checkpoint
from train_common_EL import get_device

# ─────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────
inputfile  = 'input_tensor.csv'
outputfile = 'output_tensor.csv'
dataname   = 'JX_RFRC'

# ── Reservoir hyper-parameters ────────────────
reservoir_dim  = 2048   # width of each reservoir layer
num_layers     = 4     # depth of stacked reservoir
activation     = 'tanh' # 'tanh' | 'relu' | 'sigmoid' | 'elu'
input_scaling  = 1.0    # half-range of the uniform weight initialiser
ridge_alpha    = 1e-5   # L2 regularisation for the ridge regression fit

# ── Non-constant diffusivity feature handling ─────────────────────────────
standardize_inputs = True   # per-feature zero-mean/unit-var; baked into ONNX
expect_ahat_feature = False  # if True, assert N_in is divisible by 3 (3 per node)

# ── Optional fine-tuning (set to 0 to disable) ────────────────────────────
finetune_epochs        = 300     # gradient-descent epochs on the readout only
finetune_lr            = 1e-3
finetune_b_size        = 256
grad_clip_norm         = 0.5
early_stop_patience    = 50

# ── Data / output files ───────────────────────
b_size         = 256
modelname      = dataname + '_model.pth'
best_modelname = dataname + '_best_model.pth'
training_error = dataname + '_error.mat'
onnx_name      = dataname + '_model.onnx'

# ─────────────────────────────────────────────
# Device
# ─────────────────────────────────────────────
device = get_device()

# ─────────────────────────────────────────────
# Data
# ─────────────────────────────────────────────
dataloader_train, dataloader_test, N_samp, N_in, N_out = \
    csv2pyt_fc(inputfile, outputfile, device, b_size, test_split=0.2)

print(f"\nDataset  : {N_samp} samples  |  {N_in} inputs  |  {N_out} outputs")

# ── Report the detected feature layout ────────────────────────────────────
if N_in % 3 == 0:
    n_nodes = N_in // 3
    k_plus1 = int(round(n_nodes ** 0.5))
    if k_plus1 * k_plus1 == n_nodes:
        print(f"Feature  : per-node 2×2 SPD â  →  {n_nodes} nodes = (k+1)² with "
              f"k+1={k_plus1} (polynomial order k={k_plus1 - 1}), 3 entries/node")
    else:
        print(f"Feature  : {N_in} inputs, divisible by 3 ({n_nodes} nodes) but "
              f"node count is not a perfect square — check the sampler.")
else:
    print(f"Feature  : {N_in} inputs — NOT divisible by 3 → legacy scalar "
          f"diffusivity feature assumed.")
    if expect_ahat_feature:
        raise ValueError(
            f"expect_ahat_feature=True but N_in={N_in} is not divisible by 3; "
            f"the non-constant â feature must have 3 entries per node.")

# ─────────────────────────────────────────────
# Per-feature input standardisation (train-split statistics)
# ─────────────────────────────────────────────
# Collected once by iterating the loaders (which yield (inputs, targets) — the
# same contract the fine-tune loop already relies on). The RAW tensors are kept
# so the standardisation can be undone/baked into the ONNX export, and so the
# ONNX dummy input is a RAW feature vector.
def _collect(loader):
    xs, ys = [], []
    for xb, yb in loader:
        xs.append(xb)
        ys.append(yb)
    return torch.cat(xs, 0), torch.cat(ys, 0)

X_train_raw, Y_train = _collect(dataloader_train)
X_test_raw,  Y_test  = _collect(dataloader_test)

if standardize_inputs:
    feat_mean = X_train_raw.mean(dim=0, keepdim=True)
    feat_std  = X_train_raw.std(dim=0, unbiased=False, keepdim=True).clamp_min(1e-8)
    print(f"Standardise: per-feature zero-mean/unit-var from {X_train_raw.shape[0]} "
          f"train rows  (baked into ONNX)")
else:
    feat_mean = torch.zeros(1, N_in, device=X_train_raw.device, dtype=X_train_raw.dtype)
    feat_std  = torch.ones( 1, N_in, device=X_train_raw.device, dtype=X_train_raw.dtype)
    print("Standardise: disabled (raw inputs)")

# Rebuild normalised dataloaders. The reservoir (fit + fine-tune + evaluate) now
# always sees standardised inputs; the export wrapper re-applies the same
# transform to RAW inputs so inference is unchanged.
X_train_n = (X_train_raw - feat_mean) / feat_std
X_test_n  = (X_test_raw  - feat_mean) / feat_std

dataloader_train = DataLoader(TensorDataset(X_train_n, Y_train),
                              batch_size=b_size, shuffle=True)
dataloader_test  = DataLoader(TensorDataset(X_test_n,  Y_test),
                              batch_size=b_size, shuffle=False)

# ─────────────────────────────────────────────
# Model
# ─────────────────────────────────────────────
model = RFRC(
    input_size    = N_in,
    output_size   = N_out,
    reservoir_dim = reservoir_dim,
    num_layers    = num_layers,
    activation    = activation,
    input_scaling = input_scaling,
    ridge_alpha   = ridge_alpha,
).to(device)

trainable = sum(p.numel() for p in model.parameters() if p.requires_grad)
total     = sum(p.numel() for p in model.parameters())
print(f"Model    : {total:,} total params  |  {trainable:,} trainable (readout only)")
print(f"Reservoir: dim={reservoir_dim}  layers={num_layers}  "
      f"activation={activation}  alpha={ridge_alpha}\n")

criterion = torch.nn.L1Loss()

# ─────────────────────────────────────────────
# Real-time plot setup
# ─────────────────────────────────────────────
plt.ion()
fig, ax = plt.subplots(figsize=(10, 6))
line_train, = ax.plot([], [], 'b-', label='Train Loss (fine-tune)', linewidth=1.5)
line_test,  = ax.plot([], [], 'r-', label='Test Loss  (fine-tune)', linewidth=1.5)
ax.set_yscale('log')
ax.set_xlabel('Fine-tune Epochs')
ax.set_ylabel('Loss (Log Scale)')
ax.set_title(f'RF-RC Training Monitor: {dataname}')
ax.legend()
ax.grid(True, which="both", ls="-", alpha=0.3)
plt.tight_layout()

# ─────────────────────────────────────────────
# Phase 1 — One-shot ridge regression
# ─────────────────────────────────────────────
print("=" * 60)
print("Phase 1  — One-shot ridge regression fit")
print("=" * 60)

t0   = time.perf_counter()
info = model.fit_readout(dataloader_train, device)
t1   = time.perf_counter()

train_loss_ridge = model.evaluate(dataloader_train, criterion)
test_loss_ridge  = model.evaluate(dataloader_test,  criterion)

print(f"  H matrix : {info['H_shape']}  |  Y matrix : {info['Y_shape']}")
print(f"  W_out    : {info['W_shape']}  |  alpha={info['ridge_alpha']}")
print(f"  Fit time : {t1 - t0:.3f} s")
print(f"  Train L1 : {train_loss_ridge:.10f}")
print(f"  Test  L1 : {test_loss_ridge:.10f}\n")

train_loss_history = []
test_loss_history  = []

# ─────────────────────────────────────────────
# Phase 2 — Optional gradient-descent fine-tune
# ─────────────────────────────────────────────
if finetune_epochs > 0:
    print("=" * 60)
    print(f"Phase 2  — Fine-tuning readout for up to {finetune_epochs} epochs")
    print("=" * 60)

    # Only optimise the readout layer
    optimizer = torch.optim.Adam(model.readout.parameters(), lr=finetune_lr)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min', factor=0.5, patience=20, min_lr=1e-8, threshold=1e-6
    )

    best_test_loss   = test_loss_ridge
    patience_counter = 0

    # Save the ridge solution as the initial best checkpoint
    torch.save({
        'epoch':                0,
        'model_state_dict':     {k: v.cpu() for k, v in model.state_dict().items()},
        'optimizer_state_dict': optimizer.state_dict(),
        'best_test_loss':       best_test_loss,
    }, best_modelname)

    for epoch in range(1, finetune_epochs + 1):
        # ── Train ──────────────────────────────────────────────────────────
        model.train()
        epoch_train_loss = 0.0
        for inputs, targets in dataloader_train:
            optimizer.zero_grad()
            loss = criterion(model(inputs), targets)
            loss.backward()
            if grad_clip_norm is not None:
                torch.nn.utils.clip_grad_norm_(
                    model.readout.parameters(), max_norm=grad_clip_norm)
            optimizer.step()
            epoch_train_loss += loss.item()
        avg_train = epoch_train_loss / len(dataloader_train)

        # ── Evaluate ───────────────────────────────────────────────────────
        avg_test = model.evaluate(dataloader_test, criterion)

        train_loss_history.append(avg_train)
        test_loss_history.append(avg_test)

        scheduler.step(avg_test)
        current_lr = optimizer.param_groups[0]['lr']

        # ── Best checkpoint ────────────────────────────────────────────────
        if avg_test < best_test_loss:
            best_test_loss   = avg_test
            patience_counter = 0
            torch.save({
                'epoch':                epoch,
                'model_state_dict':     {k: v.cpu()
                                         for k, v in model.state_dict().items()},
                'optimizer_state_dict': optimizer.state_dict(),
                'best_test_loss':       best_test_loss,
            }, best_modelname)
        else:
            patience_counter += 1

        print(
            f"  Epoch [{epoch:4d}/{finetune_epochs}] | "
            f"Train: {avg_train:.10f} | Test: {avg_test:.10f} | "
            f"LR: {current_lr:.2e} | "
            f"Patience: {patience_counter}/{early_stop_patience}"
        )

        # ── Live plot ──────────────────────────────────────────────────────
        if epoch % 10 == 0:
            epochs_range = list(range(len(train_loss_history)))
            line_train.set_data(epochs_range, train_loss_history)
            line_test.set_data(epochs_range,  test_loss_history)
            ax.relim()
            ax.autoscale_view()
            fig.canvas.draw()
            fig.canvas.flush_events()
            plt.pause(0.001)

        # ── Early stopping ─────────────────────────────────────────────────
        if patience_counter >= early_stop_patience:
            print(f"\n  Early stopping at epoch {epoch} "
                  f"(no improvement for {early_stop_patience} epochs).")
            break

    print(f"\n  Best test loss after fine-tuning: {best_test_loss:.10f}")

else:
    # No fine-tuning — treat the ridge solution as the best model
    optimizer = torch.optim.SGD(model.readout.parameters(), lr=0.0)
    torch.save({
        'epoch':                0,
        'model_state_dict':     {k: v.cpu() for k, v in model.state_dict().items()},
        'optimizer_state_dict': optimizer.state_dict(),
        'best_test_loss':       test_loss_ridge,
    }, best_modelname)
    print("Fine-tuning disabled  (finetune_epochs=0).")

# ─────────────────────────────────────────────
# Save last checkpoint + loss history
# ─────────────────────────────────────────────
save_checkpoint(model, optimizer,
                len(train_loss_history),
                train_loss_history, test_loss_history,
                modelname)
print(f"\nSaved last checkpoint : {modelname}")
print(f"Saved best checkpoint : {best_modelname}")

savemat(training_error, {
    'train_err': train_loss_history,
    'test_err':  test_loss_history,
    'ridge_train_loss': train_loss_ridge,
    'ridge_test_loss':  test_loss_ridge,
})

# ─────────────────────────────────────────────
# Final plot — info box
# ─────────────────────────────────────────────
if train_loss_history:
    epochs_range = list(range(len(train_loss_history)))
    line_train.set_data(epochs_range, train_loss_history)
    line_test.set_data(epochs_range,  test_loss_history)
    ax.relim()
    ax.autoscale_view()

hparam_text = (
    f"input_size          = {N_in}\n"
    f"standardize_inputs  = {standardize_inputs}\n"
    f"reservoir_dim       = {reservoir_dim}\n"
    f"num_layers          = {num_layers}\n"
    f"activation          = {activation}\n"
    f"ridge_alpha         = {ridge_alpha:.1e}\n"
    f"finetune_epochs     = {finetune_epochs}\n"
    f"finetune_lr         = {finetune_lr:.1e}\n"
    f"grad_clip_norm      = {grad_clip_norm}\n"
    f"early_stop_patience = {early_stop_patience}\n"
    f"─────────────────────────────\n"
    f"ridge train L1      = {train_loss_ridge:.6f}\n"
    f"ridge test  L1      = {test_loss_ridge:.6f}"
)
ax.text(
    0.98, 0.97,
    hparam_text,
    transform           = ax.transAxes,
    fontsize            = 8,
    verticalalignment   = 'top',
    horizontalalignment = 'right',
    bbox                = dict(boxstyle='round,pad=0.4',
                               facecolor='white',
                               alpha=0.8,
                               edgecolor='gray'),
    fontfamily          = 'monospace',
)
fig.canvas.draw()
fig.canvas.flush_events()

# ─────────────────────────────────────────────
# ONNX export — from best weights, on CPU
# ─────────────────────────────────────────────
# The exported graph consumes the RAW â feature vector: a fixed affine pre-layer
# (x - feat_mean)/feat_std reproduces the training-time standardisation inside the
# ONNX, so the Julia inference side feeds raw features exactly as before. When
# standardize_inputs=False the pre-layer is the identity (mean=0, std=1).
class _ExportWrapper(nn.Module):
    """Bakes the input standardisation into the model for inference/ONNX."""
    def __init__(self, base, feat_mean, feat_std):
        super().__init__()
        self.base = base
        self.register_buffer('feat_mean', feat_mean)
        self.register_buffer('feat_std',  feat_std)

    def forward(self, x):
        return self.base((x - self.feat_mean) / self.feat_std)


cpu_device = torch.device("cpu")
best_ckpt  = torch.load(best_modelname, map_location=cpu_device)
model_cpu  = model.to(cpu_device)
model_cpu.load_state_dict(best_ckpt['model_state_dict'])
model_cpu.eval()

export_model = _ExportWrapper(
    model_cpu,
    feat_mean.to(cpu_device),
    feat_std.to(cpu_device),
).to(cpu_device)
export_model.eval()

# Dummy input is a RAW feature row (NOT standardised), since the wrapper
# standardises internally.
dummy_input = X_train_raw[:1].to(cpu_device)
torch.onnx.export(
    export_model,
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
print(f"Saved ONNX model : {onnx_name}  (exported from best checkpoint, "
      f"standardisation {'baked in' if standardize_inputs else 'disabled'})")

# ─────────────────────────────────────────────
# Finalise
# ─────────────────────────────────────────────
plt.ioff()
plt.show()
print("Done.")
