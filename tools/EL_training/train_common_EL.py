# train_common_EL.py
import os
import torch
from NN_EL      import FCNN
from IO_EL      import csv2pyt_fc
from SLmodel_EL import save_checkpoint
from scipy.io   import savemat

# ── CPU thread count (M2 Air has 8 performance cores) ───────────────────────
_NUM_THREADS = 8
torch.set_num_threads(_NUM_THREADS)
os.environ["OMP_NUM_THREADS"] = str(_NUM_THREADS)


def get_device():
    """
    Select compute device.
    For networks with ~42k parameters MPS kernel-launch overhead dominates,
    so CPU is hardcoded as the optimal choice on Apple Silicon for this workload.
    Change FORCE_CPU = False to re-enable automatic GPU selection.
    """
    FORCE_CPU = True

    if not FORCE_CPU:
        if torch.cuda.is_available():
            device = torch.device("cuda:0")
        elif torch.backends.mps.is_available():
            device = torch.device("mps")
        else:
            device = torch.device("cpu")
    else:
        device = torch.device("cpu")

    print(f"Training on : {device}  "
          f"(threads={torch.get_num_threads()}, FORCE_CPU={FORCE_CPU})")
    return device


def setup_problem(inputfile, outputfile, b_size, device):
    """
    Create dataloaders, model, and loss function.
    Returns:
        dataloader_train, dataloader_test, model, criterion
    """
    dataloader_train, dataloader_test, N_samp, N_in, N_out = \
        csv2pyt_fc(inputfile, outputfile, device, b_size, test_split=0.2)
    model     = FCNN(input_size=N_in, output_size=N_out).to(device)
    criterion = torch.nn.L1Loss()
    return dataloader_train, dataloader_test, model, criterion


def train_and_eval(
    model,
    optimizer,
    dataloader_train,
    dataloader_test,
    criterion,
    num_epochs,
    start_epoch          = 0,
    train_loss_history   = None,
    test_loss_history    = None,
    modelname            = None,
    best_modelname       = None,
    training_error       = None,
    # ── new keyword arguments ───────────────────────────────────────────────
    scheduler            = None,   # e.g. ReduceLROnPlateau instance
    early_stop_patience  = None,   # int; None = disabled
    grad_clip_norm       = 1.0,    # max_norm for clip_grad_norm_; None = disabled
    plot_callback        = None,   # callable(epoch_idx, train_hist, test_hist)
):
    """
    Core training + evaluation loop with optional:
      - LR scheduling (ReduceLROnPlateau or any torch scheduler)
      - Early stopping
      - Gradient clipping
      - Best-model checkpoint (best_modelname)

    Returns:
        model, optimizer, epoch_idx,
        train_loss_history, test_loss_history
    """
    if train_loss_history is None:
        train_loss_history = []
    if test_loss_history is None:
        test_loss_history = []

    total_epochs     = start_epoch + num_epochs
    epoch_idx        = start_epoch
    best_test_loss   = float('inf')
    patience_counter = 0
    stopped_early    = False

    for _ in range(num_epochs):
        epoch_idx += 1

        # ── Train ────────────────────────────────────────────────────────────
        model.train()
        epoch_train_loss = 0.0
        for inputs, targets in dataloader_train:
            optimizer.zero_grad()
            outputs = model(inputs)
            loss    = criterion(outputs, targets)
            loss.backward()
            if grad_clip_norm is not None:
                torch.nn.utils.clip_grad_norm_(
                    model.parameters(), max_norm=grad_clip_norm)
            optimizer.step()
            epoch_train_loss += loss.item()
        avg_train_loss = epoch_train_loss / len(dataloader_train)
        train_loss_history.append(avg_train_loss)

        # ── Evaluate ─────────────────────────────────────────────────────────
        model.eval()
        epoch_test_loss = 0.0
        with torch.no_grad():
            for inputs, targets in dataloader_test:
                outputs = model(inputs)
                epoch_test_loss += criterion(outputs, targets).item()
        avg_test_loss = epoch_test_loss / len(dataloader_test)
        test_loss_history.append(avg_test_loss)

        # ── LR scheduler ─────────────────────────────────────────────────────
        if scheduler is not None:
            if isinstance(scheduler,
                          torch.optim.lr_scheduler.ReduceLROnPlateau):
                scheduler.step(avg_test_loss)
            else:
                scheduler.step()

        # ── Best-model checkpoint ─────────────────────────────────────────────
        if avg_test_loss < best_test_loss:
            best_test_loss   = avg_test_loss
            patience_counter = 0
            if best_modelname is not None:
                cpu_state = {k: v.cpu()
                             for k, v in model.state_dict().items()}
                torch.save({
                    'epoch':                epoch_idx,
                    'model_state_dict':     cpu_state,
                    'optimizer_state_dict': optimizer.state_dict(),
                    'best_test_loss':       best_test_loss,
                }, best_modelname)
        else:
            patience_counter += 1

        # ── Console output ───────────────────────────────────────────────────
        current_lr = optimizer.param_groups[0]['lr']
        print(
            f"Epoch [{epoch_idx}/{total_epochs}] | "
            f"Train Loss: {avg_train_loss:.10f} | "
            f"Test Loss: {avg_test_loss:.10f} | "
            f"LR: {current_lr:.2e}" +
            (f" | Patience: {patience_counter}/{early_stop_patience}"
             if early_stop_patience is not None else "")
        )

        # ── Plot callback ─────────────────────────────────────────────────────
        if plot_callback is not None:
            plot_callback(epoch_idx, train_loss_history, test_loss_history)

        # ── Early stopping ────────────────────────────────────────────────────
        if early_stop_patience is not None and \
                patience_counter >= early_stop_patience:
            print(f"\nEarly stopping at epoch {epoch_idx} "
                  f"(no improvement for {early_stop_patience} epochs).")
            stopped_early = True
            break

    # ── Save last-epoch checkpoint ────────────────────────────────────────────
    if modelname is not None:
        cpu_state = {k: v.cpu() for k, v in model.state_dict().items()}
        save_checkpoint(model, optimizer, epoch_idx,
                        train_loss_history, test_loss_history, modelname)
        if best_modelname is not None:
            tag = "best" if stopped_early else "last==best"
            print(f"Saved last checkpoint : {modelname}")
            print(f"Saved best checkpoint : {best_modelname} "
                  f"(test loss = {best_test_loss:.10f}, {tag})")

    # ── Save loss history ─────────────────────────────────────────────────────
    if training_error is not None:
        savemat(training_error, {
            'train_err': train_loss_history,
            'test_err':  test_loss_history,
        })

    return model, optimizer, epoch_idx, train_loss_history, test_loss_history
