import torch
import numpy as np
from sklearn.model_selection import train_test_split

def csv2pyt_fc(
    input_file,
    output_file,
    device,
    b_size=100,
    test_split=0.2,
):
    # --- read CSVs, skipping header row with column names ---
    X = np.loadtxt(input_file, delimiter=',', skiprows=1)
    Y = np.loadtxt(output_file, delimiter=',', skiprows=1)

    # Expect shapes: (N_in, N_samp) and (N_out, N_samp)
    N_in,  N_samp  = X.shape
    N_out, N_samp2 = Y.shape
    assert N_samp == N_samp2, "Input/output sample counts mismatch"

    # columns = samples → transpose
    X = X.T   # (N_samp, N_in)
    Y = Y.T   # (N_samp, N_out)

    X_train, X_test, Y_train, Y_test = train_test_split(
        X, Y, test_size=test_split, random_state=42
    )

    X_train = torch.from_numpy(X_train).float().to(device)
    Y_train = torch.from_numpy(Y_train).float().to(device)
    X_test  = torch.from_numpy(X_test ).float().to(device)
    Y_test  = torch.from_numpy(Y_test ).float().to(device)

    train_ds = torch.utils.data.TensorDataset(X_train, Y_train)
    test_ds  = torch.utils.data.TensorDataset(X_test,  Y_test)

    train_loader = torch.utils.data.DataLoader(
        train_ds, batch_size=b_size, shuffle=True
    )
    test_loader = torch.utils.data.DataLoader(
        test_ds, batch_size=b_size, shuffle=False
    )

    return train_loader, test_loader, N_samp, N_in, N_out
