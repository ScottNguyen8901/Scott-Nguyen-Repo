import json
import random
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import TensorDataset, DataLoader

DATA_NPZ_PATH = "ml_data_all.npz"

N_EPOCHS = 200
BATCH_SIZE = 128
LR = 3e-4
PATIENCE = 20

DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def load_ml_data(npz_path=DATA_NPZ_PATH):
    data = np.load(npz_path)
    X_norm = data["X_norm"]
    y = data["y"]
    dep_j = data["dep_julian"]
    dep_year = data["dep_year"]

    idx = np.argsort(dep_j)
    X_norm = X_norm[idx]
    y = y[idx]
    dep_j = dep_j[idx]
    dep_year = dep_year[idx]

    return X_norm, y, dep_j, dep_year


def train_val_test_split(X, y, dep_j, train_frac=0.7, val_frac=0.15):
    N = len(dep_j)
    n_train = int(train_frac * N)
    n_val = int(val_frac * N)

    X_train = X[:n_train]
    y_train = y[:n_train]

    X_val = X[n_train:n_train + n_val]
    y_val = y[n_train:n_train + n_val]

    X_test = X[n_train + n_val:]
    y_test = y[n_train + n_val:]

    dep_train = dep_j[:n_train]
    dep_val = dep_j[n_train:n_train + n_val]
    dep_test = dep_j[n_train + n_val:]

    print("Split sizes:")
    print(f"  Train: {X_train.shape[0]}")
    print(f"  Val:   {X_val.shape[0]}")
    print(f"  Test:  {X_test.shape[0]}")

    return (
        X_train, y_train, dep_train,
        X_val, y_val, dep_val,
        X_test, y_test, dep_test,
    )


def make_loaders(X_train, y_train, X_val, y_val, batch_size=BATCH_SIZE):
    X_train_t = torch.from_numpy(X_train).float()
    y_train_t = torch.from_numpy(y_train).float()

    X_val_t = torch.from_numpy(X_val).float()
    y_val_t = torch.from_numpy(y_val).float()

    train_ds = TensorDataset(X_train_t, y_train_t)
    val_ds = TensorDataset(X_val_t, y_val_t)

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size, shuffle=False)

    return train_loader, val_loader


class DVSurrogate(nn.Module):
    def __init__(self, in_dim=2, hidden_dim=128, dropout_p=0.1):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_p),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Dropout(dropout_p),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, x):
        return self.net(x)


def mae(y_true, y_pred):
    return torch.mean(torch.abs(y_true - y_pred))


def rmse(y_true, y_pred):
    return torch.sqrt(torch.mean((y_true - y_pred) ** 2))


def r2_score(y_true, y_pred):
    y_mean = torch.mean(y_true)
    ss_tot = torch.sum((y_true - y_mean) ** 2)
    ss_res = torch.sum((y_true - y_pred) ** 2)
    return 1.0 - ss_res / ss_tot


def train_model(model, train_loader, val_loader, n_epochs=N_EPOCHS, lr=LR, device=DEVICE, patience=PATIENCE):
    model.to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=1e-4)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=5, min_lr=1e-6
    )

    best_val_loss = float("inf")
    best_state_dict = None
    epochs_no_improve = 0

    for epoch in range(1, n_epochs + 1):
        model.train()
        train_loss = 0.0
        for Xb, yb in train_loader:
            Xb = Xb.to(device)
            yb = yb.to(device)

            optimizer.zero_grad()
            preds = model(Xb)
            loss = criterion(preds, yb)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()
            train_loss += loss.item() * Xb.size(0)

        train_loss /= len(train_loader.dataset)

        model.eval()
        val_loss = 0.0
        val_mae = 0.0
        val_rmse = 0.0
        val_r2 = 0.0
        n_val = 0

        with torch.no_grad():
            for Xb, yb in val_loader:
                Xb = Xb.to(device)
                yb = yb.to(device)

                preds = model(Xb)
                loss = criterion(preds, yb)

                bs = Xb.size(0)
                val_loss += loss.item() * bs
                val_mae += mae(yb, preds).item() * bs
                val_rmse += rmse(yb, preds).item() * bs
                val_r2 += r2_score(yb, preds).item() * bs
                n_val += bs

        val_loss /= n_val
        val_mae /= n_val
        val_rmse /= n_val
        val_r2 /= n_val

        scheduler.step(val_loss)

        if epoch % 5 == 0 or epoch == 1:
            current_lr = optimizer.param_groups[0]["lr"]
            print(
                f"Epoch {epoch:03d}: "
                f"train MSE={train_loss:.4f}, "
                f"val MSE={val_loss:.4f}, "
                f"val MAE={val_mae:.4f}, "
                f"val RMSE={val_rmse:.4f}, "
                f"val R2={val_r2:.4f}, "
                f"LR={current_lr:.2e}"
            )

        if val_loss < best_val_loss - 1e-4:
            best_val_loss = val_loss
            best_state_dict = model.state_dict()
            epochs_no_improve = 0
        else:
            epochs_no_improve += 1
            if epochs_no_improve >= patience:
                print(f"Early stopping at epoch {epoch}")
                break

    if best_state_dict is not None:
        model.load_state_dict(best_state_dict)

    return model


def evaluate_on_test(model, X_test, y_test, label="Test", device=DEVICE, y_mean=None, y_std=None):
    model.to(device)
    model.eval()

    X_t = torch.from_numpy(X_test).float().to(device)
    y_t = torch.from_numpy(y_test).float().to(device)

    with torch.no_grad():
        preds = model(X_t)

    if (y_mean is not None) and (y_std is not None):
        preds_denorm = preds * y_std + y_mean
        y_denorm = y_t * y_std + y_mean
    else:
        preds_denorm = preds
        y_denorm = y_t

    test_mae = mae(y_denorm, preds_denorm).item()
    test_rmse = rmse(y_denorm, preds_denorm).item()
    test_r2 = r2_score(y_denorm, preds_denorm).item()

    print(f"{label} performance:")
    print(f"  MAE  = {test_mae:.4f} km/s")
    print(f"  RMSE = {test_rmse:.4f} km/s")
    print(f"  R2   = {test_r2:.4f}")

    return preds_denorm.cpu().numpy()


def train_from_npz(npz_path=DATA_NPZ_PATH, save_path="dv_surrogate_model_2030_2032_supervised.pt"):
    print(f"Using device: {DEVICE}")

    torch.manual_seed(42)
    np.random.seed(42)
    random.seed(42)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(42)

    X_norm, y, dep_j, _ = load_ml_data(npz_path)
    print(f"Total samples: {X_norm.shape[0]}")

    y_mean = y.mean()
    y_std = y.std()
    y_norm = (y - y_mean) / y_std

    (
        X_train, y_train, _,
        X_val, y_val, _,
        X_test, y_test, _,
    ) = train_val_test_split(X_norm, y_norm, dep_j, train_frac=0.7, val_frac=0.15)

    train_loader, val_loader = make_loaders(X_train, y_train, X_val, y_val, batch_size=BATCH_SIZE)

    in_dim = X_norm.shape[1]
    model = DVSurrogate(in_dim=in_dim, hidden_dim=128, dropout_p=0.1)

    model = train_model(model, train_loader, val_loader, n_epochs=N_EPOCHS, lr=LR, device=DEVICE, patience=PATIENCE)

    evaluate_on_test(model, X_val, y_val, label="Validation", device=DEVICE, y_mean=y_mean, y_std=y_std)
    evaluate_on_test(model, X_test, y_test, label="Test", device=DEVICE, y_mean=y_mean, y_std=y_std)
    evaluate_on_test(model, X_norm, y_norm, label="All", device=DEVICE, y_mean=y_mean, y_std=y_std)

    torch.save(
        {
            "model_state_dict": model.state_dict(),
            "y_mean": float(y_mean),
            "y_std": float(y_std),
            "in_dim": in_dim,
            "hidden_dim": 128,
            "dropout_p": 0.1,
        },
        save_path,
    )
    print(f"Saved model: {save_path}")
    return model


if __name__ == "__main__":
    train_from_npz(DATA_NPZ_PATH)
