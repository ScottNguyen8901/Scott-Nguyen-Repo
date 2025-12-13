# Here we import the necessary libraries and the graph plotting function:
import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel, RBF, ExpSineSquared
import os

# Create folder "plots" if it doesn't exist
os.makedirs("plots", exist_ok=True)

def plot_and_save(y, mean, std, title, filename):
    """Saves the GP prediction figure into /plots WITHOUT displaying it."""
    plt.figure()
    plt.plot(y, label="Observations", color="black")
    plt.plot(mean, label="Mean prediction", color="red")
    plt.fill_between(
        np.arange(y.size).ravel(),
        mean - 1.96 * std,
        mean + 1.96 * std,
        alpha=0.5,
        label=r"95% confidence interval",
    )
    plt.legend()
    plt.xlabel("$t$")
    plt.ylabel("$f(t)$")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(f"plots/{filename}", dpi=300)
    plt.close()   # <-- closes figure, prevents display


# === Generate sine wave and training data exactly as before ===

rng = np.random.RandomState(0)

# 1000 samples between 0 and 30 seconds
x = np.linspace(0, 30, num=1_000).reshape(-1, 1)
y = np.sin(x).ravel()

# Save true signal
plt.figure()
plt.plot(x, y)
plt.ylabel('y (regressor)')
plt.xlabel('x (time)')
plt.title("True sine wave")
plt.tight_layout()
plt.savefig("plots/true_signal.png", dpi=300)
plt.close()


# Select 40 time instants at random from the first 400
training_sample_indices = rng.choice(np.arange(0, 400), size=40, replace=False)
x_train = x[training_sample_indices]
y_train = y[training_sample_indices] + 0.1 * np.random.randn(40)

# Save noisy training sample plot
plt.figure()
plt.plot(x, y)
plt.scatter(x_train, y_train, color="red")
plt.ylabel('y (regressor)')
plt.xlabel('x (time)')
plt.title("Training samples (noisy)")
plt.tight_layout()
plt.savefig("plots/training_samples.png", dpi=300)
plt.close()


# ---------------- 8.1(a) Train + Test Each Kernel ----------------

# ----- Linear kernel: DotProduct + WhiteKernel -----
kernel_lin = 1 * DotProduct(sigma_0=1.0, sigma_0_bounds=(1e-1, 10.0)) + WhiteKernel(0.1)

gp_lin = GaussianProcessRegressor(
    kernel=kernel_lin,
    alpha=0.0,
    normalize_y=True,
    n_restarts_optimizer=10,
    random_state=0,
)

gp_lin.fit(x_train, y_train)
y_mean_lin, y_std_lin = gp_lin.predict(x, return_std=True)

print("Learned linear kernel:", gp_lin.kernel_)

plot_and_save(
    y, y_mean_lin, y_std_lin,
    "Linear kernel (DotProduct + WhiteKernel)",
    "gp_linear.png"
)


# ----- Squared Exponential kernel: RBF + WhiteKernel -----
kernel_rbf = 1 * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 10.0)) + WhiteKernel(0.1)

gp_rbf = GaussianProcessRegressor(
    kernel=kernel_rbf,
    alpha=0.0,
    normalize_y=True,
    n_restarts_optimizer=10,
    random_state=0,
)

gp_rbf.fit(x_train, y_train)
y_mean_rbf, y_std_rbf = gp_rbf.predict(x, return_std=True)

print("Learned RBF kernel:", gp_rbf.kernel_)

plot_and_save(
    y, y_mean_rbf, y_std_rbf,
    "Squared Exponential (RBF + WhiteKernel)",
    "gp_rbf.png"
)


# ----- Exp-Sine-Squared (periodic) kernel -----
# Give reasonable bounds and more restarts so the optimizer can find
# a good period (~ 2*pi) and length-scale.
kernel_per = 1 * ExpSineSquared(
    length_scale=1.0,
    length_scale_bounds=(1e-2, 10.0),
    periodicity=2 * np.pi,
    periodicity_bounds=(1.0, 10.0),
)

gp_per = GaussianProcessRegressor(
    kernel=kernel_per,
    alpha=0.0,
    normalize_y=True,
    n_restarts_optimizer=50,  # <<-- much larger to avoid bad local minima
    random_state=0,
)

gp_per.fit(x_train, y_train)
y_mean_per, y_std_per = gp_per.predict(x, return_std=True)

print("Learned periodic kernel:", gp_per.kernel_)

plot_and_save(
    y, y_mean_per, y_std_per,
    "Exp-Sine-Squared (periodic) kernel",
    "gp_periodic.png"
)

# ------------------------
# 8.1(c) Inspect kernel parameters
# ------------------------

print("\n========== Learned Kernel Parameters ==========\n")

print("Linear kernel learned parameters:")
print(gp_lin.kernel_)
print()

print("RBF kernel learned parameters:")
print(gp_rbf.kernel_)
print()

print("Periodic kernel learned parameters:")
print(gp_per.kernel_)
print()

# ------------------------
# 8.1(d) Periodic kernel WITHOUT WhiteKernel
# ------------------------

kernel_per_no_white = ExpSineSquared(
    length_scale=1.0,
    length_scale_bounds=(1e-2, 10.0),
    periodicity=2*np.pi,
    periodicity_bounds=(1.0, 10.0),
)

gp_per_no_white = GaussianProcessRegressor(
    kernel=kernel_per_no_white,
    alpha=0.0,
    normalize_y=True,
    n_restarts_optimizer=50,
    random_state=0,
)

gp_per_no_white.fit(x_train, y_train)
y_mean_per_no_white, y_std_per_no_white = gp_per_no_white.predict(x, return_std=True)

print("Learned periodic kernel (no WhiteKernel):")
print(gp_per_no_white.kernel_)
print()

plot_and_save(
    y, y_mean_per_no_white, y_std_per_no_white,
    "Periodic kernel WITHOUT WhiteKernel",
    "gp_periodic_no_white.png"
)