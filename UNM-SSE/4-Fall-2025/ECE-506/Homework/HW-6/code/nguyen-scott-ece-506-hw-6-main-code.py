import torch
import torch.optim as optim
import torch.nn as nn
import matplotlib.pyplot as plt
import numpy as np

def get_loss(c_pred, c_true):
    return torch.mean((c_pred - c_true)**2)

def predict(alpha0, alpha1, f):
    return alpha0 + alpha1 * f

# ============================================================
# Problem 3 – Automatic Differentiation
# ============================================================

x = torch.ones(6, requires_grad=True)
f = (1/6) * torch.sum((x - 2)**2) + torch.sum(x)
f.backward()
print("Problem 3 gradient:", x.grad)

# ============================================================
# Problem 4(a) – Fahrenheit → Celsius using Autograd
# ============================================================

f_vals = 100 * torch.rand(10) + 32
c_true = (f_vals - 32) / 1.8
f_mean = f_vals.mean()
f_std = f_vals.std()
f_norm = (f_vals - f_mean) / f_std

beta0 = torch.zeros(1, requires_grad=True)
beta1 = torch.zeros(1, requires_grad=True)
optimizer = optim.SGD([beta0, beta1], lr=0.1)

max_iters = 2000
loss = None

for iteration in range(max_iters):
    optimizer.zero_grad()
    c_pred_norm = beta0 + beta1 * f_norm
    loss = get_loss(c_pred_norm, c_true)
    loss.backward()
    optimizer.step()

with torch.no_grad():
    alpha1 = beta1 / f_std
    alpha0 = beta0 - beta1 * f_mean / f_std

true_formula = "True:    C = -17.78 + 0.56·F"
learned_formula = f"Learned: C = {alpha0.item():.2f} + {alpha1.item():.2f}·F"

print("\nFahrenheit to Celsius")
print("\t" + true_formula)
print("\t" + learned_formula)
print("Final loss:", loss.item())
print("Iterations:", max_iters)

# ============================================================
# Problem 4(b) – Neural-Net Framework with SGD and Momentum
# ============================================================

class LinearNet(nn.Module):
    def __init__(self):
        super(LinearNet, self).__init__()
        self.linear = nn.Linear(1, 1)

    def forward(self, x):
        return self.linear(x)

f_input_norm = f_norm.unsqueeze(1)
c_target = c_true.unsqueeze(1)
criterion = nn.MSELoss()

net_sgd = LinearNet()
optimizer_sgd = optim.SGD(net_sgd.parameters(), lr=0.1, momentum=0.0)
num_epochs = 2000

for epoch in range(num_epochs):
    optimizer_sgd.zero_grad()
    c_pred = net_sgd(f_input_norm)
    loss_sgd = criterion(c_pred, c_target)
    loss_sgd.backward()
    optimizer_sgd.step()

with torch.no_grad():
    w_sgd = net_sgd.linear.weight.item()
    b_sgd = net_sgd.linear.bias.item()
    alpha1_sgd = w_sgd / f_std
    alpha0_sgd = b_sgd - w_sgd * f_mean / f_std

net_mom = LinearNet()
optimizer_mom = optim.SGD(net_mom.parameters(), lr=0.1, momentum=0.9)

for epoch in range(num_epochs):
    optimizer_mom.zero_grad()
    c_pred_m = net_mom(f_input_norm)
    loss_mom = criterion(c_pred_m, c_target)
    loss_mom.backward()
    optimizer_mom.step()

with torch.no_grad():
    w_mom = net_mom.linear.weight.item()
    b_mom = net_mom.linear.bias.item()
    alpha1_mom = w_mom / f_std
    alpha0_mom = b_mom - w_mom * f_mean / f_std

print("\nProblem 4(b): Neural network with SGD")
print("SGD (no momentum)")
print(f"\tLearned: C = {alpha0_sgd.item():.2f} + {alpha1_sgd.item():.2f}·F")
print(f"\tFinal loss: {loss_sgd.item()}")

print("\nSGD (with momentum = 0.9)")
print(f"\tLearned: C = {alpha0_mom.item():.2f} + {alpha1_mom.item():.2f}·F")
print(f"\tFinal loss: {loss_mom.item()}")

# ============================================================
# Problem 5(e) – Q1, Q2, Q3 Sketch
# ============================================================

c = 0.5
Q1 = np.array([[3, 1],[1, 2]])
Q2 = Q1 / c
Q3 = Q1 / (c**2)

def plot_ellipse(Q, ax, label):
    theta = np.linspace(0, 2*np.pi, 400)
    unit_circle = np.vstack([np.cos(theta), np.sin(theta)])
    eigvals, eigvecs = np.linalg.eigh(Q)
    D_inv_sqrt = np.diag(1/np.sqrt(eigvals))
    transform = eigvecs @ D_inv_sqrt @ eigvecs.T
    ellipse = transform @ unit_circle
    ax.plot(ellipse[0], ellipse[1], label=label)

fig, ax = plt.subplots(figsize=(6,6))
plot_ellipse(Q1, ax, "Q1")
plot_ellipse(Q2, ax, "Q2")
plot_ellipse(Q3, ax, "Q3")
ax.set_aspect('equal')
ax.set_title("2D Sketch of Q1, Q2, Q3 for c = 0.5")
ax.legend()
plt.savefig("plots/fig_5.png")
plt.show()