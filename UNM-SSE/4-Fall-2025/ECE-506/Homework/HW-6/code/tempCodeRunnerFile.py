import torch

# ============================================================
# Function Definitions (used by later problems)
# ============================================================

def get_loss(c_pred, c_true):
    """Mean-squared error loss."""
    return torch.mean((c_pred - c_true)**2)

def predict(alpha0, alpha1, f):
    """Linear model: c = alpha0 + alpha1 * f."""
    return alpha0 + alpha1 * f


# ============================================================
# Problem 3 – Automatic Differentiation
# ============================================================

# Create x with requires_grad=True so autograd tracks operations
x = torch.ones(6, requires_grad=True)

# Define the function f(x) = (1/6) * sum (x_i - 2)^2 + x_i
f = (1/6) * torch.sum((x - 2)**2) + torch.sum(x)

# Backpropagate to compute gradients df/dx_i
f.backward()

# Print the gradient vector
print("Problem 3 gradient:", x.grad)


# ============================================================
# Problem 4(a) – Fahrenheit → Celsius using Autograd
# ============================================================

# Generate random Fahrenheit values
f_vals = 100 * torch.rand(10)

# Ground truth Celsius values
c_true = (f_vals - 32) / 1.8

# Initialize parameters
alpha0 = torch.zeros(1, requires_grad=True)
alpha1 = torch.zeros(1, requires_grad=True)

learning_rate = 0.001

# Training loop
for iteration in range(2000):

    # Forward pass
    c_pred = predict(alpha0, alpha1, f_vals)
    loss = get_loss(c_pred, c_true)

    # Backward pass
    loss.backward()

    # Parameter update
    with torch.no_grad():
        alpha0 -= learning_rate * alpha0.grad
        alpha1 -= learning_rate * alpha1.grad

    # Zero gradients for next iteration
    alpha0.grad.zero_()
    alpha1.grad.zero_()

# Display results
true_formula = "True:    C = -17.78 + 0.56·F"
learned_formula = f"Learned: C = {alpha0.item():.2f} + {alpha1.item():.2f}·F"

print("\nFahrenheit to Celsius")
print("\t" + true_formula)
print("\t" + learned_formula)
print("Final loss:", loss.item())
