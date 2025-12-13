import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Problem 4: Multiobjective Optimization (Figure for Problem #4)
#   - Generates plots/fig_1.png
#   - objective2 on x-axis, objective1 on y-axis
# ============================================================

points = [
    (2, 5),
    (2, 7),
    (3, 3),
    (3, 4),
    (3, 6),
    (3, 8),
    (4, 1),
    (4, 6),
    (4, 8),
    (6, 1),
    (6, 5),
    (7, 3),
    (8, 1),
]

x = [p[0] for p in points]  # objective 2
y = [p[1] for p in points]  # objective 1

plt.figure(figsize=(6, 4))
plt.scatter(x, y, marker='x', s=80)

plt.xlabel("objective 2")
plt.ylabel("objective 1")

plt.xlim(0, 9)
plt.ylim(0, 9)

plt.xticks(range(0, 10))
plt.yticks(range(0, 10))

plt.grid(False)
plt.tight_layout()

# Save exactly where LaTeX expects it
plt.savefig("plots/fig_1.png", dpi=300)
plt.close()


# ============================================================
# Problem 6(a): Convex, nonnegative function region
#   - Region where convex f with f(x) >= 0 must lie
#   - Generates plots/fig_5.png
# ============================================================

# Example endpoints and function values
x_k = 1.0
x_k1 = 4.0

f_k = 2.0
f_k1 = 4.5

# Slope of the secant line between (x_k, f_k) and (x_k1, f_k1)
secant_slope = (f_k1 - f_k) / (x_k1 - x_k)

# Choose a tangent slope at x_k that is <= secant_slope (for convexity)
grad_k = 0.5 * secant_slope

# Build grid on [x_k, x_k1]
xs = np.linspace(x_k, x_k1, 200)

# Tangent line at x_k: lower bound
tangent = f_k + grad_k * (xs - x_k)

# Secant line between (x_k, f_k) and (x_k1, f_k1): upper bound
secant = f_k + secant_slope * (xs - x_k)

# Non-negativity bound
zero_line = np.zeros_like(xs)

# Region where convex, nonnegative f must lie
lower_bound = np.maximum(zero_line, tangent)
upper_bound = secant

plt.figure(figsize=(6, 4))

# Fill feasible region for f between x_k and x_k1
plt.fill_between(xs, lower_bound, upper_bound,
                 alpha=0.2, hatch='///', edgecolor='k')

# Plot tangent, secant, and f >= 0
plt.plot(xs, tangent, label=r'Tangent at $x_k$', linestyle='--')
plt.plot(xs, secant, label='Line segment', linestyle='-')
plt.axhline(0, color='black', linewidth=1, linestyle=':')

# Mark the endpoints (x_k, f_k) and (x_{k+1}, f_{k+1})
plt.scatter([x_k, x_k1], [f_k, f_k1], color='red', zorder=5)
plt.text(x_k, f_k + 0.15, r'$(x_k, f(x_k))$', ha='center')
plt.text(x_k1, f_k1 + 0.15, r'$(x_{k+1}, f(x_{k+1}))$', ha='center')

plt.xlabel(r'$x$')
plt.ylabel(r'$f(x)$')
plt.title(r'Region where convex, nonnegative $f$ must lie on $[x_k, x_{k+1}]$')
plt.legend(loc='upper left')

plt.xlim(x_k - 0.5, x_k1 + 0.5)
plt.ylim(0, max(f_k1, np.max(secant)) + 1)

plt.tight_layout()
plt.savefig('plots/fig_5.png', dpi=300)
plt.close()


# ============================================================
# Problem 6(b): Concave function with 0 <= f(x) <= U
#   - Region where concave f with bounds must lie
#   - Generates plots/fig_6.png
# ============================================================

# Example endpoints and values
x_k = 1.0
x_k1 = 4.0

f_k = 3.0
f_k1 = 2.0
U = 5.0  # upper bound on f(x)

xs = np.linspace(x_k, x_k1, 200)

# Secant line: lower bound for concave f
secant = f_k + (f_k1 - f_k) * (xs - x_k) / (x_k1 - x_k)

# Tangent line at x_k: upper bound (slope chosen to be consistent with concavity)
grad_k = -0.3
tangent = f_k + grad_k * (xs - x_k)

# Horizontal bounds 0 <= f <= U
zero_line = np.zeros_like(xs)
upper_cap = U * np.ones_like(xs)

upper_bound = np.minimum(tangent, upper_cap)
lower_bound = np.maximum(secant, zero_line)

plt.figure(figsize=(6, 4))

plt.fill_between(xs, lower_bound, upper_bound,
                 alpha=0.25, hatch='\\\\', edgecolor='k')

plt.plot(xs, secant, label='Secant line (lower bound)', linestyle='-')
plt.plot(xs, tangent, label=r'Tangent at $x_k$ (upper bound)', linestyle='--')
plt.axhline(0, color='black', linestyle=':')
plt.axhline(U, color='black', linestyle=':')

plt.scatter([x_k, x_k1], [f_k, f_k1], color='red', zorder=5)
plt.text(x_k, f_k + 0.15, r'$(x_k, f(x_k))$', ha='center')
plt.text(x_k1, f_k1 + 0.15, r'$(x_{k+1}, f(x_{k+1}))$', ha='center')

plt.xlabel(r'$x$')
plt.ylabel(r'$f(x)$')
# IMPORTANT: no math here, just plain text to avoid \le parsing issues
plt.title('Concave function region with bounds 0 <= f <= U')
plt.legend()

plt.ylim(0, U + 1)
plt.xlim(x_k - 0.5, x_k1 + 0.5)

plt.tight_layout()
plt.savefig('plots/fig_6.png', dpi=300)
plt.close()
