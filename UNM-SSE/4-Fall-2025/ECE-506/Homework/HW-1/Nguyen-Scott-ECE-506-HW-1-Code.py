import numpy as np
import matplotlib.pyplot as plt

# 1c

t = np.linspace(0, 1, 200)
x1 = 2 - 2*t
x2 = t

plt.figure(figsize=(5,5))

# Shade feasible region (polygon: (0,0) -> (2,0) -> (0,1))
feasible_x = [0, 2, 0]
feasible_y = [0, 0, 1]
plt.fill(feasible_x, feasible_y, color='lightblue', alpha=0.5, label='Feasible Region')

# Plot the boundary line
plt.plot(x1, x2, linewidth=2, color='blue', label=r'$x_1 + 2x_2 = 2,\; x_1, x_2 \geq 0$')
plt.scatter([2, 0], [0, 1], s=60, color='blue')

# Labels for endpoints
plt.text(2, 0, r'  $(2,0)$', va='center')
plt.text(0, 1, r'  $(0,1)$', va='center')

# Axis labels and limits
plt.xlabel(r'$x_1$')
plt.ylabel(r'$x_2$')
plt.xlim(0, 2.1)
plt.ylim(0, 1.1)
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True, linestyle='--', alpha=0.5)

# Legend
plt.legend(loc='upper right')
plt.tight_layout()

# Save and show
plt.savefig('1c.png', dpi=200)
plt.show()