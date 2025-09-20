import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linprog

# --- Console: allow Unicode (for ∇) and set 2-sig-fig formatter ---
sys.stdout.reconfigure(encoding="utf-8")
def fmt(x): return f"{float(x):.2g}"
def fmt_vec(v): return "[" + ", ".join(fmt(t) for t in np.atleast_1d(v)) + "]"

# --- Plot helpers (consistent style across problems) ---
def setup_square(xmin, xmax, ymin, ymax, title="", xlabel=r"$x_1$", ylabel=r"$x_2$"):
    plt.figure(figsize=(6, 6))
    if title:
        plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.grid(True, linestyle="--", alpha=0.5)

def save_png(name):
    plt.tight_layout()
    plt.savefig(name, dpi=200)

# =========================
# Problem 1(c)
# =========================
def problem_1c():
    # Feasible region: x1 >= 0, x2 >= 0, x1 + 2 x2 <= 2 (triangle with (0,0),(2,0),(0,1))
    setup_square(
        0, 2.1, 0, 1.1,
        title=r"Feasible region for $x_1\geq 0,\ x_2\geq 0,\ x_1+2x_2\leq 2$"
    )
    # Shade feasible polygon
    plt.fill([0, 2, 0], [0, 0, 1], color="lightblue", alpha=0.5, label="Feasible region")
    # Boundary line x1 + 2 x2 = 2 (endpoints)
    plt.plot([2, 0], [0, 1], linewidth=2, label=r"$x_1+2x_2=2$")
    plt.scatter([2, 0], [0, 1], s=60)
    plt.text(2, 0, r"  $(2,0)$", va="center")
    plt.text(0, 1, r"  $(0,1)$", va="center")
    plt.legend(loc="upper right")
    save_png("1c.png")

# =========================
# Problem 3(c)
# =========================
def problem_3c():
    # s(x) = c^T ReLU(Ax+b)
    A = np.array([[-1.0,  1.0],
                  [ 0.5,  0.5]])
    b = np.array([0.0, 1.0])
    c = np.array([1.0, 2.0])

    x1_min, x1_max = -2.0, 2.0
    x2_min, x2_max = -2.0, 2.0
    res = 400
    x1 = np.linspace(x1_min, x1_max, res)
    x2 = np.linspace(x2_min, x2_max, res)
    X1, X2 = np.meshgrid(x1, x2)
    P = np.stack([X1.ravel(), X2.ravel()], axis=1)

    H = (P @ A.T) + b       # pre-activations
    Y = np.maximum(0.0, H)  # ReLU
    S = (Y @ c).reshape(res, res)
    H_maps = [H[:, i].reshape(res, res) for i in range(A.shape[0])]

    # (i) Hidden-region boundaries and active sides
    setup_square(x1_min, x1_max, x2_min, x2_max,
                 title=r"Hidden regions: $a_i^{\top}x+b_i=0$ and active sides (ReLU)")
    for i, H_i in enumerate(H_maps):
        plt.contourf(X1, X2, (H_i > 0).astype(float), levels=[-0.5, 0.5, 1.5], alpha=0.15)
        cs = plt.contour(X1, X2, H_i, levels=[0.0], linewidths=2)
        if cs.allsegs[0]:
            plt.clabel(cs, fmt={0.0: f"h_{i+1}=0"}, inline=True, fontsize=9)
    save_png("3c_i.png")

    # (ii) Output surface s(x) heatmap + contours, with hidden boundaries overlaid
    plt.figure(figsize=(7, 6))
    plt.title(r"Output surface $s(x)=c^{\top}\mathrm{ReLU}(Ax+b)$")
    im = plt.imshow(S, extent=[x1_min, x1_max, x2_min, x2_max],
                    origin="lower", aspect="equal")
    plt.colorbar(im, label=r"$s(x)$")
    CS = plt.contour(X1, X2, S, colors="k", linewidths=0.8, levels=10)
    plt.clabel(CS, inline=True, fontsize=8, fmt="%.2f")
    for H_i in H_maps:
        plt.contour(X1, X2, H_i, levels=[0.0], colors="white", linewidths=1.2, alpha=0.9)
    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.xlim(x1_min, x1_max)
    plt.ylim(x2_min, x2_max)
    plt.grid(True, linestyle="--", alpha=0.4)
    plt.tight_layout()
    plt.savefig("3c_ii.png", dpi=200)

# =========================
# Problem 3(d)
# =========================
def problem_3d():
    A = np.array([[-1.0,  1.0],
                  [ 0.5,  0.5]])
    b = np.array([0.0, 1.0])

    x1_min, x1_max = -4.0, 4.0
    x2_min, x2_max = -4.0, 4.0
    res = 401
    x1 = np.linspace(x1_min, x1_max, res)
    x2 = np.linspace(x2_min, x2_max, res)
    X1, X2 = np.meshgrid(x1, x2)

    Y1 = -X1 + X2 + b[0]           # y1 >= 0 -> x2 >= x1
    Y2 = 0.5*X1 + 0.5*X2 + b[1]    # y2 >= 0 -> x1 + x2 >= -2

    feasible = (Y1 >= 0) & (Y2 >= 0)

    setup_square(x1_min, x1_max, x2_min, x2_max,
                 title=r"3(d) Feasible region: $x_2\geq x_1$ and $x_1+x_2\geq -2$")
    plt.contourf(X1, X2, feasible.astype(float), levels=[-0.5, 0.5, 1.5], alpha=0.3)
    C1 = plt.contour(X1, X2, Y1, levels=[0.0], linewidths=2)
    C2 = plt.contour(X1, X2, Y2, levels=[0.0], linewidths=2, linestyles="--")
    if C1.allsegs[0]:
        plt.clabel(C1, fmt={0.0: r"$x_2=x_1$"}, inline=True, fontsize=9)
    if C2.allsegs[0]:
        plt.clabel(C2, fmt={0.0: r"$x_1+x_2=-2$"}, inline=True, fontsize=9)
    corner = (-1.0, -1.0)
    plt.scatter([corner[0]], [corner[1]], s=60)
    plt.text(corner[0], corner[1], r"  $(-1,-1)$", va="center")
    save_png("3d.png")

# =========================
# Problem 3(e)
# =========================
def problem_3e():
    # Minimize [0,2]·x s.t. x1 - x2 <= 0, -x1 - x2 <= 2
    A_ub = np.array([[ 1.0, -1.0],
                     [-1.0, -1.0]])
    b_ub = np.array([0.0, 2.0])
    c = np.array([0.0, 2.0])
    bounds = [(None, None), (None, None)]

    res = linprog(c, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method="highs")
    A = np.array([[-1.0,  1.0],
                  [ 0.5,  0.5]])
    b = np.array([0.0, 1.0])
    y = A @ res.x + b
    orig_obj = np.array([1.0, 2.0]) @ y

    print("=== Problem 3(e) Answers ===")
    print(f"status: {res.message}")
    print(f"x*     : {fmt_vec(res.x)}")
    print(f"obj c^T x: {fmt(res.fun)}")
    print(f"y      : {fmt_vec(y)}")
    print(f"[1 2]^T y: {fmt(orig_obj)}")
    print(f"y >= 0?  {np.all(y >= -1e-10)}")

# =========================
# Problem 4(a)
# =========================
def problem_4a():
    x1_min, x1_max = -1.5, 1.5
    x2_min, x2_max = -1.5, 1.5
    res = 801
    x1 = np.linspace(x1_min, x1_max, res)
    x2 = np.linspace(x2_min, x2_max, res)
    X1, X2 = np.meshgrid(x1, x2)

    L1 = np.abs(X1) + np.abs(X2)
    feasible = (L1 <= 1.0)

    setup_square(x1_min, x1_max, x2_min, x2_max,
                 title=r"$\ell_1$ feasible region: $|x_1|+|x_2|\leq 1$")
    plt.contourf(X1, X2, feasible.astype(float), levels=[-0.5, 0.5, 1.5], alpha=0.30)
    C = plt.contour(X1, X2, L1, levels=[1.0], linewidths=2)
    if C.allsegs[0]:
        plt.clabel(C, fmt={1.0: r"$|x_1|+|x_2|=1$"}, inline=True, fontsize=9)

    verts = np.array([[ 1, 0], [ 0, 1], [-1, 0], [ 0,-1]])
    labels = [r"$(1,0)$", r"$(0,1)$", r"$(-1,0)$", r"$(0,-1)$"]
    plt.scatter(verts[:, 0], verts[:, 1], s=50)
    for (xv, yv), lab in zip(verts, labels):
        plt.text(xv, yv, "  " + lab, va="center", ha="left")
    save_png("4a.png")

# =========================
# Problem 5
# =========================
def problem_5():
    def f(x): x1, x2 = x; return (x1 - 2.0)**2 + 10.0*(x2 - 3.0)**2
    def grad_f(x): x1, x2 = x; return np.array([2.0*(x1 - 2.0), 20.0*(x2 - 3.0)])
    H = np.array([[2.0, 0.0], [0.0, 20.0]])

    x_star_analytic = np.array([2.0, 3.0])
    grad_at_star = grad_f(x_star_analytic)
    f_at_star = f(x_star_analytic)
    eigvals = np.linalg.eigvalsh(H)
    is_pd = np.all(eigvals > 0)

    # Simple GD (demo)
    x = np.array([-3.0, 5.0]); alpha = 0.1; max_iters = 200; tol = 1e-10
    for k in range(max_iters):
        g = grad_f(x)
        if np.linalg.norm(g) < tol: break
        x -= alpha * g
    x_star_gd = x; f_star_gd = f(x_star_gd); grad_at_gd = grad_f(x_star_gd)

    print("=== Problem 5 Answers ===")
    print("Analytic stationary point (∇f = 0):")
    print(f"x* (analytic)       = {fmt_vec(x_star_analytic)}")
    print(f"∇f(x*) (analytic)   = {fmt_vec(grad_at_star)}")
    print(f"f(x*) (analytic)    = {fmt(f_at_star)}")
    print(f"H PD? {is_pd} (eigs = {fmt_vec(eigvals)})\n")
    print("Gradient-descent solution (demo):")
    print(f"x* (GD)             = {fmt_vec(x_star_gd)}  (iters <= {k+1})")
    print(f"∇f(x*) (GD)         = {fmt_vec(grad_at_gd)}")
    print(f"f(x*) (GD)          = {fmt(f_star_gd)}")

    # Contours + gradient field
    x1 = np.linspace(-2.0, 6.0, 200)
    x2 = np.linspace(-1.0, 7.0, 200)
    X1, X2 = np.meshgrid(x1, x2)
    Z = (X1 - 2.0)**2 + 10.0*(X2 - 3.0)**2

    skip = 8
    X1s, X2s = X1[::skip, ::skip], X2[::skip, ::skip]
    U, V = 2.0*(X1s - 2.0), 20.0*(X2s - 3.0)

    plt.figure(figsize=(7, 6))
    cs = plt.contour(X1, X2, Z, levels=15)
    plt.clabel(cs, inline=True, fontsize=8)
    plt.quiver(X1s, X2s, U, V, angles="xy", scale_units="xy", scale=60)
    plt.scatter([x_star_analytic[0]], [x_star_analytic[1]], s=70, marker="o",
                label="Analytic minimizer (2, 3)")
    plt.scatter([x_star_gd[0]], [x_star_gd[1]], s=50, marker="x", label="GD solution")
    plt.title("Contours of f(x1,x2) and gradient field")
    plt.xlabel("x1"); plt.ylabel("x2")
    plt.xlim(x1.min(), x1.max()); plt.ylim(x2.min(), x2.max())
    plt.gca().set_aspect("equal", adjustable="box")
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig("5.png", dpi=200)

def problem_6a():
    # -------- Problem 6(a): f(x1,x2) = (x1 - 1)^2 + (x2 - 4)^2 --------

    # Two-sig-fig formatting (consistent with other problems)
    def fmt(x): return f"{float(x):.2g}"
    def fmt_vec(v): return "[" + ", ".join(fmt(t) for t in np.atleast_1d(v)) + "]"

    # Objective
    def f(x):
        x1, x2 = x
        return (x1 - 1.0)**2 + (x2 - 4.0)**2

    # Euclidean projection onto the l1 ball of radius 1
    def project_onto_l1(u, radius=1.0):
        u = np.asarray(u, dtype=float)
        if np.sum(np.abs(u)) <= radius:
            return u.copy()
        a = np.abs(u)
        s = np.sort(a)[::-1]
        cssv = np.cumsum(s)
        rho = np.max(np.where(s - (cssv - radius) / (np.arange(1, len(s)+1)) > 0)[0]) + 1
        tau = (cssv[rho-1] - radius) / rho
        return np.sign(u) * np.maximum(a - tau, 0.0)

    # Unconstrained minimizer: ∇f=0 -> (1,4)
    x_u = np.array([1.0, 4.0])

    # Constrained minimizer: projection of x_u onto |x1|+|x2| <= 1
    x_c = project_onto_l1(x_u, radius=1.0)

    # Print answers
    print("=== Problem 6(a) Answers ===")
    print("Unconstrained minimizer:")
    print(f"x_u          = {fmt_vec(x_u)}")
    print(f"f(x_u)       = {fmt(f(x_u))}")
    print(f"||x_u||_1    = {fmt(np.sum(np.abs(x_u)))}")
    print("\nConstrained minimizer (projection onto |x1|+|x2|<=1):")
    print(f"x_c          = {fmt_vec(x_c)}")
    print(f"f(x_c)       = {fmt(f(x_c))}")
    print(f"||x_c||_1    = {fmt(np.sum(np.abs(x_c)))}  (should be 1 if on boundary)")

    # -------- Plot: contours + l1 constraint + points --------
    x1_min, x1_max = -1.5, 5.0
    x2_min, x2_max = -1.5, 5.0
    res = 400
    x1 = np.linspace(x1_min, x1_max, res)
    x2 = np.linspace(x2_min, x2_max, res)
    X1, X2 = np.meshgrid(x1, x2)

    Z  = (X1 - 1.0)**2 + (X2 - 4.0)**2                          # f(x)
    L1 = np.abs(X1) + np.abs(X2)                                # |x1|+|x2|
    feasible = (L1 <= 1.0)

    plt.figure(figsize=(7, 6))
    plt.title(r"6(a): Contours of $f(x_1,x_2)$ and $\ell_1$ constraint $|x_1|+|x_2|\leq 1$")

    # Contours of f
    cs = plt.contour(X1, X2, Z, levels=20)
    plt.clabel(cs, inline=True, fontsize=8)

    # Shade feasible l1 region and draw boundary
    plt.contourf(X1, X2, feasible.astype(float), levels=[-0.5, 0.5, 1.5], alpha=0.25)
    C = plt.contour(X1, X2, L1, levels=[1.0], colors='k', linewidths=2)
    if C.allsegs[0]:
        plt.clabel(C, fmt={1.0: r"$|x_1|+|x_2|=1$"}, inline=True, fontsize=9)

    # Mark points and projection segment
    plt.scatter([x_u[0]], [x_u[1]], s=70, marker='o', label=r'Unconstrained $(1,4)$')
    plt.scatter([x_c[0]], [x_c[1]], s=70, marker='x', label='Constrained optimum')
    plt.plot([x_u[0], x_c[0]], [x_u[1], x_c[1]], linestyle='--', linewidth=1.2)

    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.xlim(x1_min, x1_max)
    plt.ylim(x2_min, x2_max)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig("6a.png", dpi=200)

def problem_6b():
    # -------- Problem 6(b): f(x1,x2) = (x1 - 5)^2 + x2^2 --------

    # Two-sig-fig formatting
    def fmt(x): return f"{float(x):.2g}"
    def fmt_vec(v): return "[" + ", ".join(fmt(t) for t in np.atleast_1d(v)) + "]"

    # Objective
    def f(x):
        x1, x2 = x
        return (x1 - 5.0)**2 + (x2 - 0.0)**2

    # Euclidean projection onto the l1 ball of radius 1
    def project_onto_l1(u, radius=1.0):
        u = np.asarray(u, dtype=float)
        if np.sum(np.abs(u)) <= radius:
            return u.copy()
        a = np.abs(u)
        s = np.sort(a)[::-1]
        cssv = np.cumsum(s)
        rho = np.max(np.where(s - (cssv - radius) / (np.arange(1, len(s)+1)) > 0)[0]) + 1
        tau = (cssv[rho-1] - radius) / rho
        return np.sign(u) * np.maximum(a - tau, 0.0)

    # Unconstrained minimizer: ∇f=0 -> (5,0)
    x_u = np.array([5.0, 0.0])

    # Constrained minimizer: projection of x_u onto |x1|+|x2| <= 1
    x_c = project_onto_l1(x_u, radius=1.0)

    # Print answers
    print("=== Problem 6(b) Answers ===")
    print("Unconstrained minimizer:")
    print(f"x_u          = {fmt_vec(x_u)}")
    print(f"f(x_u)       = {fmt(f(x_u))}")
    print(f"||x_u||_1    = {fmt(np.sum(np.abs(x_u)))}")
    print("\nConstrained minimizer (projection onto |x1|+|x2|<=1):")
    print(f"x_c          = {fmt_vec(x_c)}")
    print(f"f(x_c)       = {fmt(f(x_c))}")
    print(f"||x_c||_1    = {fmt(np.sum(np.abs(x_c)))}  (should be 1 if on boundary)")

    # -------- Plot: contours + l1 constraint + points --------
    x1_min, x1_max = -1.5, 6.0
    x2_min, x2_max = -1.5, 1.5
    res = 400
    x1 = np.linspace(x1_min, x1_max, res)
    x2 = np.linspace(x2_min, x2_max, res)
    X1, X2 = np.meshgrid(x1, x2)

    Z  = (X1 - 5.0)**2 + (X2 - 0.0)**2                     # f(x)
    L1 = np.abs(X1) + np.abs(X2)                            # |x1|+|x2|
    feasible = (L1 <= 1.0)

    plt.figure(figsize=(7, 6))
    plt.title(r"6(b): Contours of $f(x_1,x_2)$ and $\ell_1$ constraint $|x_1|+|x_2|\leq 1$")

    # Contours of f
    cs = plt.contour(X1, X2, Z, levels=20)
    plt.clabel(cs, inline=True, fontsize=8)

    # Shade feasible l1 region and draw boundary
    plt.contourf(X1, X2, feasible.astype(float), levels=[-0.5, 0.5, 1.5], alpha=0.25)
    C = plt.contour(X1, X2, L1, levels=[1.0], colors='k', linewidths=2)
    if C.allsegs[0]:
        plt.clabel(C, fmt={1.0: r"$|x_1|+|x_2|=1$"}, inline=True, fontsize=9)

    # Mark points and projection segment
    plt.scatter([x_u[0]], [x_u[1]], s=70, marker='o', label=r'Unconstrained $(5,0)$')
    plt.scatter([x_c[0]], [x_c[1]], s=70, marker='x', label='Constrained optimum')
    plt.plot([x_u[0], x_c[0]], [x_u[1], x_c[1]], linestyle='--', linewidth=1.2)

    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.xlim(x1_min, x1_max)
    plt.ylim(x2_min, x2_max)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig("6b.png", dpi=200)

def problem_6c():
    # -------- Problem 6(c): f(x1,x2) = (x1 - 0.5)^2 + x2^2 --------

    # Two-sig-fig formatting
    def fmt(x): return f"{float(x):.2g}"
    def fmt_vec(v): return "[" + ", ".join(fmt(t) for t in np.atleast_1d(v)) + "]"

    # Objective
    def f(x):
        x1, x2 = x
        return (x1 - 0.5)**2 + (x2 - 0.0)**2

    # Unconstrained minimizer: ∇f=0 -> (0.5, 0)
    x_u = np.array([0.5, 0.0])
    l1_u = np.sum(np.abs(x_u))
    on_boundary = np.isclose(l1_u, 1.0, atol=1e-12)

    # Since ||x_u||_1 = 0.5 < 1, the constrained optimum equals x_u (interior point)
    x_c = x_u.copy()

    # Print answers
    print("=== Problem 6(c) Answers ===")
    print("Unconstrained minimizer:")
    print(f"x_u          = {fmt_vec(x_u)}")
    print(f"f(x_u)       = {fmt(f(x_u))}")
    print(f"||x_u||_1    = {fmt(l1_u)}  (<= 1, interior)")
    print("\nConstrained minimizer:")
    print(f"x_c          = {fmt_vec(x_c)}  (same as x_u; interior)")
    print(f"f(x_c)       = {fmt(f(x_c))}")
    print(f"on boundary? = {on_boundary}")

    # -------- Plot: contours + l1 constraint + point --------
    x1_min, x1_max = -1.5, 2.0
    x2_min, x2_max = -1.5, 1.5
    res = 400
    x1 = np.linspace(x1_min, x1_max, res)
    x2 = np.linspace(x2_min, x2_max, res)
    X1, X2 = np.meshgrid(x1, x2)

    Z  = (X1 - 0.5)**2 + (X2 - 0.0)**2                 # f(x)
    L1 = np.abs(X1) + np.abs(X2)                        # |x1|+|x2|
    feasible = (L1 <= 1.0)

    plt.figure(figsize=(7, 6))
    plt.title(r"6(c): Contours of $f(x_1,x_2)$ and $\ell_1$ constraint $|x_1|+|x_2|\leq 1$")

    # Contours of f
    cs = plt.contour(X1, X2, Z, levels=20)
    plt.clabel(cs, inline=True, fontsize=8)

    # Shade feasible l1 region and draw boundary
    plt.contourf(X1, X2, feasible.astype(float), levels=[-0.5, 0.5, 1.5], alpha=0.25)
    C = plt.contour(X1, X2, L1, levels=[1.0], colors='k', linewidths=2)
    if C.allsegs[0]:
        plt.clabel(C, fmt={1.0: r"$|x_1|+|x_2|=1$"}, inline=True, fontsize=9)

    # Mark interior optimum
    plt.scatter([x_u[0]], [x_u[1]], s=70, marker='o', label=r'Optimum $(0.5,0)$ (interior)')

    plt.xlabel(r"$x_1$")
    plt.ylabel(r"$x_2$")
    plt.xlim(x1_min, x1_max)
    plt.ylim(x2_min, x2_max)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig("6c.png", dpi=200)

# =========================
# Run all
# =========================
if __name__ == "__main__":
    problem_1c()
    problem_3c()
    problem_3d()
    problem_3e()
    problem_4a()
    problem_5()
    problem_6a()
    problem_6b()
    problem_6c()