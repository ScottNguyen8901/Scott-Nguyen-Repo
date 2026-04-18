# -------------------------------------------------------------
# Homework 4 - Plot Generator
# Problems: 1(a), 1(b), 1(c), 1(e), 1(f), 1(h), 1(i), 1(j)
# -------------------------------------------------------------

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle

os.makedirs("plots", exist_ok=True)

# -------------------------------------------------------------
# Shared parameters
# -------------------------------------------------------------

T0 = 1.0
T1 = 1.0
T2 = 1.0

B1 = 0.2
B2 = 0.15

# -------------------------------------------------------------
# Helper functions
# -------------------------------------------------------------

def draw_1d_impulses(x_vals, xlabel, title, filename):
    plt.figure(figsize=(10, 3))
    for x in x_vals:
        plt.arrow(
            x, 0, 0, 1,
            head_width=0.08,
            head_length=0.15,
            length_includes_head=True
        )

    plt.axhline(0, linewidth=1)
    plt.xlim(min(x_vals) - 1, max(x_vals) + 1)
    plt.ylim(0, 1.2)
    plt.xlabel(xlabel)
    plt.title(title)
    plt.grid(True)
    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def draw_2d_impulse_lattice(T1, T2, filename):
    m_vals = np.arange(-5, 6)
    n_vals = np.arange(-5, 6)

    fig, ax = plt.subplots(figsize=(6, 6))

    for m in m_vals:
        for n in n_vals:
            x = m * T1
            y = n * T2
            ax.arrow(
                x, y, 0, 0.8,
                head_width=0.15,
                head_length=0.2,
                length_includes_head=True
            )

    ax.axhline(0, linewidth=1)
    ax.axvline(0, linewidth=1)
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("2D Impulse Train")
    ax.grid(True)

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def draw_spectral_copies_ellipses(dx, dy, filename):
    fig, ax = plt.subplots(figsize=(8, 8))
    spec_width = 0.55
    spec_height = 0.35

    for m in [-1, 0, 1]:
        for n in [-1, 0, 1]:
            cx = m * dx
            cy = n * dy
            ax.add_patch(Ellipse(
                (cx, cy),
                width=spec_width,
                height=spec_height,
                fill=False,
                linewidth=2.5 if (m, n) == (0, 0) else 1.8
            ))
            ax.plot(cx, cy, "k.", markersize=6)

    ax.axhline(0, linewidth=1)
    ax.axvline(0, linewidth=1)

    ax.text(0.03, -0.08, r'$(0,0)$', fontsize=10)
    ax.text(dx + 0.03, -0.08, r'$(\frac{1}{T_1},0)$', fontsize=10)
    ax.text(-dx - 0.38, -0.08, r'$(-\frac{1}{T_1},0)$', fontsize=10)
    ax.text(0.03, dy + 0.05, r'$(0,\frac{1}{T_2})$', fontsize=10)
    ax.text(0.03, -dy - 0.18, r'$(0,-\frac{1}{T_2})$', fontsize=10)
    ax.text(0, 0.12, r'$\hat{f}(\xi_1,\xi_2)$', ha='center', fontsize=11)

    ax.set_xlim(-1.6 * dx, 1.6 * dx)
    ax.set_ylim(-1.6 * dy, 1.6 * dy)
    ax.set_aspect("equal")
    ax.set_xlabel(r'$\xi_1$')
    ax.set_ylabel(r'$\xi_2$')
    ax.set_title("Frequency-Domain Spectral Copies")
    ax.grid(True)

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def draw_spectral_rectangles(dx, dy, B1, B2, title, filename):
    fig, ax = plt.subplots(figsize=(8, 8))

    for m in [-1, 0, 1]:
        for n in [-1, 0, 1]:
            cx = m * dx
            cy = n * dy
            ax.add_patch(Rectangle(
                (cx - B1, cy - B2),
                2 * B1,
                2 * B2,
                fill=False,
                linewidth=2 if (m, n) == (0, 0) else 1.5
            ))

    ax.axhline(0, linewidth=1)
    ax.axvline(0, linewidth=1)
    ax.set_xlim(-1.6 * dx, 1.6 * dx)
    ax.set_ylim(-1.6 * dy, 1.6 * dy)
    ax.set_aspect("equal")
    ax.set_xlabel(r'$\xi_1$')
    ax.set_ylabel(r'$\xi_2$')
    ax.set_title(title)
    ax.grid(True)

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()


def draw_filter_passband(filename):
    fig, ax = plt.subplots(figsize=(8, 8))

    orig_width = 1.0
    orig_height = 1.0
    pass_width = orig_width / 2
    pass_height = orig_height / 3

    outer = Rectangle(
        (-orig_width / 2, -orig_height / 2),
        orig_width,
        orig_height,
        fill=False,
        linewidth=2
    )
    ax.add_patch(outer)

    inner = Rectangle(
        (-pass_width / 2, -pass_height / 2),
        pass_width,
        pass_height,
        fill=False,
        linewidth=2,
        linestyle="--"
    )
    ax.add_patch(inner)

    ax.text(0, 0.05, "Kept Passband", ha="center")
    ax.text(0.35, 0.35, "Removed", fontsize=10)
    ax.text(-0.45, -0.35, "Removed", fontsize=10)

    ax.axhline(0, linewidth=1)
    ax.axvline(0, linewidth=1)
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.6, 0.6)
    ax.set_aspect("equal")
    ax.set_xlabel(r'$\xi_1$')
    ax.set_ylabel(r'$\xi_2$')
    ax.set_title("FFT Passband After Anti-Aliasing Filter")
    ax.grid(True)

    plt.savefig(filename, dpi=300, bbox_inches="tight")
    plt.close()

# -------------------------------------------------------------
# Problem 1(a): Impulse Train
# -------------------------------------------------------------

x_samples = np.arange(-12, 13) * T0
draw_1d_impulses(
    x_vals=x_samples,
    xlabel="x",
    title="Impulse Train",
    filename="plots/problem_1a_impulse_train.png"
)

# -------------------------------------------------------------
# Problem 1(b): 2D Impulse Train
# -------------------------------------------------------------

draw_2d_impulse_lattice(
    T1=T1,
    T2=T2,
    filename="plots/problem_1b_impulse_lattice.png"
)

# -------------------------------------------------------------
# Problem 1(c): Fourier-domain Impulse Train
# -------------------------------------------------------------

xi_samples = np.arange(-12, 13) / T0
draw_1d_impulses(
    x_vals=xi_samples,
    xlabel=r'$\xi$',
    title="Fourier Impulse Train",
    filename="plots/problem_1c_fourier_impulses.png"
)

# -------------------------------------------------------------
# Problem 1(e): Frequency-Domain Spectral Copies
# -------------------------------------------------------------

draw_spectral_copies_ellipses(
    dx=1 / T1,
    dy=1 / T2,
    filename="plots/problem_1e_spectral_copies.png"
)

# -------------------------------------------------------------
# Problem 1(f): Non-overlapping spectral copies
# -------------------------------------------------------------

draw_spectral_rectangles(
    dx=1 / T1,
    dy=1 / T2,
    B1=B1,
    B2=B2,
    title="Non-overlapping Spectral Copies",
    filename="plots/problem_1f_no_overlap.png"
)

# -------------------------------------------------------------
# Problem 1(h): Interpolation
# -------------------------------------------------------------

draw_spectral_rectangles(
    dx=2 / T1,
    dy=3 / T2,
    B1=B1,
    B2=B2,
    title="Interpolation: Spectral Copies Move Farther Apart",
    filename="plots/problem_1h_interpolation_spacing.png"
)

# -------------------------------------------------------------
# Problem 1(i): Downsampling / Aliasing
# -------------------------------------------------------------

draw_spectral_rectangles(
    dx=1 / (2 * T1),
    dy=1 / (3 * T2),
    B1=B1,
    B2=B2,
    title="Aliasing: Overlapping Spectral Copies",
    filename="plots/problem_1i_aliasing.png"
)

# -------------------------------------------------------------
# Problem 1(j): Digital filtering before downsampling
# -------------------------------------------------------------

draw_filter_passband(
    filename="plots/problem_1j_filter_passband.png"
)