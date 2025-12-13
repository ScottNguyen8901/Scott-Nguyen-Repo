import os
import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import make_blobs
from sklearn import svm
from sklearn.svm import LinearSVC
from numpy.linalg import eig

# Problem 1

N_dimensions=2
center_1=np.zeros(N_dimensions)
center_2=np.random.randn(N_dimensions)
center_2=center_2/np.linalg.norm(center_2)

X, y = make_blobs(n_samples=200, centers=[center_1, center_2], cluster_std=0.25)
y=2*(y-1/2)

plt.scatter(X[:, 0], X[:, 1], c=y)
plt.title('Two Gaussian Clusters in 2D  ')
plt.xlabel('Feature 1')
plt.ylabel('Feature 2')

os.makedirs("plots", exist_ok=True)
plt.savefig("plots/p1.png", dpi=300, bbox_inches="tight")
plt.show()

# Problem 2

clf = svm.SVC(kernel='linear')
clf.fit(X, y)
y_pred = clf.predict(X)

w = clf.coef_[0]
b = clf.intercept_[0]

plt.figure(figsize=(6,5))
plt.scatter(X[:, 0], X[:, 1], c=y, cmap='bwr', edgecolors='k', s=40)

xx = np.linspace(X[:, 0].min()-0.5, X[:, 0].max()+0.5, 200)
yy = (-w[0] * xx - b) / w[1]
yy_plus  = (-w[0] * xx - b + 1) / w[1]
yy_minus = (-w[0] * xx - b - 1) / w[1]

plt.plot(xx, yy, 'k-', label='Decision boundary')
plt.plot(xx, yy_plus, 'k--', label='Margin +1')
plt.plot(xx, yy_minus, 'k--', label='Margin -1')

plt.scatter(
    clf.support_vectors_[:, 0],
    clf.support_vectors_[:, 1],
    s=120, facecolors='none', edgecolors='k', label='Support vectors'
)

plt.title('Linear SVM with Boundary, Margins, and Support Vectors')
plt.xlabel('Feature 1')
plt.ylabel('Feature 2')
plt.legend(loc='best')
plt.axis('equal')

plt.savefig("plots/p2.png", dpi=300, bbox_inches="tight")
plt.show()

# Problem 3.1

def data(T, N_data, N_dimensions):
    sigma = 0.15
    centers = np.concatenate(
        (np.ones((1, T.shape[0])), np.eye(T.shape[0])),
        axis=0
    ) @ T
    X, y = make_blobs(
        n_samples=N_data,
        centers=centers,
        cluster_std=np.sqrt(sigma) / np.sqrt(N_dimensions)
    )
    y[np.where(np.mod(y, 2) == 1)] = -1
    y[np.where(np.mod(y, 2) == 0)] = 1
    return X, y

def run_once(N_clusters, N_dim, N_train, N_val, C_grid):
    A = np.random.randn(N_dim, N_dim)
    _, V = eig(A @ A.T)
    T = V[0:N_clusters-1, :]

    X_train, y_train = data(T, N_train, N_dim)
    X_val,   y_val   = data(T, N_val, N_dim)

    best_C, best_err = None, 1.0

    for C in C_grid:
        clf = LinearSVC(C=C, max_iter=2000, dual=False)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_val)
        err = np.mean(y_pred != y_val)
        if err < best_err:
            best_err, best_C = err, C

    return best_C, best_err

N_clusters = 8
N_val = 100 * N_clusters
C_grid = np.logspace(-1.5, 1, 15)
N_repeats = 100

setups = [
    (8, 7,   2*8),
    (8, 70,  2*8),
    (8, 7,   100*8),
    (8, 70,  100*8),
]

results = []

for (Nc, Nd, Ntr) in setups:
    print(f"\nRunning setup: N_dim={Nd}, N_train={Ntr} ...")
    C_list, E_list = [], []

    for i in range(N_repeats):
        C_opt, E_opt = run_once(Nc, Nd, Ntr, N_val, C_grid)
        C_list.append(C_opt)
        E_list.append(E_opt)

        pct = ((i + 1) / N_repeats) * 100
        print(f"Progress: {pct:6.2f}% ({i+1}/{N_repeats})", end="\r")

    print()
    results.append({
        "N_dim": Nd,
        "N_train": Ntr,
        "C_opt_avg": np.mean(C_list),
        "E_opt_avg": np.mean(E_list)
    })

print("\nAverage results over", N_repeats, "runs:")
print(" N_dim | N_train |  <C_opt>   |  <E_opt>")
print("------------------------------------------")
for r in results:
    print(f"{r['N_dim']:6d} | {r['N_train']:8d} | {r['C_opt_avg']:9.4f} | {r['E_opt_avg']:.4f}")