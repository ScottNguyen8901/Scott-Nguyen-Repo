from scipy import signal
import matplotlib.pyplot as plt
from sklearn.svm import SVR
import numpy as np
import os

os.makedirs("plots", exist_ok=True)
np.random.seed(0)

def buffer(x, n, hop):  # Generated with AI. This returns an array of buffered signal windows
  """
  Args:
    x: Input array.
    n: Length of each buffer.
    hop: Hop size (The overlap is n-hop).

  Returns:
    A 2D array where each row is a buffer.
  """
  if len(x) < n:
    raise ValueError("Input array length must be greater than or equal to the buffer length.")
  num_buffers = (len(x) - n) // hop + 1
  result = np.zeros((num_buffers, n))
  for i in range(num_buffers):
    start = i * hop
    end = start + n
    result[i, :] = x[start:end]
  return result

def signal_gen(sigma, N):
  T = 200
  sos = signal.butter(10, 0.2, 'low', output='sos')
  z = np.random.randn(N+T)
  x = signal.sosfilt(sos, z)
  x = x + sigma*np.random.randn(x.shape[0])
  x = x[T:]
  return x

def add_outliers(x, p, a):
  x = x + a*(np.random.rand(x.shape[0]) > (1-p))*np.random.randn(x.shape[0])
  return x

def build_datasets(N=100, W=10, H=1, sigma=0.0):
  xtrain = signal_gen(sigma, N)
  xval   = signal_gen(sigma, N)
  xtest  = signal_gen(sigma, N)

  Xtrain_full = buffer(xtrain, W, 1)
  Xval_full   = buffer(xval,   W, 1)
  Xtest_full  = buffer(xtest,  W, 1)

  ytrain_clean = xtrain[W+H-1:]
  yval_clean   = xval[W+H-1:]
  ytest_clean  = xtest[W+H-1:]

  Xtrain = Xtrain_full[0:-H]
  Xval   = Xval_full[0:-H]
  Xtest  = Xtest_full[0:-H]

  return Xtrain, ytrain_clean, Xval, yval_clean, Xtest, ytest_clean

def corrcoef(y_true, y_pred):
  return np.corrcoef(y_true, y_pred)[0, 1]

# 5.1A
N = 100
W = 10
H = 1
sigma = 0.0

xtrain = signal_gen(sigma, N)
xtest  = signal_gen(sigma, N)

Xtrain_full = buffer(xtrain, W, 1)
ytrain_clean = xtrain[W+H-1:]
Xtrain = Xtrain_full[0:-H]

print("\n5.1A")
print("Xtrain shape:", Xtrain.shape)
print("ytrain_clean shape:", ytrain_clean.shape)

Xtest_full = buffer(xtest, W, 1)
ytest_clean = xtest[W+H-1:]
Xtest = Xtest_full[0:-H]

print("Xtest shape:", Xtest.shape)
print("ytest_clean shape:", ytest_clean.shape)

# 5.1B
p = 0.1
a = 1.0
ytrain = add_outliers(ytrain_clean.copy(), p, a)

X = Xtrain.T
y = ytrain
XXT = X @ X.T
Xy = X @ y

w_mmse = np.linalg.solve(XXT, Xy)
w_ridge = np.linalg.solve(XXT + 0.2*np.eye(XXT.shape[0]), Xy)

svr_model = SVR(C=10.0, epsilon=0.01, kernel='linear')
svr_model.fit(Xtrain, ytrain)

yhat_mmse  = Xtest @ w_mmse
yhat_ridge = Xtest @ w_ridge
yhat_svr   = svr_model.predict(Xtest)

print("\n5.1B")
print("Models trained: MMSE, Ridge, SVR")

# 5.1C
t = np.arange(len(ytest_clean))
plt.figure()
plt.plot(t, ytest_clean, label='Clean', linewidth=2)
plt.plot(t, yhat_mmse,  '--', label='MMSE')
plt.plot(t, yhat_ridge, '--', label='Ridge')
plt.plot(t, yhat_svr,   '--', label='SVR')
plt.xlabel('Time')
plt.ylabel('Signal')
plt.title('Test Predictions')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("plots/predictions.png")
plt.close()
print("5.1C: saved plots/predictions.png")

# 5.1D
se_mmse  = (ytest_clean - yhat_mmse)**2
se_ridge = (ytest_clean - yhat_ridge)**2
se_svr   = (ytest_clean - yhat_svr)**2

plt.figure()
plt.plot(t, se_mmse,  label='MMSE')
plt.plot(t, se_ridge, label='Ridge')
plt.plot(t, se_svr,   label='SVR')
plt.xlabel('Time')
plt.ylabel('Squared Error')
plt.title('Squared Error Over Time')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("plots/squared_error.png")
plt.close()
print("5.1D: saved plots/squared_error.png")

# 5.2A
print("\n5.2A")

p = 0.1
a = 1.0
W = 10
sigma = 0.0
N = 100
num_trials = 20
Hs = np.arange(1, 11)

corr_mmse = np.zeros_like(Hs, dtype=float)
corr_svr  = np.zeros_like(Hs, dtype=float)

for idx, H in enumerate(Hs):
  cm = []
  cs = []
  for _ in range(num_trials):
    Xtrain, ytrain_clean, _, _, Xtest, ytest_clean = build_datasets(N, W, H, sigma)
    ytrain = add_outliers(ytrain_clean.copy(), p, a)

    X = Xtrain.T
    y = ytrain
    XXT = X @ X.T
    Xy  = X @ y

    w_mmse = np.linalg.solve(XXT, Xy)
    yhat_mmse = Xtest @ w_mmse

    svr = SVR(C=10.0, epsilon=0.01, kernel='linear')
    svr.fit(Xtrain, ytrain)
    yhat_svr = svr.predict(Xtest)

    cm.append(corrcoef(ytest_clean, yhat_mmse))
    cs.append(corrcoef(ytest_clean, yhat_svr))

  corr_mmse[idx] = np.mean(cm)
  corr_svr[idx]  = np.mean(cs)

plt.figure()
plt.plot(Hs, corr_mmse, 'o-', label='MMSE')
plt.plot(Hs, corr_svr,  'o-', label='SVR')
plt.xlabel('Prediction horizon H')
plt.ylabel('Correlation coefficient')
plt.title('Correlation vs Horizon (averaged)')
plt.grid()
plt.legend()
plt.tight_layout()
plt.savefig("plots/corr_vs_horizon.png")
plt.close()
print("5.2A: saved plots/corr_vs_horizon.png")

# 5.2B
print("\n5.2B")

def get_run(H):
  Xtrain, ytrain_clean, _, _, Xtest, ytest_clean = build_datasets(N, W, H, sigma)
  ytrain = add_outliers(ytrain_clean.copy(), p, a)

  X = Xtrain.T
  y = ytrain
  XXT = X @ X.T
  Xy  = X @ y
  w_mmse = np.linalg.solve(XXT, Xy)
  yhat_mmse = Xtest @ w_mmse

  svr = SVR(C=10.0, epsilon=0.01, kernel='linear')
  svr.fit(Xtrain, ytrain)
  yhat_svr = svr.predict(Xtest)

  return ytest_clean, yhat_mmse, yhat_svr

Hs_scatter = [1, 5]

plt.figure(figsize=(10,8))
for i, H in enumerate(Hs_scatter, 1):
  ytrue, yhat_m, yhat_s = get_run(H)
  ymin = min(ytrue.min(), yhat_m.min(), yhat_s.min())
  ymax = max(ytrue.max(), yhat_m.max(), yhat_s.max())
  line = np.linspace(ymin, ymax, 100)

  plt.subplot(2,2,2*i-1)
  plt.scatter(ytrue, yhat_m, s=10, alpha=0.7)
  plt.plot(line, line, 'r--')
  plt.xlabel('True')
  plt.ylabel('MMSE prediction')
  plt.title(f'H = {H} (MMSE)')
  plt.grid()

  plt.subplot(2,2,2*i)
  plt.scatter(ytrue, yhat_s, s=10, alpha=0.7)
  plt.plot(line, line, 'r--')
  plt.xlabel('True')
  plt.ylabel('SVR prediction')
  plt.title(f'H = {H} (SVR)')
  plt.grid()

plt.tight_layout()
plt.savefig("plots/corr_scatter_H1_H5.png")
plt.close()
print("5.2B: saved plots/corr_scatter_H1_H5.png")

# 5.2D-F
print("\n5.2D-F")

def svr_grid_search(H, epsilons, cs):
  Xtrain, ytrain_clean, Xval, yval_clean, Xtest, ytest_clean = build_datasets(N, W, H, sigma)
  ytrain = add_outliers(ytrain_clean.copy(), p, a)

  corr_grid = np.zeros((len(epsilons), len(cs)))
  for i, eps in enumerate(epsilons):
    for j, C in enumerate(cs):
      svr = SVR(C=C, epsilon=eps, kernel='linear')
      svr.fit(Xtrain, ytrain)
      yhat_val = svr.predict(Xval)
      corr_grid[i, j] = corrcoef(yval_clean, yhat_val)

  best_idx = np.unravel_index(np.argmax(corr_grid), corr_grid.shape)
  best_eps = epsilons[best_idx[0]]
  best_C   = cs[best_idx[1]]

  svr_best = SVR(C=best_C, epsilon=best_eps, kernel='linear')
  svr_best.fit(Xtrain, ytrain)
  yhat_test = svr_best.predict(Xtest)
  corr_test = corrcoef(ytest_clean, yhat_test)

  return Xtrain, ytrain_clean, Xval, yval_clean, Xtest, ytest_clean, corr_grid, best_eps, best_C, yhat_test, corr_test

max_epsilons = 10
max_cs = 10
epsilons = np.logspace(-2, 0, max_epsilons)
cs       = np.logspace(-2, 2, max_cs)

for H in [1, 2]:
  print(f"\nH = {H}")
  (Xtrain, ytrain_clean, Xval, yval_clean,
   Xtest, ytest_clean, corr_grid, best_eps, best_C,
   yhat_test, corr_test) = svr_grid_search(H, epsilons, cs)

  print(f"Best epsilon: {best_eps:.4f}, Best C: {best_C:.4f}, Test corr: {corr_test:.4f}")

  Xg, Yg = np.meshgrid(np.log10(cs), np.log10(epsilons))
  plt.figure()
  cp = plt.contourf(Xg, Yg, corr_grid, levels=20)
  plt.colorbar(cp, label='Correlation')
  plt.xlabel(r'$\log_{10}(C)$')
  plt.ylabel(r'$\log_{10}(\varepsilon)$')
  plt.title(f'SVR correlation on validation set (H={H})')
  plt.tight_layout()
  plt.savefig(f"plots/contour_corr_H{H}.png")
  plt.close()
  print(f"Saved plots/contour_corr_H{H}.png")

  t = np.arange(len(ytest_clean))
  plt.figure()
  plt.plot(t, ytest_clean, label='Clean test', linewidth=2)
  plt.plot(t, yhat_test, '--', label='SVR (best)')
  plt.xlabel('Time')
  plt.ylabel('Signal')
  plt.title(f'Test prediction with best SVR hyperparameters (H={H})')
  plt.grid()
  plt.legend()
  plt.tight_layout()
  plt.savefig(f"plots/pred_best_SVR_H{H}.png")
  plt.close()
  print(f"Saved plots/pred_best_SVR_H{H}.png")
