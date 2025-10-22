import numpy as np

# ma = F - kx - bv
# a = (F - kx - bv)
# x1 = dx/dt
# x2 = dv/dt

def spm(t, X):
    x = X[0]
    v = X[1]
    dxdt = X[1]
    dvdt = (F - k*x - b*v)
    dxdt = X[1]
    X_dot = np.array([dvdt, dxdt])
    return X_dot

x0 = 0
v0 = 0
X0 = np.array([x0, v0])

F = 4
k = 2
b = 2

t0 = 0
tf = 10
dt = 0.01

t_vec = np.arange(t0, tf, dt)

X = X0

x_vec = [x0]
v_vec = [v0]

for t in t_vec:
    X = X + spm(t, X)*dt
    x_vec.append(X[0])
    v_vec.append(X[1])

print(x_vec[-1])
print(v_vec[-1])