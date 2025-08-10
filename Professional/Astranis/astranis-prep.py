import matplotlib.pyplot as plt

def accel_smd(F, x, v, m, b, k):
    """
    Compute acceleration for spring mass damper system

    F = 
    """
    a = (F - b*v - k*x) / m

    return a


m = 1
b = 2
k = 5
x0 = 0
v0 = 0
F = 4

dt = 0.01
tf = 10

x = x0
v = v0

t_val = []
x_val = []
v_val = []
a_val = []

for i in range(int(tf/dt) + 1):
    t = i*dt
    a = accel_smd(F, x, v, m, b, k)
    v += a*dt
    x += v*dt

    t_val.append(t)
    x_val.append(x)
    v_val.append(v)
    a_val.append(a)

for i in range(len(t_val)):
    print(f"{t_val[i]:.2f}, {x_val[i]:.2f}, {v_val[i]:.2f}, {a_val[i]:.2f}")

plt.figure

plt.subplot(3, 1, 1)
plt.plot(t_val, x_val, label = 'Position x(t)')
plt.ylabel('Position (m)')
plt.grid()

plt.subplot(3, 1, 2)
plt.plot(t_val, v_val, label = 'Position x(t)')
plt.ylabel('Velocity (m/s)')
plt.grid()

plt.subplot(3, 1, 3)
plt.plot(t_val, a_val, label = 'Position x(t)')
plt.ylabel('Acceleration (m/s^2)')
plt.grid()

plt.show()