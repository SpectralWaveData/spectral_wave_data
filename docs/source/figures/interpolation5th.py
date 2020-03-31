import numpy as np
import matplotlib.pyplot as plt

# data points
dt  = 1.0

fim1 = 1.0
dfim1 = -1.3
fi = 2.0
dfi = -0.5
fip1 = 0.7
dfip1 = -0.3
fip2 = 0.7
dfip2 = 3.3

# Final curve between [-dt, dt]
q1 = dfi * dt
q2 = fim1 - 2 * fi + fip1  + (dfim1 - dfip1) * dt/4
q3 = -3 * fim1 - 3 * fi + 5 * fip1 + fip2 - (3 * dfim1  + 23 * dfi + 13 * dfip1 + dfip2) * dt / 4
q4 = 3 * fim1 + 7 * fi - 8 * fip1 - 2 * fip2 + (3 * dfim1  + 30 * dfi + 25 * dfip1 + 2 * dfip2) * dt / 4
q5 = - fim1 - 3 * fi + 3 * fip1 + fip2 - (dfim1  + 11 * dfi + 11 * dfip1 + dfip2) * dt / 4


# curve defining ddfi
a0 = fi
a1 = dfi
a2 = -(8 * fi - 4 * fim1 - 4 * fip1 - dfim1 * dt + dfip1 * dt)/(4 * dt)
a3 = -(5 * fim1 - 5 * fip1 + 8 * dfi * dt + dfim1 * dt + dfip1 * dt) / (4 * dt)
a4 = -(-4 * fi + 2 * fim1 + 2 * fip1 + dfim1 * dt - dfip1 * dt)/ ( 4 * dt)
a5 = -(-3 * fim1 + 3 * fip1 - 4 * dfi * dt - dfim1 * dt - dfip1 * dt)/ (4 * dt)

# curve defining ddfip1
b0 = fip1
b1 = dfip1
b2 = -(8 * fip1 - 4 * fi - 4 * fip2 - dfi * dt + dfip2 * dt)/(4 * dt)
b3 = -(5 * fi - 5 * fip2 + 8 * dfip1 * dt + dfi * dt + dfip2 * dt) / (4 * dt)
b4 = -(-4 * fip1 + 2 * fi + 2 * fip2 + dfi * dt - dfip2 * dt)/ ( 4 * dt)
b5 = -(-3 * fi + 3 * fip2 - 4 * dfip1 * dt - dfi * dt - dfip2 * dt)/ (4 * dt)


def fun(t):
    d = (t - 0.0) / dt
    return fi + q1 * d  + q2 * d**2 + q3 * d**3 + q4 * d**4 + q5 * d**5
t = np.linspace(0.0, dt, 500)
f = fun(t)

def fun_a(t):
    return a0 + a1 * t + a2 * t ** 2 + a3 * t ** 3 + a4 * t ** 4 + a5 * t ** 5
ta = np.linspace(-dt, dt, 500)
fa = fun_a(ta)

def fun_b(t):
    T = t - dt
    return b0 + b1 * T + b2 * T ** 2 + b3 * T ** 3 + b4 * T ** 4 + b5 * T ** 5
tb = np.linspace(0.0, 2*dt, 500)
fb = fun_b(tb)

def point(pos):
    if pos=='xim1':
        x = -dt
        y = fim1
    if pos=='xi':
        x = 0.0
        y = fi
    if pos=='xip1':
        x = dt
        y = fip1
    if pos=='xip2':
        x = 2*dt
        y = fip2
    return (x, y)


def slope(pos):
    ds = 0.3 * dt
    if pos=='xim1':
        x1 = -dt - ds
        x2 = -dt + ds
        y1 = fim1 - dfim1 * ds
        y2 = fim1 + dfim1 * ds
    if pos=='xi':
        x1 = 0.0 - ds
        x2 = 0.0 + ds
        y1 = fi - dfi * ds
        y2 = fi + dfi * ds
    if pos=='xip1':
        x1 = dt - ds
        x2 = dt + ds
        y1 = fip1 - dfip1 * ds
        y2 = fip1 + dfip1 * ds
    if pos=='xip2':
        x1 = 2*dt - ds
        x2 = 2*dt + ds
        y1 = fip2 - dfip2 * ds
        y2 = fip2 + dfip2 * ds
    return ((x1, x2), (y1, y2))


p1 = point('xim1')
p2 = point('xi')
p3 = point('xip1')
p4 = point('xip2')
s1 = slope('xim1')
s2 = slope('xi')
s3 = slope('xip1')
s4 = slope('xip2')

fig, ax = plt.subplots()
ax.plot(t, f, label='5th order spline $f(t), t\in[t_i, t_{i+1}]$', lw=3, color='black')
ax.plot(ta, fa, label=r'5th order spline defining $d^2f_i / dt^2$', color='green')
ax.plot(tb, fb, label=r'5th order spline defining $d^2f_{i+1} / dt^2$', color='magenta')
ax.plot(p1[0], p1[1], marker='o',  markersize=10, color='red', label=r'Specified values $f_i$')
ax.plot(p2[0], p2[1], marker='o',  markersize=10, color='red')
ax.plot(p3[0], p3[1], marker='o',  markersize=10, color='red')
ax.plot(p4[0], p4[1], marker='o',  markersize=10, color='red')
ax.plot(s1[0], s1[1], color='blue', label=r'Specified slopes $df_i/dt$')
ax.plot(s2[0], s2[1], color='blue')
ax.plot(s3[0], s3[1], color='blue')
ax.plot(s4[0], s4[1], color='blue')
ax.legend()
#ax.set_aspect('equal')
plt.xticks([-1, 0, 1, 2], [r'$t_{i-1}$', r'$t_{i}$', r'$t_{i+1}$', r'$t_{i+2}$'])
ax.set_xlabel('Time')
ax.set_ylabel('f-values')
plt.show()

