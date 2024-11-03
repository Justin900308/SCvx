from array import array

import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt
from numpy import linalg as LA
from scipy import signal

# A=np.random.rand(4,4)
# B=np.random.rand(4,2)
######################################################################
n = 2
m = 2
T = 100
N = m * (T - 1) + n * T
Ncon = n * (T - 1) + T
Ts = 0.1
A = np.array([[0, 0], [0, 0]])
# A=A+np.eye(n)
B = np.array([[1, 0], [0, 1]])
C = np.eye(2)
D = np.zeros((2, 2))
sys = signal.StateSpace(A, B, C, D)
sysd = sys.to_discrete(Ts)
Ad = sysd.A
Bd = sysd.B

C = B
n = A.shape[0]

#######################################################################
u = 1 * np.ones([2, 99])
u[0, 0:50] = 10
u[0, 50:99] = 10
u[1, 50:99] = -0.00000000
x = np.zeros([m, 100])

x[:, 0] = [20, 20]  # comes out

u = -1.0 * u

for i in range(len(u[0])):
    x[:, i + 1] = (Ad @ x[:, i]) + (Bd @ u[:, i])

# print(x[:,-1])
plt.scatter(x[0, :], x[1, :])
# plt.show()

#######################################################################
centre1 = np.zeros((T, 2))
a = 5
for i in range(T):
    centre1[i, :] = np.array((10 + i * 0., 10 + i * 0.))
    t = np.linspace(0, 2 * np.pi, 100)
    xe = centre1[i, 0] + a * np.cos(t)
    ye = centre1[i, 1] + a * np.sin(t)

    plt.plot(xe, ye)
    plt.pause(0.001)

centre = [10.0, 10.0]
cen = np.array([10, 10, 0, 0])

t = np.linspace(0, 2 * np.pi, 100)
xe = centre[0] + a * np.cos(t)
ye = centre[1] + a * np.sin(t)

plt.plot(xe, ye)

#######################################################################

x_stacked = x.T.reshape(-1)
u_stacked = u.T.reshape(-1)
z = np.hstack((x_stacked, u_stacked))  # initial z
z = np.reshape(z, (N, 1))


# plt.gca().set_aspect('equal')
# plt.show()

#######################################################################

def z_bar(z, j):
    zbar = cp.Variable((N, 1))
    x = zbar[:n * T]
    u = zbar[n * T:]

    objective = cp.Minimize(cp.norm(zbar - z, 2))
    if j < n * (T - 1):
        t = int(j / n)
        x_t = x[n * t:n * (t + 1)]
        x_tp1 = x[n * (t + 1):n * (t + 2)]
        u_t = u[m * t:m * (t + 1)]
        g = Ad @ x_t + Bd @ u_t - x_tp1
        const = [g[np.mod(j, n)] <= 0]

    else:
        t = j - n * (T - 1)
        x = zbar[:n * T]
        x_t = x[n * t:n * (t + 1)]
        h = (x_t[0] - centre1[t, 0]) ** 2 + (x_t[1] - centre1[t, 1]) ** 2 - a ** 2
        const = [h <= 0]

    prb = cp.Problem(objective, const)
    prb.solve(solver=cp.CLARABEL)
    # print(prb.value)
    return zbar.value


def z2xu(z, n, m, T):
    x_st = z[:n * T]
    u_st = z[n * T:]

    X = x_st.reshape((T, n)).T
    U = u_st.reshape((T - 1, m)).T
    return X, U


def gfun(y, Ad, Bd, T, n, m):
    x = y[:n * T]
    u = y[n * T:]

    out = 0
    for t in range(T - 1):
        x_t = x[n * t:n * (t + 1)]
        x_tp1 = x[n * (t + 1):n * (t + 2)]
        u_t = u[m * t:m * (t + 1)]
        gi = Ad @ x_t + Bd @ u_t - x_tp1
        out += cp.norm(gi, 1)
    return out


def grad_g_t(y, t):
    out = np.zeros((n, N))
    out[:, n * t:n * (t + 1)] = Ad
    out[:, n * (t + 1):n * (t + 2)] = -np.eye(n)
    out[:, n * T + m * t: n * T + m * (t + 1)] = Bd
    return out


def grad_h_t(y, t):
    out = np.zeros((N, 1))
    x = y[:n * T]
    x_t = x[n * t:n * t + 2]
    out[n * t:n * t + 2] = 0.5 * (x_t - centre1[t, 0]) / LA.norm(x_t - centre1[t, 1])
    return out


def grad_q_j(y, j):
    if j < n * (T - 1):
        t = int(j / n)
        out = grad_g_t(y, t)
        out = out[np.mod(j, n), :]
    else:
        t = j - n * (T - 1)
        out = grad_h_t(y, t)
    return out

velocity = np.zeros((m * (T - 1), 1))
Xinit, Unit = z2xu(z, n, m, T)
print(x[0, 0], x[1, 0])
for it in range(20):
    y2 = cp.Variable((N, 1))

    lbda = 100
    gn = gfun(y2, Ad, Bd, T, n, m)


    # v_diff=cp.Variable(m*(T-1))
    v_diff = cp.Variable(m * (T - 1))
    # for i in range(m*(T-1)):
    #     velocity[i,0]=2

    for i in range((T - 1)):
        X=z[i*m,0]
        Y=z[i*m+1,0]

        P_error = (np.array(([-25], [-25])) - ([z[i*m,0]],[z[i*m+1,0]]))
        velocity[i*n:i*n+2, 0] = 0.1*P_error.flatten()

    v_diff = cp.norm(y2[n * T:] - velocity, 2)

    Objective2 = cp.Minimize(v_diff + lbda * gn + 0 * cp.norm(y2[n * T - 4:n * T - 2], 2))
    constraint2 = [
        y2[0] == x[0, 0],
        y2[1] == x[1, 0],
        # y2[n * T - 2] == 0,
        # y2[n * T - 1] == 0,

    ]

    # n * (T - 1)
    # Ncon
    for j in range(Ncon):
        # print(j)
        z_bar_j = z_bar(z, j)
        # print(z_bar_j[0])
        dq_j = grad_q_j(z_bar_j, j)
        constraint2.append(np.transpose(dq_j) @ (y2 - z_bar_j) >= 0)

    problem2 = cp.Problem(Objective2, constraint2)
    problem2.solve(solver=cp.CLARABEL)
    z_new = y2.value
    print(np.linalg.norm(z_new - z, 2))
    if np.linalg.norm(z_new - z, 2) < 1:
        print("Converged")
        break
    z = z_new
    # print("z:",z)
    print("it:", it)
    Xn, Un = z2xu(z_new, n, m, T)
    plt.plot(xe, ye, color='black')
    plt.scatter(Xn[0, :], Xn[1, :])
    plt.pause(0.001)
    plt.clf()

for i in range(T):
    centre1[i, :] = np.array((10 + i * 0., 10 + i * 0.))
    t = np.linspace(0, 2 * np.pi, 100)
    xe = centre1[i, 0] + a * np.cos(t)
    ye = centre1[i, 1] + a * np.sin(t)

    plt.plot(xe, ye)
    plt.scatter(Xn[0, 1:i], Xn[1, 1:i])
    plt.pause(0.001)
    plt.clf()

aaa = 5
