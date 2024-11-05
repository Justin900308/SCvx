from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp
from numpy import linalg as LA

# Define system matrices
A = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [0, 0, 0, 0], [0, 0, 0, 0]])
B = np.array([[0, 0], [0, 0], [1, 0], [0, 1]])
C = np.eye(4)
D = np.zeros((4, 2))

# Continuous-time system
sys = signal.StateSpace(A, B, C, D)

Ts = 0.1  # Sampling time
T_end = 8  # End time

T = T_end / Ts + 1

# Discretize the system
sysd = sys.to_discrete(Ts)
Ad = sysd.A
Bd = sysd.B

# Initialize state and input arrays
count = 0
X = np.zeros((4, 1))  # State trajectory
X[:, 0] = np.array([0, 0, 0, 0])  # Initial state

u = np.empty((2, 0))  # Initialize u as an empty 2x0 array to store control inputs

# Compute the initial trajectory
for t in np.arange(0, T_end + Ts, Ts):
    if t < 2:
        u = np.concatenate((u, np.array([[0.1], [0.1]])), axis=1)  # Append control input [0, 1]
    else:
        u = np.concatenate((u, np.array([[0.1], [0.1]])), axis=1)  # Append control input [1, 1]

    # Compute the next state
    X = np.hstack((X, Ad @ X[:, count].reshape(-1, 1) + Bd @ u[:, count].reshape(-1, 1)))
    count += 1

# Reshape the state and input vectors
XX = X[:, :-1].reshape(-1, 1)
uu = np.array(u).reshape(-1, 1)


# Function to stack state and control vectors (equivalent to stack_vec in MATLAB)
def stack_vec(X, u):
    return np.vstack((X, u))


z = stack_vec(XX, uu)

# Plot the trajectory
plt.plot(X[0, :], X[1, :], '.')
plt.xlabel('X1')
plt.ylabel('X2')
plt.title('Trajectory Plot')
plt.grid(True)
# plt.show()


from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp

# Initialize variables

# obs_center = np.array([4, 6])
obs_center = np.zeros((int(T), 2))
for i in range(int(T)):
    obs_center[i, :] = np.array([12 - 0.02 * i, 8 - 0.00 * i])

R = 2
alpha = 1
r_default = 1

lambda_param = 10000
rho0 = 0.01
rho1 = 0.2
rho2 = 0.9
tol = 0.001

# Time length N and trajectory state X
N = X.shape[1]

# Convert u to a numpy array for proper matrix operations
u = np.array(u)
u = np.reshape(u, (2, int(T_end / Ts + 1)))

# Plot the obstacle (circle)
theta = np.linspace(0, 2 * np.pi, 201)
x_theta = R * np.cos(theta)
y_theta = R * np.sin(theta)
plt.figure(1)

for i in range(int(T)):
    plt.plot(X[0, :], X[1, :], '.')
    plt.plot(obs_center[i,0] + x_theta, obs_center[i,1]  + y_theta)



# Start the iterative optimization process
for k in range(201):

    # Define variables for optimization
    w = cp.Variable((2, N - 1))
    v = cp.Variable((4, N - 1))
    d = cp.Variable((4, N))

    s = cp.Variable((1, N - 1))

    # Define the cost function
    Linear_cost = 1 * cp.norm(((u + w) * Ts), 1) + 1 * lambda_param * cp.sum(
        cp.sum(cp.abs(v))) + 1 * lambda_param * cp.sum(cp.pos(s))

    # cost = 500 * cp.sum(U * Ts) + 1 * lambda_param * cp.sum(cp.sum(cp.abs(v)))
    # 500 * np.sum(U_val)
    # 1 * lambda_param * np.sum(np.sum(np.abs(v_val)))
    # Define constraints
    constraints = [d[:, 0] == np.zeros(4)]

    E = np.eye(4)

    for i in range(N - 1):
        constraints.append(
            X[:, i + 1] + d[:, i + 1] == (Ad @ X[:, i] + Ad @ d[:, i]) + (Bd @ u[:, i] + Bd @ w[:, i]) + E @ v[:, i])

        # constraints.append(cp.abs(w[0, i]) <= r_default)
        constraints.append(w[0, i] <= r_default)
        constraints.append(-r_default <= w[0, i])
        constraints.append(w[1, i] <= r_default)
        constraints.append(-r_default <= w[1, i])
        # if i<N-2:
        #     constraints.append(cp.abs(v[1:4,i]) <= 0.0000001)
        # constraints.append( cp.abs(v[i]) <= w[1, i])
        # constraints.append(cp.abs(w[1, i]) <= r_default)

        # Obstacle avoidance constraint
        constraints.append(2 * R - cp.norm(X[0:2, i] - obs_center[i, :], 2) - (X[0:2, i] - obs_center[i, :]).T @ (
                X[0:2, i] + d[0:2, i] - obs_center[i, :]) / cp.norm(X[0:2, i] - obs_center[i, :], 2) <= s[:, i])
        constraints.append(s[:, i] >= 0)
        # constraints.append(R - cp.norm(X[0:2, i] - obs_center, 2) - (X[0:2, i] - obs_center).T @ (X[0:2, i] - obs_center) / cp.norm(X[0:2, i] - obs_center, 2) == s[i])

    # Terminal condition
    constraints.append(X[:, N - 1] + d[:, N - 1] == np.array([10, 10, 0, 0]))

    # Define the problem
    problem = cp.Problem(cp.Minimize(Linear_cost), constraints)

    # Solve the optimization problem
    problem.solve(solver=cp.CLARABEL)

    # Update the variables after solving
    w_val = w.value
    v_val = v.value
    d_val = d.value
    s_val = s.value
    # U_val = U.value
    j = 1
    linear_cost = np.zeros((201, 1))
    linear_cost[k, 0] = (1 * LA.norm(((u + w_val) * Ts), ord=1) +
                         1 * lambda_param * np.sum(np.sum(np.abs(v_val))) +
                         1 * lambda_param * np.sum(s_val))

    rho0 = 0
    rho1 = 0.25
    rho2 = 0.7
    delta_L = (linear_cost[k + 1, 0] - linear_cost[k, 0]) / linear_cost[k, 0]
    if delta_L <= rho0:
        r_default = np.max((r_default, 0.8))
        X = X + d_val
        u = u + w_val
    elif delta_L <= rho1:
        r_default = r_default
        X = X + d_val
        u = u + w_val
    elif delta_L <= rho2:
        r_default = r_default / 3.2
        X = X + d_val
        u = u + w_val
    else:
        X = X
        u = u
    r_default = 1
    # X[:, j + 1] + d_val[:, j + 1]
    # (Ad @ X[:, j] + Ad @ d_val[:, j]) + (Bd @ u[:, j] + Bd @ w_val[:, j]) + E @ v_val[:, j]
    # Update the trajectory

    # Plot the updated trajectory
    plt.plot(X[0, :], X[1, :], '.')
    plt.pause(0.001)
    # plt.clf()
    ss = np.zeros((N - 1))

    ss_max = np.array([0])
    for i in range(N - 1):
        ss[i] = R - LA.norm(X[0:2, i] - obs_center, 2)

    if (np.max(ss) < 0) and (k > 20):
        break

    # # Break if the max value of s is less than 0 (constraint satisfied)
    # if np.max(s.value) < 0:
    #     break

# Final trajectory plot
plt.clf()

plt.figure(2)
for i in range(int(T)):

    plt.plot(X[0, i], X[1, i], '.')
    plt.plot(obs_center[i,0] + x_theta, obs_center[i,1]  + y_theta)
    plt.xlim((0,15))
    plt.ylim((0,15))
    plt.pause(0.001)
    plt.clf()


# plt.figure(2)
# plt.plot(X[0, :], X[1, :], '.')
# plt.plot(obs_center[0] + x_theta, obs_center[1] + y_theta)
# plt.show()
