# satellite Parameter File
import numpy as np
import control as cnt
import sys
from scipy import signal
sys.path.append('..')  # add parent directory
import vtolParam as P

# import variables from satelliteParam for later import through current file
Ts = P.Ts
sigma = P.sigma
beta = P.beta
u_max = P.F_max

# tuning parameters
tr_h = 0.15
wn_h = 2.2/tr_h
zeta_h = 0.707
h_integrator_pole = -3
dist_obsv_pole = -1.0         # pole for disturbance observer

# pick observer poles
tr_h_obs = tr_h/5.0
wn_h_obs = 2.2/tr_h_obs

# Longitudinal Parameters
A = np.matrix([[0.0, 1.0],
                [0.0, 0.0]])

B = np.matrix([[0],
                [1/(P.mc + 2*P.mr)]])

C = np.matrix([[1.0, 0.0]])

# Augmented Longitudinal Parameters
A1 = np.matrix([[0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0]])

B1 = np.matrix([[0],
                [1/(P.mc + 2*P.mr)],
                  [0.0]])

# gain calculation
des_char_poly = np.convolve([1, 2*zeta_h*wn_h, wn_h**2],
    np.poly([h_integrator_pole]))
des_poles = np.roots(des_char_poly)

# Compute the gains if the system is controllable
if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 3:
    print("The longitudinal system is not controllable")
else:
    K1 = cnt.acker(A1, B1, des_poles)
    K = np.matrix([K1.item(0), K1.item(1)])
    ki = K1.item(2)

# ----------------------------------------- OBSERVER GAINS -----------------------------------------
# compute longitudinal observer gains
des_obs_char_poly = [1, 2*zeta_h*wn_h_obs, wn_h_obs**2]
des_obs_poles = np.roots(des_obs_char_poly)

# Compute the gains if the system is observable
if np.linalg.matrix_rank(cnt.ctrb(A.T, C.T)) != 2:
    print("The longitudinal observer system is not observable")
else:
    # place_poles returns an object with various properties.  The gains are accessed through .gain_matrix
    # .T transposes the matrix
    L = signal.place_poles(A.T, C.T, des_obs_poles).gain_matrix.T

print('----- Longitudinal -----')
print('K1', K1)
print('K: ', K)
print('ki: ', ki)
print('L^T: ', L.T)


# computer observer gains
# Augmented Matrices
A2 = np.concatenate((
        np.concatenate((A, B), axis=1),
        np.zeros((1, 3))),
        axis=0)
C2 = np.concatenate((C, np.zeros((1, 1))), axis=1)


des_obs_char_poly = np.convolve([1, 2*zeta_h*wn_h_obs, wn_h_obs**2],
    np.poly([dist_obsv_pole]))
des_obs_poles = np.roots(des_obs_char_poly)

# Compute the gains if the system is observable
if np.linalg.matrix_rank(cnt.ctrb(A2.T, C2.T)) != 3:
    print("The disturbance observer system is not observable")
else:
    # place_poles returns an object with various properties.  The gains are accessed through .gain_matrix
    # .T transposes the matrix
    L2 = signal.place_poles(A2.T, C2.T, des_obs_poles).gain_matrix.T
    L = L2[0:2]
    Ld = L2[2:3]

print('L2', L2)
print('L: ', L)
print('Ld: ', Ld)



