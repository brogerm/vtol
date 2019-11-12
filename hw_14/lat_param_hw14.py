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
u_max = P.tau_max

# tuning parameters
tr_th = 0.5
tr_z = 0.9
wn_th = 2.2/tr_th
wn_z = 2.2/tr_z
zeta_th = 0.707
zeta_z = 0.707
z_integrator_pole = -3.0
dist_obsv_pole = -1.8         # pole for disturbance observer

# pick observer poles
tr_th_obs = tr_th/5.0
tr_z_obs = tr_z/5.0
wn_th_obs = 2.2/tr_th_obs
wn_z_obs = 2.2/tr_z_obs

# Lateral parameters
Fe = P.m * P.g
A = np.matrix([[0.0, 0.0,                 1.0,                   0.0],
                [0.0, 0.0,                 0.0,                   1.0],
                [0.0, -Fe/(P.mc + 2*P.mr), -P.mu/(P.mc + 2*P.mr), 0.0],
                [0.0, 0.0,                 0.0,                   0.0]])

B = np.matrix([[0.0],
                [0.0],
                [0.0],
                [1/(P.Jc + 2*P.mr*P.d**2)]])

C = np.matrix([[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0]])

# Augmented lateral parameters
A1 = np.matrix([[0.0, 0.0,                 1.0,                   0.0, 0.0],
                [0.0, 0.0,                 0.0,                   1.0, 0.0],
                [0.0, -Fe/(P.mc + 2*P.mr), -P.mu/(P.mc + 2*P.mr), 0.0, 0.0],
                [0.0, 0.0,                 0.0,                   0.0, 0.0],
                [-1.0, 0.0, 0.0, 0.0, 0.0]])

B1 = np.matrix([[0.0],
                [0.0],
                [0.0],
                [1/(P.Jc + 2*P.mr*P.d**2)],
                [0.0]])

C1 = np.matrix([[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0]])

# gain calculation
des_char_poly = np.convolve(
    np.convolve([1, 2*zeta_th*wn_th, wn_th**2],
                [1, 2*zeta_z*wn_z, wn_z**2]),
    np.poly([z_integrator_pole]))
des_poles = np.roots(des_char_poly)

# Compute the gains if the system is controllable
if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 5:
    print("The lateral system is not controllable")
else:
    K1 = cnt.acker(A1, B1, des_poles)
    K = np.matrix([K1.item(0), K1.item(1), K1.item(2), K1.item(3)])
    ki = K1.item(4)

# ----------------------------------------- OBSERVER GAINS -----------------------------------------
# compute lateral observer gains
des_obs_char_poly = np.convolve([1, 2*zeta_z*wn_z_obs, wn_z_obs**2],
                                 [1, 2*zeta_th*wn_th_obs, wn_th_obs**2])
des_obs_poles = np.roots(des_obs_char_poly)

# Compute the gains if the system is observable
if np.linalg.matrix_rank(cnt.ctrb(A.T, C.T)) != 4:
    print("The lateral observer system is not observable")
else:
    # place_poles returns an object with various properties.  The gains are accessed through .gain_matrix
    # .T transposes the matrix
    L = signal.place_poles(A.T, C.T, des_obs_poles).gain_matrix.T


print('----- Lateral-----')
print('K1', K1)
print('K: ', K)
print('ki: ', ki)
print('L^T: ', L.T)

# computer observer gains
# Augmented Matrices
A2 = np.concatenate((
        np.concatenate((A, B), axis=1),
        np.zeros((1, 5))),
        axis=0)
C2 = np.concatenate((C, np.zeros((2, 1))), axis=1)

des_obs_char_poly = np.convolve(
    np.convolve([1, 2*zeta_z*wn_z_obs, wn_z_obs**2], [1, 2*zeta_th*wn_th_obs, wn_th_obs**2]),
    np.poly([dist_obsv_pole]))
des_obs_poles = np.roots(des_obs_char_poly)

# Compute the gains if the system is observable
if np.linalg.matrix_rank(cnt.ctrb(A2.T, C2.T)) != 5:
    print("The system is not observable")
else:
    # place_poles returns an object with various properties.  The gains are accessed through .gain_matrix
    # .T transposes the matrix
    L2 = signal.place_poles(A2.T, C2.T, des_obs_poles).gain_matrix.T
    L = L2[0:4, 0:2]
    Ld = L2[4:5, 0:2]

print('L2', L2)
print('L: ', L)
print('Ld: ', Ld)



