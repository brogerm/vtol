# satellite Parameter File
import numpy as np
import control as cnt
import sys
sys.path.append('..')  # add parent directory
import vtolParam as P

# import variables from satelliteParam for later import through current file
Ts = P.Ts
sigma = P.sigma
beta = P.beta
tau_max = P.tau_max

# tuning parameters
wn_h = 0.6
wn_th = 0.15
wn_z = 1
zeta_h = 0.707
zeta_th = 0.707
zeta_z = 0.707

# Longitudinal Parameters
Al = np.matrix([[0.0, 1.0],
                [0.0, 0.0]])

Bl = np.matrix([[0],
                [1/(P.mc + 2*P.mr)]])

Cl = np.matrix([[1.0, 0.0]])

# gain calculation
des_char_poly =[1, 2*zeta_h*wn_h, wn_h**2]
des_poles = np.roots(des_char_poly)

# Compute the gains if the system is controllable
if np.linalg.matrix_rank(cnt.ctrb(Al, Bl)) != 2:
    print("The system is not controllable")
else:
    K_l = cnt.acker(Al, Bl, des_poles)
    kr_l = -1.0/(Cl[0]*np.linalg.inv(Al-Bl*K_l)*Bl)

# Lateral parameters

# State Space Equations
# xdot = A*x + B*u
# y = C*x
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

# gain calculation
des_char_poly = np.convolve([1, 2*zeta_th*wn_th, wn_th**2],
                            [1, 2*zeta_z*wn_z, wn_z**2])
des_poles = np.roots(des_char_poly)

# Compute the gains if the system is controllable
if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
    print("The system is not controllable")
else:
    K_lat = cnt.acker(A, B, des_poles)
    kr_lat = -1.0/(C[0]*np.linalg.inv(A-B*K_lat)*B)

print('K_l: ', K_l)
print('kr_l: ', kr_l)
print('K_lat: ', K_lat)
print('kr_lat: ', kr_lat)



