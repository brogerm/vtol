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
F_max = P.F_max
tau_max = P.tau_max

# tuning parameters
wn_h = 2
wn_th = 1
wn_z = 3
zeta_h = 0.707
zeta_th = 0.707
zeta_z = 0.707

# Longitudinal Parameters
Alon = np.matrix([[0.0, 1.0],
                [0.0, 0.0]])

Blon = np.matrix([[0],
                [1/(P.mc + 2*P.mr)]])

Clon = np.matrix([[1.0, 0.0]])

# gain calculation
des_char_poly =[1, 2*zeta_h*wn_h, wn_h**2]
des_poles = np.roots(des_char_poly)

# Compute the gains if the system is controllable
if np.linalg.matrix_rank(cnt.ctrb(Alon, Blon)) != 2:
    print("The system is not controllable")
else:
    K_lon = cnt.acker(Alon, Blon, des_poles)
    kr_lon = -1.0/(Clon[0]*np.linalg.inv(Alon-Blon*K_lon)*Blon)

# Lateral parameters

# State Space Equations
# xdot = A*x + B*u
# y = C*x
Fe = P.m * P.g
Alat = np.matrix([[0.0, 0.0,                 1.0,                   0.0],
               [0.0, 0.0,                 0.0,                   1.0],
               [0.0, -Fe/(P.mc + 2*P.mr), -P.mu/(P.mc + 2*P.mr), 0.0],
               [0.0, 0.0,                 0.0,                   0.0]])

Blat = np.matrix([[0.0],
               [0.0],
               [0.0],
               [1/(P.Jc + 2*P.mr*P.d**2)]])

Clat = np.matrix([[1.0, 0.0, 0.0, 0.0],
               [0.0, 1.0, 0.0, 0.0]])

# gain calculation
des_char_poly = np.convolve([1, 2*zeta_th*wn_th, wn_th**2],
                            [1, 2*zeta_z*wn_z, wn_z**2])
des_poles = np.roots(des_char_poly)

# Compute the gains if the system is controllable
if np.linalg.matrix_rank(cnt.ctrb(Alat, Blat)) != 4:
    print("The system is not controllable")
else:
    K_lat = cnt.acker(Alat, Blat, des_poles)
    kr_lat = -1.0/(Clat[0]*np.linalg.inv(Alat-Blat*K_lat)*Blat)

print('K_lon: ', K_lon)
print('kr_lon: ', kr_lon)
print('K_lat: ', K_lat)
print('kr_lat: ', kr_lat)



