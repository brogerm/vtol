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
tr_h = 0.15
tr_th = 0.5
tr_z = 0.9
wn_h = 2.2/tr_h
wn_th = 2.2/tr_th
wn_z = 2.2/tr_z
zeta_h = 0.707
zeta_th = 0.707
zeta_z = 0.707
h_integrator_pole = -1
z_integrator_pole = -2

# Longitudinal Parameters
Alon = np.matrix([[0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [-1.0, 0.0, 0.0]])

Blon = np.matrix([[0],
                [1/(P.mc + 2*P.mr)],
                  [0.0]])

Clon = np.matrix([[1.0, 0.0]])

# gain calculation
des_char_poly = np.convolve([1, 2*zeta_h*wn_h, wn_h**2],
    np.poly([h_integrator_pole]))
des_poles = np.roots(des_char_poly)

# Compute the gains if the system is controllable
if np.linalg.matrix_rank(cnt.ctrb(Alon, Blon)) != 3:
    print("The system is not controllable")
else:
    K1_lon = cnt.acker(Alon, Blon, des_poles)
    K_lon = np.matrix([K1_lon.item(0), K1_lon.item(1)])
    ki_lon = K1_lon.item(2)

# Lateral parameters

# State Space Equations
# xdot = A*x + B*u
# y = C*x
Fe = P.m * P.g
Alat = np.matrix([[0.0, 0.0,                 1.0,                   0.0, 0.0],
                [0.0, 0.0,                 0.0,                   1.0, 0.0],
                [0.0, -Fe/(P.mc + 2*P.mr), -P.mu/(P.mc + 2*P.mr), 0.0, 0.0],
                [0.0, 0.0,                 0.0,                   0.0, 0.0],
                [-1.0, 0.0, 0.0, 0.0, 0.0]])

Blat = np.matrix([[0.0],
                [0.0],
                [0.0],
                [1/(P.Jc + 2*P.mr*P.d**2)],
                [0.0]])

Clat = np.matrix([[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0]])

# gain calculation
des_char_poly = np.convolve(
    np.convolve([1, 2*zeta_th*wn_th, wn_th**2],
                [1, 2*zeta_z*wn_z, wn_z**2]),
    np.poly([z_integrator_pole]))
des_poles = np.roots(des_char_poly)

# Compute the gains if the system is controllable
if np.linalg.matrix_rank(cnt.ctrb(Alat, Blat)) != 5:
    print("The system is not controllable")
else:
    K1_lat = cnt.acker(Alat, Blat, des_poles)
    K_lat = np.matrix([K1_lat.item(0), K1_lat.item(1), K1_lat.item(2), K1_lat.item(3)])
    ki_lat = K1_lat.item(4)

print('K1_lon', K1_lon)
print('K_lon: ', K_lon)
print('ki_lon: ', ki_lon)
print('K1_lat', K1_lat)
print('K_lat: ', K_lat)
print('ki_lat: ', ki_lat)



