# Inverted Pendulum Parameter File
import numpy as np
import control as cnt

# Physical parameters of the vtol
ml = 0.25       # Mass of left motor, kg
mc = 1.0        # Mass of center pod, kg
mr = ml         # Mass of right motor, kg
d = 0.3        # Length between center pod and each motor, m
cl = 0.25       # Length of center pod, m
ch = 0.5        # Height of center pod, m
Jc = 0.0042     # Rotational inertia of center pod, kg m^2
radius = 0.125  # radius of circular rotor

# General parameters
g = 9.81        # Gravity, m/s^2
mu = 0.1        # Viscosity, kg/s

# Initial Conditions
z0 = 0.0                # ,m
theta0 = 0.0*np.pi/180  # ,rads
h0 = 0.0                # ,m
zdot0 = 0.0             # ,m/s
thetadot0 = 0.0         # ,rads/s
hdot0 = 0.0             # ,m/s

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 50.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate

# saturation limits
F_max = 5.0                # Max Force, N

####################################################
#       PD Control: Time Design Strategy
####################################################
tr_th = 0.8
zeta_th = 0.707
tr_z = 10 * tr_th
zeta_z = 0.707

# --------------------------
#           Outer Loop
# --------------------------
# coefficients for desired inner loop
wn_th = 2.2/tr_th
alpha1_th = 2.0 * zeta_th * wn_th
alpha0_th = wn_th**2

kd_th = alpha1_th * (Jc + 2 * mr * d**2)
kp_th = alpha0_th * (Jc + 2 * mr * d**2)
DC_gain = 1

# --------------------------
#           Inner Loop
# --------------------------
# coefficients for desired inner loop
wn_z = 2.2/tr_z
alpha1_z = 2.0 * zeta_z * wn_z
alpha0_z = wn_z**2

m = mc + 2 * mr
F0 = m * g
kd_z = -(alpha1_z * m - mu) / (F0 * DC_gain)
kp_z = -alpha0_z * m / (F0 * DC_gain)

#---------------------------------------------
#               Vertical Loop
#---------------------------------------------
tr_h = 0.8
zeta_h = 0.707
wn_h = 2.2/tr_h

alpha1_h = 2.0 * zeta_h * wn_h
alpha0_h = wn_h**2

kd_h = alpha1_h * (mc + 2 * mr)
kp_h = alpha0_h * (mc + 2 * mr)

# saturation limits
F_max = 100.0  # Max Force, N
tau_max = 100.0  # Max Torque, N-m

# dirty derivative parameters
sigma = 0.05  # cutoff freq for dirty derivative
beta = (2.0*sigma-Ts)/(2.0*sigma+Ts)  # dirty derivative gain

print('DC_gain', DC_gain)
print('kp_th: ', kp_th)
print('kd_th: ', kd_th)
print('kp_z: ', kp_z)
print('kd_z: ', kd_z)
print('kp_h: ', kp_h)
print('kd_h: ', kd_h)

