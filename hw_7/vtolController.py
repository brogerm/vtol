import numpy as np
import vtolParam as P

class vtolController:
    def __init__(self):
        self.kp_z = P.kp_z
        self.kd_z = P.kd_z
        self.kp_th = P.kp_th
        self.kd_th = P.kd_th
        self.kp_h = P.kp_h
        self.kd_h = P.kd_h
        self.DC_gain = P.DC_gain

    def update(self, z_r, h_r, state):
        z = state[0]
        h = state[1]
        theta = state[2]
        zdot = state[3]
        hdot = state[4]
        thetadot = state[5]

        F_tilde = P.kp_h*(h_r-h) - P.kd_h*hdot
        # the reference angle for theta comes from the outer loop PD control
        theta_r = (self.kp_z * (z_r - z) - self.kd_z * zdot) * self.DC_gain
        # The force applied to the rod comes from the inner loop PD control
        tau = self.kp_th * (theta_r - theta) - self.kd_th * thetadot
        return [F_tilde, tau]
