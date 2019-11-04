import numpy as np
import paramHw11 as P

class Controller:
    # state feedback control using dirty derivatives to estimate zdot and thetadot
    def __init__(self):
        self.z_dot = 0.0
        self.h_dot = 0.0
        self.theta_dot = 0.0
        self.z_d1 = 0.
        self.h_d1 = 0.0
        self.theta_d1 = 0.0
        self.K_lon = P.K_lon                 # state feedback gain
        self.kr_lon = P.kr_lon               # Input gain
        self.K_lat = P.K_lat               # state feedback gain
        self.kr_lat = P.kr_lat               # Input gain
        self.Flimit = P.F_max       # Maxiumum Force
        self.Taulimit = P.tau_max       # Maxiumum torque
        self.beta = P.beta           # dirty derivative gain
        self.Ts = P.Ts               # sample rate of controller

    def u(self, y_r, y):
        # y_r is the referenced input
        # y is the current state
        z_r = y_r[0]
        h_r = y_r[1]
        z = y[0]
        h = y[1]
        theta = y[2]

        # differentiate z and theta
        self.differentiateZ(z)
        self.differentiateTheta(theta)
        self.differentiateH(h)

        # Construct the longitudinal and lateral states
        xlon = np.matrix([[h], [self.h_dot]])
        xlat = np.matrix([[z], [theta], [self.z_dot], [self.theta_dot]])

        # Compute the longitudinal state feedback controller
        f_unsat = -self.K_lon*xlon + self.kr_lon*h_r + P.Fe

        # Compute the lateral state feedback controller
        tau_unsat = -self.K_lat*xlat + self.kr_lat*z_r

        f = self.saturate(f_unsat, self.Flimit)
        tau = self.saturate(tau_unsat, self.Taulimit)
        return [f.item(0), tau.item(0)]

    def differentiateZ(self, z):
        '''
            differentiate z
        '''
        self.z_dot = self.beta*self.z_dot + (1-self.beta)*((z - self.z_d1) / self.Ts)
        self.z_d1 = z

    def differentiateTheta(self, theta):
        '''
            differentiate theta
        '''
        self.theta_dot = self.beta*self.theta_dot + (1-self.beta)*((theta - self.theta_d1) / self.Ts)
        self.theta_d1 = theta

    def differentiateH(self, h):
        '''
            differentiate h
        '''
        self.h_dot = self.beta*self.h_dot + (1-self.beta)*((h - self.h_d1) / self.Ts)
        self.h_d1 = h

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u

