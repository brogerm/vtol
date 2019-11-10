import numpy as np
import lon_param_hw13 as P

class LonController:
    # state feedback control using dirty derivatives to estimate zdot and thetadot
    def __init__(self):
        self.x_hat = np.matrix([
            [0.0],  # initial estimate for h_hat
            [0.0]])  # initial estimate for h_hat_dot
        self.f_d1 = 0.0  # Computed torque, delayed by one sample
        self.h_d1 = 0.0
        self.K = P.K  # state feedback gain
        self.ki = P.ki  # Input gain
        self.limit = P.u_max  # Maxiumum Force
        self.h_integrator = 0.0
        self.h_error_d1 = 0.0  # error signal delayed by 1 sample
        self.L = P.L  # observer gain
        self.A = P.A  # system model
        self.B = P.B
        self.C = P.C
        self.beta = P.beta  # dirty derivative gain
        self.Ts = P.Ts  # sample rate of controller

    def u(self, y_r, y):
        # y_r is the referenced input
        # y is the current state
        h_r = y_r[0]

        # update the observer and extract h_hat
        self.updateObserver(y)
        h_hat = self.x_hat[0]

        # integrate error
        h_error = h_r - h_hat
        self.h_integrateError(h_error)

        # Compute the longitudinal state feedback controller
        f_unsat = -self.K * self.x_hat - self.ki * self.h_integrator
        f = self.saturate(f_unsat)
        self.updateForce(f)

        return [f.item(0)]

    def updateObserver(self, y_m):
        N = 10
        y = np.matrix([
            [y_m[0]]])
        for i in range(0, N):
            self.x_hat = self.x_hat + self.Ts/float(N)*(
                self.A*self.x_hat + self.B*self.f_d1 + self.L*(y-self.C*self.x_hat)
            )
            print(self.x_hat)

    def updateForce(self, f):
        self.f_d1 = f

    def h_integrateError(self, error):
        self.h_integrator = self.h_integrator + (self.Ts / 2.0) * (error + self.h_error_d1)
        self.h_error_d1 = error

    def saturate(self, u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u

