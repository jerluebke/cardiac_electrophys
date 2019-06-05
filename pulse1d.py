# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from alpha import AlphaBase


def Lap1d(a, o, dx, bc='neumann'):
    if bc == 'neumann':
        o[1:-1] = (a[2:] - 2*a[1:-1] + a[:-2]) / dx**2
        o[0], o[-1] = o[1], o[-2]
    elif bc == 'periodic':
        o[:] = (np.roll(a, 1) + np.roll(a, -1) - 2*a) / dx**2
    return o


class Pulse1d(AlphaBase):
    def __init__(self, xmax, dx, tmax, eta, plot_interval=1,
                 **alpha_params):
        super().__init__(**alpha_params)
        self.dx     = dx
        self.eta    = eta
        self.dt     = .1 * dx**2 / (2. * eta)
        self.t      = np.arange(0, tmax, self.dt)
        self.x      = np.arange(0, xmax, dx)
        self.steps  = self.t.size
        self.V      = np.zeros_like(self.x)
        self.W      = np.zeros_like(self.x)
        self.LV     = np.zeros_like(self.x)
        self.plot_interval = plot_interval


    def integrate(self):
        # pulling instance variables to local scope
        V, W, LV = self.V, self.W, self.LV
        xmax, t, dx, dt, eta = self.x[-1], self.t, self.dx, self.dt, self.eta

        # initial excitation of action potential
        V[0] = 1.
        # initial boundary condition
        bc = 'neumann'

        # variables for measuring CV, BCL, APD
        threshold   = .1
        previous_V  = 0.
        uptime      = -10.

        for i in range(self.steps-1):
            # compute laplacian and alpha step
            Lap1d(V, LV, dx, bc)
            dV, dW = self.G(V, W)

            # euler step
            V += dt * dV + dt * eta * LV
            W += dt * dW

            # current measurement
            now         = t[i]
            current_V   = V[10]

            if current_V > threshold and previous_V < threshold:
                # detected upstroke
                CV      = xmax / (now - uptime)
                BCL     = now - uptime
                uptime  = now
                print('UPSTROKE:\n\tCV = %e\n\tBCL = %e\n' % (CV, BCL))

            if current_V < threshold and previous_V > threshold:
                # detected downstroke
                APD = now - uptime
                print('DOWNSTROKE:\n\tAPD = %e\n' % APD)

            previous_V = current_V

            # switching boundary condition to periodic, after enough time for
            # the pulse to develop
            if i == 4_000:
                bc = 'periodic'
                print('switche boundary condition to periodic...\n')

            if i % self.plot_interval == 0:
                yield V



params = dict(a=.15, k=8., e0=2e-3, m1=.2, m2=.3)
p = Pulse1d(100, .2, 1000., .3, 50, **params)
pulse_gen = p.integrate()

fig, ax = plt.subplots()
ax.set(title='action potential propagation', xlabel='x', ylabel='V',
       xlim=(0, p.x[-1]), ylim=(-.1, 1.1))
V_line, = ax.plot([], [])


def init():
    V_line.set_data([], [])
    return V_line,


def anim(V):
    V_line.set_data(p.x, V)
    return V_line,


anim = animation.FuncAnimation(fig, anim, frames=pulse_gen, interval=20,
                               init_func=init, blit=True, repeat=False)


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
