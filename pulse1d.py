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
        self.x      = np.arange(0, xmax, dx)
        self.steps  = int(np.ceil(tmax / self.dt))
        self.V      = np.zeros_like(self.x)
        self.W      = np.zeros_like(self.x)
        self.LV     = np.zeros_like(self.x)
        self.plot_interval = plot_interval


    def integrate(self):
        V, W, LV = self.V, self.W, self.LV
        dx, dt, eta = self.dx, self.dt, self.eta
        V[0] = 1.
        bc = 'neumann'

        for i in range(self.steps-1):
            Lap1d(V, LV, dx, bc)
            dV, dW = self.G(V, W)
            V += dt * dV + dt * eta * LV
            W += dt * dW

            if i == 4_000:
                bc = 'periodic'

            if i % self.plot_interval == 0:
                yield V



fig, ax = plt.subplots()
ax.set(title='action potential propagation', xlabel='x', ylabel='V',
       xlim=(0, 200), ylim=(-.1, 1.1))
V_line, = ax.plot([], [])

params = dict(a=.15, k=8., e0=2e-3, m1=.2, m2=.3)
p = Pulse1d(200, .2, 5000., .3, 50, **params)
pulse_gen = p.integrate()


def init():
    V_line.set_data([], [])
    return V_line,


def anim(V):
    V_line.set_data(p.x, V)
    return V_line,


anim = animation.FuncAnimation(fig, anim, frames=pulse_gen, interval=20,
                               init_func=init, blit=True, repeat=False)


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
