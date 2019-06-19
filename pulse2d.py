# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from alpha import AlphaBase


def Lap2d(a, o, dx, xbound='neumann', ybound='neumann', both=None):
    # compute laplacian
    o[...] = (  np.roll(a, 1, 0) + np.roll(a, -1, 0)
              + np.roll(a, 1, 1) + np.roll(a, -1, 1)
              - 4 * a) / dx**2

    # apply boundary conditions seperatly for x and y
    v = 0.
    if both:
        xbound = ybound = both
    if xbound == 'neumann':
        o[0,:], o[-1,:] = o[1,:], o[-2,:]
    elif xbound == 'dirichlet':
        o[0,:], o[-1,:] = v, v
    elif xbound == 'periodic':
        pass
    if ybound == 'neumann':
        o[:,0], o[:,-1] = o[:,1], o[:,-2]
    elif ybound == 'dirichlet':
        o[:,0], o[:,-1] = v, v
    elif ybound == 'periodic':
        pass

    return o



class Pulse2d(AlphaBase):
    def __init__(self, xmax, ymax, dx, tmax, eta, plot_interval=50,
                 **alpha_params):
        super().__init__(**alpha_params)

        self.dx     = dx
        self.eta    = eta
        self.dt     = .2 * dx**2 / (2. * eta)

        self.X, self.Y  = np.mgrid[0:xmax:dx,0:ymax:dx]
        self.steps      = int(tmax // self.dt + 1)

        self.V  = np.zeros_like(self.X)
        self.W  = np.zeros_like(self.X)
        self.LV = np.zeros_like(self.X)

        self.plot_interval = plot_interval


    def integrate(self):
        V, W, LV = self.V, self.W, self.LV
        dx, dt, eta = self.dx, self.dt, self.eta

        for i in range(self.steps-1):
            Lap2d(V, LV, dx, both='neumann')
            dV, dW = self.G(V, W)

            V += dt * dV + dt * eta * LV
            W += dt * dW

            if i % self.plot_interval == 0:
                yield V, W


    def integrate_channel(self):
        V, W, LV = self.V, self.W, self.LV
        dx, dt, eta = self.dx, self.dt, self.eta

        xbound = 'neumann'
        ybound = 'dirichlet'

        for i in range(self.steps-1):
            Lap2d(V, LV, dx, xbound=xbound, ybound=ybound)
            dV, dW = self.G(V, W)

            V += dt * dV + dt * eta * LV
            W += dt * dW

            if i == 500:
                xbound = 'periodic'

            if i % self.plot_interval == 0:
                yield V, W


    def integrate_spiral(self, delay):
        V, W, LV = self.V, self.W, self.LV
        dx, dt, eta = self.dx, self.dt, self.eta

        for i in range(self.steps-1):
            Lap2d(V, LV, dx, both='neumann')
            dV, dW = self.G(V, W)

            V += dt * dV + dt * eta * LV
            W += dt * dW

            if i == delay:
                xm, ym = self.X.shape[0] // 2, self.X.shape[1] // 2
                V[xm-5:xm+5,ym-5:ym+5] = 1.

            if i % self.plot_interval == 0:
                yield V, W



params = dict(a=.15, k=8., e0=2e-3, m1=.2, m2=.3)
#  p = Pulse2d(128, 64, .5, 1000, .3, 50, **params)
p = Pulse2d(128, 256, .5, 1000, .3, 50, **params)
#  p.V[0,:] = 1.
p.V[128:130,0:128] = 1.
p.W[130:132,0:128] = 1.

#  fig, (aV, aW) = plt.subplots(1, 2)
fig, aV = plt.subplots()
aV.grid(False)
#  aW.grid(False)
V_img = aV.imshow(p.V, animated=True)
#  W_img = aW.imshow(p.V, animated=True)

def step(arg):
    #  v, w = arg
    v = arg[0]
    V_img.set_data(v)
    #  W_img.set_data(w)
    return V_img, # W_img,

# double spirals: 1417 <= delay <= 1546
anim = animation.FuncAnimation(fig, step, frames=p.integrate(),
                               interval=20, blit=True, repeat=False)



#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
