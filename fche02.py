# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def Lap2d(a, o, dx):
    # compute laplacian
    o[...] = (  np.roll(a, 1, 0) + np.roll(a, -1, 0)
              + np.roll(a, 1, 1) + np.roll(a, -1, 1)
              - 4 * a) / dx**2

    # von Neumann boundary conditions
    o[0,:], o[-1,:] = o[1,:], o[-2,:]
    o[:,0], o[:,-1] = o[:,1], o[:,-2]

    return o



class FCHE:

    def __init__(self, xmax, ymax, dx, tmax, eta, c,
                 Vc, Vc_si, Vv, k, t0, td, tr, tsi, tv1, tv2, tvp, twm, twp,
                 plot_interval=50):
        # set simulation params
        self.__dict__.update(dict(
            Vc=Vc, Vc_si=Vc_si, Vv=Vv, k=k, t0=t0, td=td, tr=tr, tsi=tsi,
            tv1=tv1, tv2=tv2, tvp=tvp, twm=twm, twp=twp
        ))

        self.c      = c     # specific capacity
        self.dx     = dx    # lattice constant
        self.eta    = eta   # diffusity
        self.dt     = .2 * dx**2 / (2. * eta)   # cfl condition
        self.steps  = tmax // self.dt + 1.

        self.X, self.Y = np.mgrid[0:xmax:dx,0:ymax:dx]
        self.V = np.zeros_like(self.X)

        self.plot_interval = plot_interval


    def G(self, V, v, w, p, q, out):
        # pull params to local scope
        I_fi, I_so, I_si    = self.I_fi, self.I_so, self.I_si
        Vc, Vc_si, k, c     = self.Vc, self.Vc_si, self.k, self.c
        t0, td, tr, tsi     = self.t0, self.td, self.tr, self.tsi

        # compute currents
        I_fi[:] =  -v * p * (V-Vc) * (1-V) / td
        I_so[:] = V * (1-p) / t0 + p / tr
        I_si[:] = -w * (1 + np.tanh(k * (V-Vc_si))) / (2 * tsi)
        out[:]  = -(I_fi + I_so + I_si) / c

        return out


    def integrate(self):
        # pull params to local scope
        Vc, Vv, dx, dt, eta     = self.Vc, self.Vv, self.dx, self.dt, self.eta
        tv1, tv2, tvp, twm, twp =  self.tv1, self.tv2, self.tvp, self.twm, self.twp

        # init arrays
        V   = self.V
        dV  = np.zeros_like(V)
        LV  = np.zeros_like(V)
        v   = np.zeros_like(V)
        w   = np.zeros_like(V)
        p   = np.ones_like(V)
        q   = np.ones_like(V)

        for i in range(self.steps-1):
            p[:] = q[:] = 1.
            p[V < Vc] = 0.
            q[V < Vv] = 0.

            LV = Lap2d(V, dx, out=LV, both='neumann')
            dV = self.G(V, out=dV)
            V += dt * eta * LV + dt * dV
            v += dt * (1-p) * (1-v) / ((1-q) * tv1 + q * tv2) - dt * p * v / tvp
            w += dt * (1-p) * (1-w) / twm - dt * p * w / twp



#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai : 
