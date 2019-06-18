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



class FCHE_Base:
    def __init__(self, Vc, Vc_si, Vv, k, c,
                 t0, td, tr, tsi, tv1, tv2, tvp, twm, twp):
        # set simulation params
        self.__dict__.update(dict(
            Vc=Vc, Vc_si=Vc_si, Vv=Vv, k=k, c=c, t0=t0, td=td, tr=tr, tsi=tsi,
            tv1=tv1, tv2=tv2, tvp=tvp, twm=twm, twp=twp
        ))

    def _I_fi(self, V, v, p):
        return -v * p * (V - self.Vc) * (1 - V) / self.td

    def _I_so(self, V, p):
        return V * (1 - p) / self.t0 + p / self.tr

    def _I_si(self, V, w):
        return -w * (1 + np.tanh(self.k * (V - self.Vc_si))) / (2 * self.tsi)

    def Gv(self, v, p, q):
        return (1 - p) * (1 - v) / ((1 - q) * self.tv1 + q * self.tv2) \
                - p * v / self.tvp

    def Gw(self, w, p):
        return (1 - p) * (1 - w) / self.twm - p * w / self.twp

    def GV(self, V, v, w, p):
        return -(self._I_fi(V, v) + self._I_so(V, p) + self._I_si(V, w)) / self.c



class FCHE_Single_Cell(FCHE_Base):
    def __init__(self, V0, v0, w0, tmax, dt, **params):
        super().__init__(**params)
        self.t      = np.arange(0, tmax, dt)
        self.dt     = dt
        self.steps  = self.t.size

        self.V = np.zeros_like(self.t)
        self.v = np.zeros_like(self.t)
        self.w = np.zeros_like(self.t)


    def integrate(self):
        V, v, w = self.V, self.v, self.w
        p = np.ones_like(V)
        q = np.ones_like(V)

        for i in range(self.steps-1):
            p[:] = q[:] = 1.
            p[V < self.Vc] = 0.
            q[V < self.Vv] = 0.

            V += self.dt * self.GV(V, v, w, p)
            v += self.dt * self.Gv(v, p, q)
            w += self.dt * self.Gw(w, p)



class FCHE_2D(FCHE_Base):

    def __init__(self, xmax, ymax, dx, tmax, eta,
                 plot_interval=50, **params):
        super().__init__(**params)

        self.dx     = dx                        # lattice constant
        self.eta    = eta                       # diffusity
        self.dt     = .2 * dx**2 / (2. * eta)   # cfl condition
        self.steps  = tmax // self.dt + 1.

        self.X, self.Y = np.mgrid[0:xmax:dx,0:ymax:dx]
        self.V  = np.zeros_like(self.X)
        self.v  = np.zeros_like(self.X)
        self.w  = np.zeros_like(self.X)
        self.LV = np.zeros_like(self.X)

        self.plot_interval = plot_interval


    def integrate(self):
        # pull params and arrays to local scope
        Vc, Vv, dx, dt, eta = self.Vc, self.Vv, self.dx, self.dt, self.eta
        V, v, w, LV = self.V, self.v, self.w. self.LV

        p = np.ones_like(V)
        q = np.ones_like(V)

        for i in range(self.steps-1):
            p[:] = q[:] = 1.
            p[V < Vc] = 0.
            q[V < Vv] = 0.

            Lap2d(V, dx, out=LV, both='neumann')
            V += dt * eta * LV + dt * self.GV(V, v, w, p)
            v += dt * self.Gv(v, p, q)
            w += dt * self.Gw(w, p)

            if i % self.plot_interval == 0:
                yield V



#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
