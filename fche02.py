# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable


def Lap2d(a, o, dx):
    # compute laplacian
    o[...] = (  np.roll(a, 1, 0) + np.roll(a, -1, 0)
              + np.roll(a, 1, 1) + np.roll(a, -1, 1)
              - 4 * a) / dx**2

    # von Neumann boundary conditions
    o[0,:], o[-1,:] = o[1,:], o[-2,:]
    #  o[:,0], o[:,-1] = o[:,1], o[:,-2]
    # Dirichlet boundary conditions
    #  o[0,:], o[-1,:] = 0., 0.
    #  o[:,0], o[:,-1] = 0., 0.

    return o



class FCHE_Base:

    def __init__(self, Vc, Vc_si, Vv, k,
                 t0, td, tr, tsi, tv1, tv2, tvp, twm, twp,
                 c=1., _2V_model=False):
        # set simulation params
        self.__dict__.update(dict(
            Vc=Vc, Vc_si=Vc_si, Vv=Vv, k=k, c=c, t0=t0, td=td, tr=tr, tsi=tsi,
            tv1=tv1, tv2=tv2, tvp=tvp, twm=twm, twp=twp
        ))

        if _2V_model:
            self.GV = self._GV_2V
            self.Gw = lambda self, *unused: 0.
        else:
            self.GV = self._GV_3V
            self.Gw = self._Gw


    def _I_fi(self, V, v, p):
        return -v * p * (V - self.Vc) * (1 - V) / self.td

    def _I_so(self, V, p):
        return V * (1 - p) / self.t0 + p / self.tr

    def _I_si(self, V, w):
        return -w * (1 + np.tanh(self.k * (V - self.Vc_si))) / (2 * self.tsi)

    def Gv(self, v, p, q):
        return (1 - p) * (1 - v) / ((1 - q) * self.tv1 + q * self.tv2) \
                - p * v / self.tvp

    def _Gw(self, w, p):
        return (1 - p) * (1 - w) / self.twm - p * w / self.twp

    def _GV_3V(self, V, v, w, p):
        return -(self._I_fi(V, v, p) + self._I_so(V, p) + self._I_si(V, w)) / self.c

    def _GV_2V(self, V, v, _unused, p):
        return -(self._I_fi(V, v, p) + self._I_so(V, p)) / self.c



class FCHE_Single_Cell(FCHE_Base):

    def __init__(self, V0, v0, w0, tmax, dt, **params):
        super().__init__(**params)
        self.t      = np.arange(0, tmax, dt)
        self.dt     = dt
        self.steps  = self.t.size

        self.V = np.zeros_like(self.t)
        self.v = np.zeros_like(self.t)
        self.w = np.zeros_like(self.t)

        self.V[0] = V0
        self.v[0] = v0
        self.w[0] = w0


    def integrate(self):
        V, v, w = self.V, self.v, self.w

        for i in range(self.steps-1):
            p = 0. if V[i] < self.Vc else 1.
            q = 0. if V[i] < self.Vv else 1.
            V[i+1] = V[i] + self.dt * self.GV(V[i], v[i], w[i], p)
            v[i+1] = v[i] + self.dt * self.Gv(v[i], p, q)
            w[i+1] = w[i] + self.dt * self.Gw(w[i], p)



class FCHE_2D(FCHE_Base):

    def __init__(self, xmax, ymax, dx, tmax, eta,
                 plot_interval=50, **params):
        super().__init__(**params)

        self.dx     = dx                        # lattice constant
        self.eta    = eta                       # diffusity
        self.dt     = .2 * dx**2 / (2. * eta)   # cfl condition
        self.steps  = int(tmax // self.dt + 1)

        self.X, self.Y = np.mgrid[0:xmax:dx,0:ymax:dx]
        self.V  = np.zeros_like(self.X)
        self.v  = np.zeros_like(self.X)
        self.w  = np.zeros_like(self.X)
        self.LV = np.zeros_like(self.X)

        self.plot_interval = plot_interval


    def integrate(self):
        # pull params and arrays to local scope
        Vc, Vv, dx, dt, eta = self.Vc, self.Vv, self.dx, self.dt, self.eta
        V, v, w, LV = self.V, self.v, self.w, self.LV

        p = np.ones_like(V)
        q = np.ones_like(V)

        for i in range(self.steps-1):
            p[:] = q[:] = 1.
            p[V < Vc] = 0.
            q[V < Vv] = 0.

            Lap2d(V, LV, dx)
                  #  both='neumann')
            V += dt * eta * LV + dt * self.GV(V, v, w, p)
            v += dt * self.Gv(v, p, q)
            w += dt * self.Gw(w, p)

            if i % self.plot_interval == 0:
                yield 100. * V - 80.



PARAM_SETS = {
    1 : dict(
        tvp     = 3.33,
        tv1     = 19.6,
        tv2     = 1000.,
        twp     = 667.,
        twm     = 11.,
        td      = .25,
        t0      = 8.3,
        tr      = 50.,
        tsi     = 45.,
        k       = 10.,
        Vc_si   = .85,
        Vc      = .13,
        Vv      = .055
    ),

    2 : dict(
        tvp     = 10.,
        tv1     = 10.,
        tv2     = 10.,
        twp     = 0.,
        twm     = 0.,
        td      = .25,
        t0      = 10.,
        tr      = 190.,
        tsi     = 0.,
        k       = 0.,
        Vc_si   = 0.,
        Vc      = .13,
        Vv      = 0.,
        _2V_model = True
    ),

    3 : dict(
        tvp     = 3.33,
        tv1     = 19.6,
        tv2     = 1250.,
        twp     = 870.,
        twm     = 41.,
        td      = .25,
        t0      = 12.5,
        tr      = 33.33,
        tsi     = 29.,
        k       = 10.,
        Vc_si   = .85,
        Vc      = .13,
        Vv      = .04
    ),

    4 : dict(
        tvp     = 3.33,
        tv1     = 15.6,
        tv2     = 5.,
        twp     = 350.,
        twm     = 80.,
        td      = .407,
        t0      = 9.,
        tr      = 34.,
        tsi     = 26.5,
        k       = 15.,
        Vc_si   = .45,
        Vc      = .15,
        Vv      = .04
    )
}



def plot_single_cell():
    fig, ((aV, av, aw), (aIfi, aIsi, aIso)) = plt.subplots(2, 3)

    aV.set(title='action potential', xlabel='t/ms', ylabel='V/mV')
    av.set(title='fast gate variable', xlabel='t/ms')
    aw.set(title='slow gate variable', xlabel='t/ms')
    aIfi.set(title='fast inward current', xlabel='t/ms')
    aIsi.set(title='slow inward current', xlabel='t/ms')
    aIso.set(title='slow outward current', xlabel='t/ms')

    for i in PARAM_SETS.keys():
        sim = FCHE_Single_Cell(.3, 1., 1., 400, .01, **PARAM_SETS[i])
        sim.integrate()

        p = np.ones_like(sim.V)
        p[sim.V < sim.Vc] = 0.
        V_rescaled = 100 * sim.V - 80

        aV.plot(sim.t, V_rescaled, label='param set %d' % i)
        av.plot(sim.t, sim.v)
        aw.plot(sim.t, sim.w)
        aIfi.plot(sim.t, sim._I_fi(V_rescaled, sim.v, p))
        aIsi.plot(sim.t, sim._I_si(V_rescaled, sim.w))
        aIso.plot(sim.t, sim._I_so(V_rescaled, p))

    aV.legend()


def plot_2d_tissue(i):
    sim = FCHE_2D(64, 256, 1., 5000, .3, **PARAM_SETS[i])
    sim.V[:,80:120]  = .3
    sim.v[:,100:130] = 1.
    sim.w[:,80:120]  = 1.

    fig, aV = plt.subplots()
    aV.axis("off")
    aV.grid(False)

    #  div = make_axes_locatable(aV)
    #  cax = div.append_axes("bottom", "5%", "5%")

    V_img = aV.imshow(sim.V, animated=True, cmap=plt.get_cmap("plasma"))
    fig.colorbar(V_img) # , cax=cax)

    def step(arg):
        V_img.set_data(arg)
        V_img.set_clim(arg.min(), arg.max())
        return V_img,

    anim = animation.FuncAnimation(fig, step, frames=sim.integrate(),
                                   interval=20, blit=True, repeat=False)

    return sim, fig, anim


if __name__ == "__main__":
    s, f, a = plot_2d_tissue(4)


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
