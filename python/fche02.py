# -*- coding: utf-8 -*-
"""Phenomenological ionic model used by Fenton, Cherry, Hastings and Evans
in 'Multiple mechanisms of spiral wave breakup in a model of cardiac
electrical activity' (2002); here referred to as 'FCHE02'.

Similar to the Hodgkin-Huxley model the membrane potential is said to obey
the cable equation:
    c * dV_m/dt = -I_m
    I_m = I_fi + I_so + I_si

The components of the membran current I_m:
    * fast inward current
        I_fi = -v * Θ(V-Vc) * (V-Vc) * (1-V) / td

     ** depolarizes membrane upon an excitation above Vc
     ** depends on fast activation gate Θ(V-Vc) and fast inactivation gate v

    * slow outward current
        I_so = V * (1-Θ(V-Vc)) / t0 + Θ(V-Vc) / tr

     ** repolarizes membrane back to resting potential
     ** depends on fast activation gate Θ(V-Vc)

    * slow inward current
        I_si = -w * d / (2 * tsi),  d -> 1 + tanh(k * (V-Vc_si))

     ** inactivation current to balance I_so and to produce the observed
     plateau in the action potential
     ** depends on the slow inactivation gate w and on the very fast
     activation gate d, which is modeled by a steady-state function


The two gate variables governing the currents:
    * fast inactivation gate
        dv/dt = (1-Θ(V-Vc)) * (1-v) / tvm - Θ(V-Vc) * v / tvp,
            tvm = (1-Θ(V-Vv)) * tvm1 + Θ(V-Vv) * tvm2

    * slow inactivation gate
        dw/dt = (1-Θ(V-Vc)) * (1-w) / twm - Θ(V-Vc) * w / twp


Parameters:
    * tvp, tvm1, tvm2: opening ([p]lus) and closing ([m]inus) times of the
    fast variable v
    * twp. twm: opening and closing times of the slow variable w
    * td, tr: de- and repolarization times
    * t0, tsi: time constants for slow currents
    * Vc, Vv, Vc_si: voltage thresholds
    * k: activation width parameter
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pulse2d import Lap2d



class FCHE_Base:
    """Base class for FCHE simulations to manage the necessary parameters"""

    def __init__(self, Vc, Vc_si, Vv, k,
                 t0, td, tr, tsi, tv1, tv2, tvp, twm, twp,
                 c=1., _2V_model=False):
        """Description of parameters: see module docstring

        _2V_model: if True, ignore slow inward current I_si and reduce
            model from three to two variables
        """
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
    """Class to examine behaviour of the action potential of a single
    cell"""

    def __init__(self, V0, v0, w0, tmax, dt, **params):
        """Set initial values for V, v and w; define tmax and integration
        time step"""
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
        """integrate system with a simple Euler step"""
        V, v, w = self.V, self.v, self.w

        for i in range(self.steps-1):
            p = 0. if V[i] < self.Vc else 1.
            q = 0. if V[i] < self.Vv else 1.
            V[i+1] = V[i] + self.dt * self.GV(V[i], v[i], w[i], p)
            v[i+1] = v[i] + self.dt * self.Gv(v[i], p, q)
            w[i+1] = w[i] + self.dt * self.Gw(w[i], p)



class FCHE_2D(FCHE_Base):
    """Class to apply the FCHE model for 2D tissue simulations"""

    def __init__(self, xmax, ymax, dx, tmax, eta,
                 plot_interval=50, **params):
        """Define shape of tissue with xmax, ymax; set tmax; compute dt
        from dx and eta in order to fulfill the cfl condition, i.e.:
            dt < dx**2 / (2*eta)

        Yield intermediate result every plot_interval steps for plotting
        """
        super().__init__(**params)

        self.dx     = dx                        # lattice constant
        self.eta    = eta                       # diffusivity
        self.dt     = .2 * dx**2 / (2. * eta)   # cfl condition
        self.steps  = int(tmax // self.dt + 1)

        self.X, self.Y = np.mgrid[0:xmax:dx,0:ymax:dx]
        self.V  = np.zeros_like(self.X)
        self.v  = np.zeros_like(self.X)
        self.w  = np.zeros_like(self.X)
        self.LV = np.zeros_like(self.X)

        self.plot_interval = plot_interval


    def integrate(self, xbound='neumann', ybound='periodic'):
        """Integrate system using a simple Euler step

        In order to consider spatial dynamics, the cable equation governing
        V_m is expanded with a diffusion term:
            c*dV_m/dt = eta * Lap(V_m) + I_m
        """
        # pull params and arrays to local scope
        Vc, Vv, dx, dt, eta = self.Vc, self.Vv, self.dx, self.dt, self.eta
        V, v, w, LV = self.V, self.v, self.w, self.LV

        p = np.ones_like(V)
        q = np.ones_like(V)

        for i in range(self.steps-1):
            p[:] = q[:] = 1.
            p[V < Vc] = 0.
            q[V < Vv] = 0.

            Lap2d(V, LV, dx, xbound=xbound, ybound=ybound)
            V += dt * eta * LV + dt * self.GV(V, v, w, p)
            v += dt * self.Gv(v, p, q)
            w += dt * self.Gw(w, p)

            if i % self.plot_interval == 0:
                yield 100. * V - 80.



# predefined parameters, taken from the FCHE paper

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
    """for a single cell, plot temporal profile of action potential, fast
    and slow variable and the currents for all four predefined parameter
    sets"""
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


def channel(i):
    """animate channel dynamics (dirichlet conditions on the sides,
    periodic conditions at start and end) for parameter set i """
    sim = FCHE_2D(64, 256, 1., 5000, .3, 30, **PARAM_SETS[i])
    sim.V[:,10:30] = .3
    sim.v[:,20:40] = 1.
    sim.w[:,10:30] = 1.

    fig, aV = plt.subplots()
    aV.axis("off")
    aV.grid(False)
    aV.set_title('action potential - param set %d' % i)

    div = make_axes_locatable(aV)
    cax = div.append_axes('right', '5%', '5%')
    cax.set_xlabel('V/mV')

    V_img = aV.imshow(sim.V, animated=True, cmap=plt.get_cmap("plasma"))
    V_img.set_clim(-80, 25)
    fig.colorbar(V_img, cax=cax)

    def step(arg):
        V_img.set_data(arg)
        V_img.set_clim(arg.min(), arg.max())
        return V_img,

    #  FFWriter = animation.FFMpegWriter(fps=10)
    anim = animation.FuncAnimation(
        fig, step, frames=sim.integrate(xbound='dirichlet'), interval=20,
        blit=True, repeat=False)
    #  anim.save('channel-vid-params-%d.mp4', writer=FFWriter, dpi=300)

    return sim, fig, anim


def spiral_excitation(i, delay):
    """excite an AP wave in a domain with neumann conditions on all sides
    and add a point-like excitations after delay steps

    if timed correctly, this induces spiral waves in the wake of the inital
    wave front"""
    sim = FCHE_2D(128, 512, 1., 500, .3, 30, **PARAM_SETS[i])
    sim.V[:,10:30] = .3
    sim.v[:,20:40] = 1.
    sim.w[:,10:30] = 1.

    fig, aV = plt.subplots()
    aV.axis("off")
    aV.grid(False)
    aV.set_title('action potential - param set %d' % i)

    div = make_axes_locatable(aV)
    cax = div.append_axes('right', '5%', '5%')
    cax.set_xlabel('V/mV')

    V_img = aV.imshow(sim.V, animated=True, cmap=plt.get_cmap("plasma"))
    V_img.set_clim(-80, 25)
    fig.colorbar(V_img, cax=cax)
    fig.tight_layout()

    gen = sim.integrate(ybound='neumann')

    def step(i):
        if i == delay:
            xm, ym = sim.X.shape[0] // 2, sim.X.shape[1] // 2
            sim.V[xm-5:xm+5,ym-5:ym+5] = 1.

        a = next(gen)
        V_img.set_data(a)
        #  V_img.set_clim(a.min(), a.max())

        return V_img,

    #  FFWriter = animation.FFMpegWriter(fps=10)
    anim = animation.FuncAnimation(
        fig, step, frames=sim.steps//sim.plot_interval, interval=20,
        blit=True, repeat=False)
    #  anim.save('spiral-excitation-vid-params%d.mp4' % i,
    #            writer=FFWriter, dpi=300)

    return sim, fig, anim


def spiral_wave(i):
    # nice results with param set 1
    sim = FCHE_2D(256, 512, 1., 10_000, .3, **PARAM_SETS[i])
    sim.V[80:120,:128]  = .3
    sim.v[100:130,:128] = 1.
    sim.w[80:120,:128]  = 1.

    fig, aV = plt.subplots()
    aV.axis("off")
    aV.grid(False)
    aV.set_title('action potential - param set %d' % i)

    div = make_axes_locatable(aV)
    cax = div.append_axes('right', '5%', '5%')
    cax.set_xlabel('V/mV')

    V_img = aV.imshow(sim.V, animated=True, cmap=plt.get_cmap('plasma'))
    V_img.set_clim(-80, 25)
    fig.colorbar(V_img, cax=cax)

    def step(arg):
        V_img.set_data(arg)
        return V_img,

    #  FFWriter = animation.FFMpegWriter(fps=10)
    anim = animation.FuncAnimation(
        fig, step,
        frames=sim.integrate(ybound='neumann'),
        #  frames=sim.integrate(xbound='periodic'),
        interval=20, blit=True, # repeat=False,
        cache_frame_data=False)
    #  anim.save('breakup-%d-04.mp4' % i, writer=FFWriter, dpi=300)

    return sim, fig, anim


if __name__ == "__main__":
    #  s, f, a = channel(4)
    #  s, f, a = spiral_excitation(4, 81)      # plot_interval = 30
    s, f, a = spiral_wave(1)

    #  plt.show()


#  vim: set ff=unix tw=79 sw=4 ts=4 et ic ai :
