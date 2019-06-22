# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from alpha import AlphaBase


def Lap2d(a, o, dx, xbound='neumann', ybound='neumann', both=None):
    """Compute Laplacian of 2d array with lattice constant dx and apply
    boundary conditions (seperatly for x and y coordinates).
    Store output in o.

    Available boundary conditions:
        * periodic      L[0] = L[end-1]
        * von Neumann   L[0] = L[1]
        * Dirichlet     L[0] = 0
    """
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
    """class for integrating 2d tissue with aliev-panfilov model"""

    def __init__(self, xmax, ymax, dx, tmax, eta, plot_interval=50,
                 **alpha_params):
        """params should be self-explanatory

        notes:
            * dt is computed to fullfill the cfl condition
                dt < dx**2 / (2 * eta)
            * plot_interval is the number of integration steps after which
              to yield an intermediate result for plotting
            * allocate three arrays:
                V, W: model variables
                LV: to store the Laplacian of V
        """
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


    def integrate(self, xbound='neumann', ybound='periodic'):
        """integrate the system using a simple euler step with chosen
        boundary conditions (applied when computing the Laplacian of V)

        perform `self.steps` integration steps, yield the arrays V, W for
        plotting every `self.plot_interval` steps
        """
        V, W, LV = self.V, self.W, self.LV
        dx, dt, eta = self.dx, self.dt, self.eta

        for i in range(self.steps-1):
            Lap2d(V, LV, dx, xbound, ybound)
            dV, dW = self.G(V, W)

            V += dt * dV + dt * eta * LV
            W += dt * dW

            if i % self.plot_interval == 0:
                yield 100. * V - 80., W



params = dict(a=.15, k=8., e0=2e-3, m1=.2, m2=.3)


def channel():
    p = Pulse2d(32, 128, .5, 1000, .3, **params)
    p.V[:,30:32] = 1.
    p.W[:,28:30] = 1.

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, width_ratios=[16, 1])
    aV = plt.subplot(gs[0,0])
    aW = plt.subplot(gs[1,0])
    ac = plt.subplot(gs[:,1])

    aV.grid(False)
    aV.axis('off')
    aV.set_title('action potential')
    aW.grid(False)
    aW.axis('off')
    aW.set_title('relaxation variable')
    ac.set_xlabel('V/mV')

    V_img = aV.imshow(p.V, animated=True, cmap=plt.get_cmap('plasma'))
    W_img = aW.imshow(p.V, animated=True, cmap=plt.get_cmap('plasma'))
    fig.colorbar(V_img, cax=ac)

    def step(a):
        V_img.set_data(a[0])
        W_img.set_data(a[1])
        V_img.set_clim(a[0].min(), a[0].max())
        return V_img, W_img

    anim = animation.FuncAnimation(fig, step,
                                   frames=p.integrate(xbound='dirichlet'),
                                   interval=20, blit=True, repeat=False)

    return p, fig, anim


def spiral_excitation(delay):
    p = Pulse2d(64, 128, .5, 1000, .3, **params)
    p.V[:,30:32] = 1.
    p.W[:,28:30] = 1.

    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, width_ratios=[16, 1])
    aV = plt.subplot(gs[0,0])
    aW = plt.subplot(gs[1,0])
    ac = plt.subplot(gs[:,1])

    aV.grid(False)
    aV.axis('off')
    aV.set_title('action potential')
    aW.grid(False)
    aW.axis('off')
    aW.set_title('relaxation variable')
    ac.set_xlabel('V/mV')

    V_img = aV.imshow(p.V, animated=True, cmap=plt.get_cmap('plasma'))
    W_img = aW.imshow(p.V, animated=True, cmap=plt.get_cmap('plasma'))
    fig.colorbar(V_img, cax=ac)

    gen = p.integrate(ybound='neumann')

    def step(i):
        if i == delay:
            xm, ym = p.X.shape[0] // 2, p.X.shape[1] // 2
            p.V[xm-5:xm+5,ym-5:ym+5] = 1.

        a = next(gen)
        V_img.set_data(a[0])
        V_img.set_clim(a[0].min(), a[0].max())
        W_img.set_data(a[1])

        return V_img, W_img

    anim = animation.FuncAnimation(fig, step,
                                   interval=20, blit=True, repeat=False)

    return p, fig, anim


def spiral_wave():
    p = Pulse2d(128, 256, .5, 1000, .3, **params)
    p.V[128:130,0:128] = 1.
    p.W[130:132,0:128] = 1.

    fig, aV = plt.subplots()
    aV.grid(False)
    aV.axis('off')
    aV.set_title('action potential')

    div = make_axes_locatable(aV)
    cax = div.append_axes('right', '5%', '5%')
    cax.set_xlabel('V/mV')

    V_img = aV.imshow(p.V, animated=True, cmap=plt.get_cmap('plasma'))
    fig.colorbar(V_img, cax=cax)

    def step(a):
        V_img.set_data(a[0])
        V_img.set_clim(a[0].min(), a[0].max())
        return V_img,

    #  FFWriter = animation.FFMpegWriter(fps=30)
    anim = animation.FuncAnimation(fig, step,
                                   frames=p.integrate(ybound='neumann'),
                                   interval=20, blit=True, repeat=False)
    #  anim.save('spiral.mp4', writer=FFWriter, dpi=300)

    return p, fig, anim


if __name__ == '__main__':
    #  p, f, a = channel()
    #  p, f, a = spiral_excitation(23)
    p, f, a = spiral_wave()

    plt.show()


#  vim: set ff=unix tw=79 sw=4 ts=8 et ic ai :
