# Content
1. Introduction
    1. Cell structure
    2. Membrane potential
    3. Continuum description
2. Three models
    1. Hudgkin & Huxley (1952)
    2. Aliev & Panfilov (1996)
    3. Fenton et al (2002)
3. Dynamics of a single cell
4. Spatial Dynamics
    1. In 1D
    2. In 2D
5. Discussion  
A. Appendix  
    1. Aliev & Panfilov: Parameter variation
    2. The code
    3. Implementation of boundary conditions


# 1. Introduction
The human heart is a fascinating apparatus, which does its work in a
constant and reliable fashion - usually without disruption - for the whole
of a persons life. To put things into perspective, the heart of an average
human being performs
```
    70 (typical rest pulse per minute)
x 1440 (minutes per day)
x  365 (days per year)
x   80 (estimated average human life time)
~ 3e9 beats per life time,
```
while a typical car engine performs
```
300,000 (km driven during the cars life time)
/ 50    (km/h, average speed)   <- cars life time in hours
x 60    (minutes per hour)
x 2,200 (revolutions per minute)
~ 8e8 duty cycles per life time.
```

Investigating the hearts physical working principle poses an interesting
challenge with undoubtedly many relevant applications such as
gaining a deeper understanding of and developing more advanced treatments
for - many times dangerous - arrhythmia.

The heart muscles basic functionality is rhythmically contracting itself
triggered by electrical signals, which are being conducted by the heart
muscle cells (the cardiomyocytes) themselves (this is the remarkable thing
here, because usually this task is being taken care of by neural cells),
which means that this kind of tissue combines the ability to both perform
mechanical work and conduct electrical signals.

Going a little more into detail: pacemaker cells (specialized
cardiomyocytes) in the SA-node rhythmically generate action potential
which travel at about .05-1 m/s to the AV-node and from there after a
delay at about 2-4 m/s through the ventricular bundles and the Purkinje
fibres. Ultimately those action potentials cause the contraction of the
regions of heart tissue controlled by the respective bundles of conduction
cells.


## 1.1. Cell structure
For a better understanding it is helpful to take a closer look at the
microscopic structure of the cardiomyocytes:
 * tubular cells containing chains of _myofibril_ (fibres composed of long
   proteins), which are responsible for contraction of the muscle tissue
 * _sarcoplastic reticulum_: membrane-enclosed regions, mainly storing
   Ca^2+ ions
 * enclosed by a double lipid-layered membrane: the _sarcolemma_
In longitudinal direction _interlacing disks_ join the cells together and
via the _gap junctions_ allow propagation of action potentials. Because of
these features the heart muscle forms a _syncytium_, i.e. the single cells
behave like a single coordinated unit.


## 1.2. Membrane potential
Now the interior and exterior (i.e. the intermediate space between
neighbouring cells) regions of a cells exhibit different concentrations of
various ion species (this imbalance is being maintained by special ion
pumps and gates in the cell membrane), which results in a voltage between
those regions: the membrane potential `V=Φ_i-Φ_e`, which in the rest case
is equal to the rest potential `V_rest`.

> TODO: table with typical values for V\_rest

If at some point the membrane potential is perturbed by a stimulus in such
a way that it exceeds some threshold, the ion channels rapidly open causing
the concentration difference of the ions between interior and exterior
cell regions to invert resulting in a large upswing of the membrane
potential. This process is called _Depolarization_ and the peak of the
membrane potential is called _action potential_.

After reaching this peak the gates close again and the pumps recreate the
prior concentration difference which causes the membrane potential to
return to the rest value (_repolarization_).

Since the membrane posses a finite specific electric capacity `C_m` [F m^-2]
the membrane potential obeys the capacitor equation:
```
C V = Q => dQ/dt = I = C dV/dt
```

> TODO: add circuit sketch

Therewith the dynamics of `V_m` can be modeled by modeling the cross
membrane current density `I`.


## 1.3. Continuum description
At a microscopic level the propagation of action potential is a discrete
process (from cell to cell). However looking at tissue at sufficiently
large scales, it can be viewed as continuous (-> functional syncytium). It
is important to note the anisotropic nature of this process: the tissue
exhibits different conductivities in longitudinal and transversal direction
with respect to the myofibril.

### 1.3.1. Bidomain model
One formulates potentials and current densities for the intra- and
extracellular regions: `Φ_i, Φ_e, J_i, J_e`
It is important to note, that formally all of these functions are defined
on the whole domain.

To set them into relation, consider Poisson's equation and Ohm's Law:
```
E=-∇Φ, J=GE=-G∇Φ
=> J_i=-G_i ∇Φ_i, J_e=-G_e ∇Φ_e
```
where E is the electrical field associated with the potential Φ and G is
the conductivity tensor accounting for the anisotropy.
Imposing conservation of current:
```
∇\dot(J_i + J_e)=0 => -∇\dot{J_i}=∇\dot{J_e}=I_m
```
where `I_m` is the transmembrane current density (units A m^-3).
From the capacitor equation
```
I_m=\beta(CdV\dt+\sum_{s}I_s), V=Φ_i-Φ_e
```
(where `\beta` is a scaling constant with units m^-1) one finds:
> TODO: add equations

Now one has a system of two coupled PDEs with (i) being parabolic and (ii)
being elliptic, which is rather difficult to solve.

## 1.3.2. Monodomain model
In order to make matters more accessible one makes the assumption that the
intra- and extracellular anisotropies are identical, i.e. the respective
conductivities are proportional:
```
G_e=\lambda G_i
=> ... =>
dV/dt=∇\dot{D}∇V-1/C\sum_{s}I_s,    D=G/(C\beta)
```
which reduces the problem to one parabolic PDE (a diffusion equation).

Now all that is left to do is to model the conductivity tensor `D` to
represent the tissue at hand.

One way of doing this is to split this object
into components parallel and perpendicular to the direction of the
myofibril \vec{f}, like so
```
D=D_{\perp} I + (D_{\parallel}-D_{\perp}) \vec{f}\vec{f}^{T}
```
where for typical cardiomyocytes one has the relation
`D_{\parallel}/D_{\perp}~2...10`.

Another way is to neglect the anisotropies all together and write the
conductivity tensor as a scalar value: `D -> η`. This is the approach
taken by the following investigations.


_Note_: Without external stimuli (enforced by von Neumann boundary
conditions) bi- and monodomain models yield almost identical results.
However when considering such external stimuli (e.g. defibrillation) the
unequal anisotropies of intra- and extracellular regions are significant.


# 2. Three models
In this section I am going to describe three approaches [(out of an
abundance of available models)](www.cellml.org) which can be used to
describe the dynamics of the membrane potentials.


## 2.1. Hudgkin & Huxley (1952)
This model was developed by A. L. Hodgkin and A. F. Huxley to fit
measurements taken on a giant squid axon prior to detailed knowledge about
the biophysical mechanisms being available.

The membrane current density is modeled as the sum of Sodium, Potassium and
a leakage current, each obeying Ohm's Law:
```
I_m = \sum_{s}I_s = I_{Na} + I_{K} + I_{l}
I_s = g_s (V-V_s)
```
> TODO: add values for V_s

The specific resistivities are described by gating variables
(dimensionless, values in (0, 1)):
```
g_{K} = \bar{g}_{K} n^4
g_{Na} = \bar{g}_{Na} m^3 n
g_{l} = \bar{g}_{l}
```
> TODO: add values for \bar{g}\_{s}

And the gating variables obey the following ODEs:
```
di/dt = \alpha_{i}(1-i)-\beta_{i}i
```
> TODO: add \alpha\_{i}, \beta\_{i}

Thus one has a system of uncoupled four 1st order ODEs (V, n, m, h) to
solve.


## 2.2. Aliev & Panfilov (1996)
While the Hodkin-Huxley model gives a very good description of
action potential dynamics, it is rather complex and therefor not preferable
for large scale computations.

An alternative model was formulated by Aliev and Panfilov, which refrains
from using a description based on biophysical details and instead uses two
variables (the potential V and a relaxation variable W) to
phenomenologically reproduce the membrane potential dynamics of
cardiomyocytes:
```
dV/dt = -k V (V-a) (V-1) - V W
dW/dt = e(V, W) (-k V (V-a-1) - W)
e = e_0 + \frac{\mu_1 W}{V + \mu_2}
```

Here one has to solve two coupled 1st order ODEs.

One can think of the relaxation variable W as summarizing and hiding all
the complex processes involving ion pumps etc. in order to cause the
membrane potential to repolarize.


## 2.3. Fenton et al (2002)
Yet another model to be introduced here is again based on the capacitor
equation. Unlike the models presented above, this approach allows to easily
investigate chaotic wave breakup mechanisms (its technical meaning and
possible physiological interpretation will be presented in section ? and
section ? respectively).

The membrane current density is modeled as the sum of the following
phenomenological current densities:
* fast inward current

    I_fi = -v * Θ(V-Vc) * (V-Vc) * (1-V) / td

    * depolarizes membrane upon an excitation above Vc
    * depends on fast activation gate Θ(V-Vc) and fast inactivation gate v

* slow outward current

    I_so = V * (1-Θ(V-Vc)) / t0 + Θ(V-Vc) / tr

    * repolarizes membrane back to resting potential
    * depends on fast activation gate Θ(V-Vc)

* slow inward current

    I_si = -w * d / (2 * tsi),  d -> 1 + tanh(k * (V-Vc_si))

    * inactivation current to balance I_so and to produce the observe
      plateau in the action potential
    * depends on the slow inactivation gate w and on the very fast
      activation gate d, which is modeled by a steady-state function

The two gate variables governing the currents:
* fast inactivation gate

    dv/dt = (1-Θ(V-Vc)) * (1-v) / tvm - Θ(V-Vc) * v / tvp,
    tvm = (1-Θ(V-Vv)) * tvm1 + Θ(V-Vv) * tvm2

* slow inactivation gate

    dw/dt = (1-Θ(V-Vc)) * (1-w) / twm - Θ(V-Vc) * w / twp


And the parameters:
* tvp, tvm1, tvm2: opening ([p]lus) and closing ([m]inus) times of the
fast variable v
* twp. twm: opening and closing times of the slow variable w
* td, tr: de- and repolarization times
* t0, tsi: time constants for slow currents
* Vc, Vv, Vc_si: voltage thresholds
* k: activation width parameter

The problem to be solved in this description composes of three uncoupled
1st order ODEs (V, v, w).


# 3. Dynamics of a single cell
In this section the results of the numerical solutions for a single cell
using the three models presented above are discussed. All results were
obtained with a simple Forward-Euler scheme.

## 3.1. Hodgkin-Huxley
The equations where integrated for 10 ms with a time step Δt=0.01 ms (i.e.
1000 integration steps) and an initial potential of `V_0=-7mV`.

> TODO: insert hh52-10ms.png

The gating variables were initially all set to 0, yet at t=0ms the Na^+ h
gate is already \~80% open (but the Na^+ current is being suppressed by the
closed m gate). The h gate slowly closes up to t=5ms and then swings down
more rapidly, while at the same time the m gate opens strongly and allows
the influx of Na^+ ions (interestingly the Na^+ current exhibits a small
seperate peak at t\~5ms before rising up to its full strength). This causes
the membrane potential to rise up and become strongly positive.
Shortly after the quick opening of the m gate, the Ka^+ n gate also opens
allowing for an opposed Ka^+ outflux causing the membrane potential to
eventually repolarize.

The resulting action potential resembles qualitatively the results depicted
in Fig. 12 in [HH52]. An even more fit is achieved when integrating the
system for 40 ms with an initial value `V_0=-30 mV`, as depicted in Fig. 22
in [HH52].

> TODO: insert hh52-40ms.png

Another interesting effect can be observed when adding an additional source
current to equation (?), i.e.
```
CdV/dt = -\sum_{s}I_s + I_{source}
```
which causes the membrane potential to depolarize again after
repolarization with a constant rate.

> TODO: insert hh52-40ms-10nA


## 3.2. Aliev-Panfilov
_Note:_ The action potential V in this model takes values in (0, 1) and
needs to be rescaled according to `V_phys/mV=100*V-80` and `t_phys/ms=12.9*t`

The integration was performed for 60 model time units (t\_{phys,max}=774 ms)
with a time step of Δt=0.01 model time units (i.e. 6000 integration steps)
and an initial excitement of the action potential to V\_0=0.2 model voltage
units (V\_{phys,0}=-60 mV).

Upon an excitement which surpasses some threshold the action potential
quickly raises up to its maximal value. The relaxation variable begins
raising, too, slowly at first, but gradually growing faster until it
steeply reaches its peak, which pulls the action potential back to its rest
value.

> TODO: add alpha-phase+dynamics.png
The phaseplot to the right shows the characteristic connection of the two
variables (unscaled).
The plot to the left shows the temporal evolution of the two variables
(scaled)

The models dynamics are governed by a set of five parameters a, k, ε\_0,
μ\_1 and μ\_2. These are phenomenological in their nature and therefor
difficult to interpret. A short study which varies each parameter while
holding the others constant is found in the appendix.


## 3.3. Fenton et al
_Note:_ In this model the action potential also takes values in (0, 1) and
needs to be rescaled with `V_phys=(100*V-80)/mV`, however the time closely
resembles the physical time.

The integration was performed for 400 ms, using a time step Δt=0.1 ms
(resulting in 4000 integration steps) and an initial excitement of the
membrane potential V\_0=0.3 model voltage units (V\_{phys,0}=-50 mV).

The resulting action potential qualitatively resembles the result from the
Aliev-Panfilov model, with the difference that the peak does not start
decreasing right after the excitement but rather shows a slight bump in the
plateau.
Just like the Hodgkin-Huxley model, the Fenton model is based on an ionic
description, which allows to comprehend the action potential dynamics based
on underlying gate and current processes.

> TODO: insert fenton-single-cell.png

Here are two things to note here:
 * the gate variables are _inactivation gates_, i.e. 1 means fully closed,
   while 0 means fully opened.
 * the currents are modeled with Heaviside-functions in order to be
   activated when the potential surpasses a certain threshold (defined by
   the model parameters V\_c and V\_v); this causes the sudden jumps in the
   current plots.


# 4. Continuous Spatial Dynamics
The spatial tissue simulations discussed in this section are based on a
continuous approximation of the tissue, which in truth is composed of discrete
cells. The mathematical approach was introduced in section (1.3). In practice -
regardless of the concrete model used - one adds a diffusion term to the ODE
describing the temporal evolution of the membrane potential and obtains a
parabolic PDE:
```
dV/dt=G(V) --> dV/dt=ηΔV+G(V)
```
In this case the - normally anisotropic - diffusivity is given
by the scalar coefficient η, which has been set to η=0.3 in all cases presented
below.
Additionally, the time step is computed to fulfill the CFL condition:
```
Δt<Δx^2/(2η)
```

## 4.1. In 1D
The simulation is performed using the Aliev-Panfilov model. The calculations
are performed with two one-dimensional arrays representing the values of the
membrane potential and the relaxation variable along a discretized line,
whose length has been varied throughout the simulations while holding the
grid spacing constant at Δx=0.2.

At the beginning von-Neumann boundary conditions are imposed and the first
value of the potential array is set to one. Then at each time step both
arrays are updated using a simple Euler-step according to
```
dV/dt=η d_xx V+G(V, W)
dW/dt=G(V, W)
```
where G(V, W) is the right hand side in equation (?).

After sufficient steps for the peak of the action potential to develop, the
von-Neumann boundary conditions are replaced by periodic boundary
conditions, which mimics a regular succession of peaks at a constant rate.
This rate can be modified by changing the length of the line along which the
values are computed.

Now, when performing the simulation for different lengths, one finds several
quantities of interest:
 * base cycle length: time for the peak to travel along the line once
   (or the time between two peaks, the period)
 * conduction velocity: velocity of the peak (length of line / base
   cycle length)
 * action potential peak width: time between upstroke and downstroke
 * rise time: duration of the upstroke
 * maximal peak height
The measurements were taken by monitoring a fixed point and checking at each
time step whether the potential crossed either a low or a high threshold (5%
and 95% of the maximal potential respectively).
The length of the line was varied between 25 and 200 length units and for
each value 10 measurements of each aforementioned quantity were taken. From
those 10 measurements, average and standard deviation were computed and
plotted.

> TODO: insert measurement plots

Of course, the base cycle length increases linearly with the maximal length.
The conduction velocity and maximal peak height have smaller
values for small lengths, which corresponds to a high rate of succeeding
action potentials. With longer lengths, i.e. lower rates, those values first
increase strongly in value and then seem to strive to asymptotically against
some upper limit (conduction velocity: 0.0582, max peak: 0.9931).
The action potentials peak width seems to behave similarly, but to verify
the asymptotic behaviour, a closer investigation is required.
The rise time is maximal for small lengths, then falls and goes towards some
lower limit. One can conclude therefrom, that the action potentials rising
edge is the most shallow for high rates and gets steeper for lower rates.

Concerning the theoretical estimated relationship between conduction
velocity and rise time
```
cv≈√(η/(2τ))(1-2α)
```
the measured values were plotted together with the curve given by the
theoretical estimate (where α=0.17 was obtained by fitting the measured data
to the formula using a non-linear least squares method; this is justified,
since α acts only as an offset).

> TODO: insert cv vs rt plot

In order to (naively) quantize this result, both curves were approximated as
linear and a linear regression was computed:

|       | slope   | r^2  |
|---    |---      |---   |
|**exp**| -3.7e-4 | 0.92 |
|**th** | -7.7e-4 | 0.99 |

with a relative error of \~48%.


## 4.2. In 2D
The setups described below were simulated with both the Aliev-Panfilov and
the Fenton model. However, the hereafter presented results are only taken
from Fenton model, since while it exhibits the same phenomena as the
first-mentioned, it also shows some additional interesting features.

The Fenton simulation is based on three two-dimensional arrays (for the
membrane potential V and for the inactivation gates v, w) holding the
spatially resolved values of the corresponding quantities. These arrays are
updated each time step with a simple Euler scheme according to
```
dV/dt = ηΔV + G_V(V, v, w)
dv/dt = G_v(V, v, w)
dw/dt = G_w(V, v, w)
```
where `G_V`, `G_v` and `G_w` are the right hand sides of equations (?), (?)
and (?) respectively.

All of the following setups used grid steps Δx=Δy=1. Varying this parameters
means simply a rescaling of the system.


### 4.2.1. Channel
One considers a 64x256 grid with Dirichlet boundary conditions imposed along
the long side (i.e. the x-direction) and periodical boundary conditions
along the short side (i.e. the y-direction).

For the initial conditions of V and w the same narrow stripe along the y
axis and width of 20 grid points, and for v a stripe of same width but
shifted 10 grid points in positive x-direction were set to 1.

Running the simulation, an action potential peak quickly develops and loops
through the domain (because of the periodic boundaries), which mimics pulses
passing through at a constant rate. Because of the Dirichlet conditions
along the long side regions of the action potential are pulled to zero,
causing some kind of tail.

> TODO: insert channel image

The depicted situation of a confined channel could be used to describe a
bundle of cells belonging to the electrical conducting system of the heart,
which transmits an action potential to some target tissue area, in order to
trigger its contraction.


### 4.2.2. Spiral excitation
Now one considers a 128x512 grid with von-Neumann boundary conditions
imposed on all sides. The initial conditions to generate the first pulse are
the same as described above (section 4.2.1.). Because of the von-Neumann
conditions the pulse is not damped at the borders and passes through the
domain only once.

Now a small area of 10x10 grid points right in the middle of the domain has
its membrane potential set to 1. In the absence of any prior excitation this
simply produces a circular pulse running away in all directions.
But if this perturbation happens in the wake of a preceding action potential
an interaction with the relaxation variables - which have not yet fully
returned to their rest values - takes place: Either the perturbation is
suppressed right away, or it is pushed backwards taking the shape of an
expanding croissant whose edges are running forward until they meet in the
center. Upon meeting again the croissant becomes a pretzel which runs away
in circular fashion, and the then slightly inwards turned edges generate a
new croissant-shaped pulse which undergoes the same process.

> TODO: insert spiral excitation series

For the given setup, the sweet spot to insert this perturbation lies around
2430 integration steps.

While with the Fenton model and the used parameters each succeeding spiral
is slightly damped causing the process to be extinguished completely after a
finite number of iterations, the spirals generated with the Aliev-Panfilov
model do not seem to exhibit such kind of behaviour and instead seem to go
on indefinitely.


### 4.2.3. Spiral wave and breakup
Finally one considers a 256x512 grid with von-Neumann boundary conditions on
all sites.

<!--
TODO: (using Fenton model)
explain setup and discuss results for:
 * channel
 * spiral excitation
 * spiral wave and breakup
-->


# 6. Discussion
## 6.1. The Models
### 6.1.1. Hodgkin & Huxley
 * very intuitive model
 * close to the physics, even though the details were unknown at that time
    * gives a good understanding about what is going on
    * even hyperpolarization is included
 * time step used in integration: Δt=0.01ms, reasonable were obtained with
   Δt=0.055ms, but at this point the results already included artifacts due
   to numerical uncertainties
    * better results with larger time steps could be obtained with a
      different integrator (see below)
 * the source current causing repeated depolarization could perhaps be used
   to model pacemaker cells

### 6.1.2. Aliev & Panfilov
 * difficult to interpret
 * not very close to the underlying physics
 * good for spatial tissue simulations

### 6.1.3. Fenton et al
 * a phenomenological model, yet the quantities (_currents_ and _gates_)
   resemble the physics closely
 * Parameters have physical meaning and are designed in such a way that
   allows to influence the models behaviour purposefully
...
 * rather complex compared to Aliev & Panfilov; tissue simulations take
   subjectively longer

## 6.2. The Method


## 6.3. The relation between conduction velocity and rise time
 * While the gross form of the theoretical estimate and the measured data is
   similar, a rather clear deviation should be acknowledged. Quantitatively
   -from the compared slopes- a relative error of 48% stands against the
   validity of the formula.


# A. Appendix
<!--
TODO:
 * briefly explain code structure
 * usage instructions
-->


>  vim: set ff=unix tw=79 sw=4 ts=4 et ic ai :
