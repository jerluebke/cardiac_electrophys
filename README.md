# Numerical Simulation of Cardiac Tissue Electrophysiology

*Note:* Have a look at the `if __name__ == '__main__'` block of each file for
details on how to use the provided code.


## Models
### Hodgkin & Huxley
**File:** `hh52.py` 

Implements the currents and gate variables and integrates the four variables
`V, n, m, h` with a simple Euler step.

Configuration is done via the simulation parameters below the function
definitions.

Visualization is provided by the two functions `plot_all` (for action
potential, currents and gates) and `plot_ap` (action potential only)


### Aliev & Panfilov
#### Single Cell
**File:** `alpha.py`

**Classes:**
 * `AlphaBase`: handles the model's parameters and provides the RHS of its
   equations
 * `Alpha`: inherits from `AlphaBase`
    * integrates equations in time
    * plots unscaled phase portrait of `V` and `W`
    * plots physical rescaled temporal evolution of `V` and `W`
    * provides method to integrate system with varying one parameter and
      plotting the result

Configuration is done via `PARAMS` dict below definitions.

Visualization is done by calling the appropriate method of an `Alpha` instance.


#### 1D Pulse
**File:** `pulse1d.py`

**Classes:**
 * `Pulse1d`: inherits from `AlphaBase`
    * `integrate` solves the system's equations on a 1-dim domain and records
      various properties of the action potential

**Functions:**
 * `collect_measurements`: perform simulation with `Pulse1d` for given values
   of given parameter and return results in dict
 * `measurements_for_varying_xmax`: for convenience; plot results
 * `animation_of_propagating_action_potential`


#### 2D Pulse
**File:** `pulse2d.py`

**Classes:**
 * `Pulse2d`: inherits from `AlphaBase`
    * `integrate` solves the system on a 2-dim domain and provides data for
      visualization

**Functions**:
 * `channel`
 * `spiral_excitation`
 * `spiral_wave`

The provided functions integrate a special case of the system given by initial
and boundary conditions. The results are animated.


### Fenton et al
**File:** `fche02.py`

**Classes:**
 * `FCHE_Base`: handles parameters, provides currents and RHS of the model's
   equations (`V, v, w`)
 * `FCHE_Single_Cell`: inherits from `FCHE_Base`
    * `integrate` solves `V, v, w` in time
 * `FCHE_2D`: inherits from `FCHE_Base`
    * `integrate` solves `V, v, w` spatially on a 2-dim domain and provides
      data for visualization

**Functions:**
 * `plot_single_cell`
 * `channel`
 * `spiral_excitation`
 * `spiral_wave`

The provided functions integrate a special case of the system given by initial
and boundary conditions and visualize the results. `plot_single_cell` plots
action potential, currents and gates for all four parameter sets; the other
functions take the desired parameter set as input and produce a spatial animation.

Configuration is done via the `PARAM_SETS` dict, which provides four parameter
sets from [Fenton et al, 2002].

Parameter sets used for the presented results:

| Setup | Parameter set |
| ----- |:-------------:|
| single cell       | 1 |
| channel           | 4 |
| spiral excitation | 4 |
| spiral breakup    | 1 |


---

## Usage
### Simulations
All necessary code is given in the `python` files. Make your configurations in
the `if __name__ == '__main__'` block of the appropriate file and call
```
$ cd python/
$ python3 ./file.py
```

or do

```
$ cd python
$ ipyhton

In []: %load file.py
```


### Latex
The files can be compiled with the perl script [latexmk](https://mg.readthedocs.io/latexmk.html):

```
$ cd latex
$ latexmk -xelatex -outdir=./latex_build
```


---

## Technical Notes
### Implementation of Boundary Conditions
The boundary conditions are imposed during the spatial simulations when
computing the laplacian. Consider the 1-dim case:

```python
def Lap1d(arr, out, dx, boundary='neumann'):
    out[:] = (np.roll(arr, 1) + np.roll(arr, -1) - 2*arr) / dx**2
    if boundary == 'neumann':
        out[0], out[-1] = out[1], out[-2]
    elif boundary == 'dirichlet':
        out[0], out[-1] = 0., 0.
    elif boundary == 'periodic':
        # do nothing, already imposed by `np.roll`
        pass
    return out
```


>  vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : 
