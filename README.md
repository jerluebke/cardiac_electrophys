# Numerical Simulation of Cardiac Tissue Electrophysiology

*Note:* Have a look at the `if __name__ == '__main__'` block of each file for
details on how to use the provided code.


## Hodgkin & Huxley
**File:** `hh52.py` 

Implements the currents and gate variables and integrates the four variables
`V, n, m, h` with a simple Euler step.

Configuration is done via the simulation parameters below the function
definitions.

Visualization is provided by the two functions `plot_all` (for action
potential, currents and gates) and `plot_ap` (action potential only)


## Aliev & Panfilov
### Single Cell
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


### 1D Pulse
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


### 2D Pulse
**File:** `pulse2d.py`

**Classes:**
 * `Pulse2d`: inherits from `AlphaBase`
    * `integrate` solves the system on a 2-dim domain

**Functions**:
 * `channel`
 * `spiral_excitation`
 * `spiral_wave`

The provided functions integrate a special case of the system given by initial
and boundary conditions. The results are animated.


## Fenton et al
**File:** `fche02.py`

**Classes:**
 * `FCHE_Base`
 * `FCHE_Single_Cell`
 * `FCHE_2D`

**Functions:**
 * `plot_single_cell`
 * `channel`
 * `spiral_excitation`
 * `spiral_wave`

Configuration ...


>  vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : 
