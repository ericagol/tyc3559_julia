7/21/2016

Okay, so this was my first successful simulation!  I ran for 10^4 steps, and
then restarted & ran for another 10^4, which looked converged.

I saved these in two JLD files, and the routine post_process.jl computes
their properties.

To make the corner plot, I had to call corner.jl from julia, and then
call_corner.py from python.

To reproduce these results:

1). Make sure Julia & required packages are installed.  The packages used
  are JLD, PyPlot (this requires python & matplotlib to be installed).

2). Start up Julia.  Run the script:

```
julia>  include("fit_tyc3559_norv.jl")
```

3). Let it run - takes about 5 hours on 2015 Macbook Pro.

4). Next, restart & run:

```
julia>  include("fit_tyc3559_restart.jl")
```

  This takes another 5 hours.

5). Plot various diagnostics by running post_process.jl:
   include("post_process.jl")

6). Create a file to read into Python:

```
julia>  include("corner.jl")
```

7). From command line

```
$   python call_corner.py
```
creates a plot named "demo.png".

Julia code file descriptions:

```

fit_tyc3559_norv.jl:				Computes markov chain for TYC 3559 fit to Kepler with correlated noise (no RV data).

fit_tyc3559_restart.jl:				Restarts initial run to compute final 10^4 step chain.

kepler.jl:					Solves Kepler's equation.

occultquad.jl:					Compoutes a quadratic limb-darkend transit/eclipse/occultation light curve.

cgs.jl:						Defines CGS (centimeter-gram-second) constants.

lorentz_likelihood_hermitian_band_init.jl:	Compute covariance matrix & do initial decomposition.

post_process.jl:				Post-processes the markov chain to compute correlation length, covariance, etc.

compute_likelihood_save.jl:			Computes likelihood with saved decomposed matrix.

lorentz_likelihood_hermitian_band_save.jl:	Solve inverse covariance for Lorentzian ACF with saved decomposed matrix.

transit_orb.jl:					Computes transit model for an orbiting body.

corner.jl:					Takes results from markov chain & writes to .csv file for reading into Python to make corner plot.

model_tyc3559_norv.jl:				TYC 3559 photometric model (with no RV fitting).

model_tyc3559_v2.jl:				TYC 3559 photometric model with new parameterization (v2).

bandec.jl:					Complex version of numerical recipes LU decomposition for banded matrix routine.

banbks.jl:					Complex version of numerical recipes LU back-substitution for banded matrix routine.
```

Python code:

```
call_corner.jl:					Calls corner.py to make plot for the paper.
```
