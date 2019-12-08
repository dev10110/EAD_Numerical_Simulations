# EAD_Numerical_Simulations

## Notes

This is a repository containing the code to be able to simulate the 2D behaviour of a electroaerodynamic thruster.

Since it was first implemented in a relatively old version of MATLAB, it uses a relatively poorly written MATLAB PDE solver. There is documentation on MATLAB's website, but it takes a while to find and understand.

There are also a number of files which I (Dev) don't know the purpose of. Below I describe the basic architecture to help the new user to get started. The original author included a pdf called `code_desc_v1` which is extremely useful to understand the underlying theory and the implementation.

I heavily modified the interface for the original code, wrapping a lot of it into functions. This was mainly done to be able to enclose the solver into a separate namespace, and therefore allow me to create a single script that can simulate a lot of cases in a for loop.

## Most Useful Functions:

*function `convergent_solver()`:*

  Takes in a single problem (geometry + electrode potentials, initial guess on emitter charge density rho0), and returns the solution (currents, thrust, voltage at each x, y location, etc). Uses a root finder to solve for the actual charge density at the emitter (`find_emitter_eField()`) and then solves the full problem using

  Note on defining initial guess: if given a single scalar, it uses a newton-esque method to root find, but this may (in some cases lead to negative or really large rho0 and it fails). it is better to define a two element vector, with a [min_rho0, max_rho0] that bound the problem, and then it uses a bracketting method to find the actual root. We know the E field as a function of rho0 is monotonously decreasing, and so this helps us.  

*function `find_emitter_eField()`:*

  Takes in a guessed value for rho0, solves the problem using `get_results_using_rho0()`, and finds the e field around the emitter. It then returns the difference between the target e field (currently using Peek's formula). This is basically a helper function for the fzero root finder.

*function `get_results_using_rho0()`:*

  Takes in a specified value for rho0, and solves the problem. It returns a solution structure that contains all the values of interest.

*script `driver_script.m`:*

  Allows user to define a set of problems to be solved, including the initial guesses. The try-except system allows MATLAB to carry on even if one of the solves fails for some reason (usually bad guess on initial rho0)

## Code flow for `driver_script.m`

`driver script.m` - define a problem or a set of problems

For each problem:

-- TRY:

-- -- `convergent_solver()` on the problem

-- -- -- generate geometry

-- -- -- use fzero to find the required rho0 at the emitter to produce the required E field

-- -- -- solve the problem with that rho0 using `get_results_using_rho0()` and return the results.

-- Except:

-- -- note the error, continue

## Code flow for `get_results_using_rho0()`

*Args:*

  rho0: space charge density at emitter

  p, e, t: mesh

  geometry: struct defining the problem geometry

  constants: struct defining constants

  plot_dedug: boolean, if true plots are shown, if disabled, plots are hidden.

*Returns:*

  results: struct containing the solution:

*Code Flow:*

define some useful constants

solve for potential field with 0 space charge

while not converged:

-- find e field at each x, y

-- create the starting location of the characteristic lines

-- propagate the lines to the boundaries

-- while the characteristics don't touch all the boundaries:

-- -- double the number of characteristics and solve again

-- interpolate the solution

-- solve e field at each emitter and collector

-- find currents at emitter and collector (two methods)

-- resolve for potential fields considering the space charge

-- compute convergence test

-- compute the thrust and powers (this part could potentially be placed outside the convergence loop)

store everything into results structure

return results

end



# PDE solvers used:
The PDE solver used is MATLAB's built in PDESolver, `pdenonlin` (documentation)[https://uk.mathworks.com/help/pde/ug/pdenonlin.html]

In the documentation it states that this is no-longer recommended, and that `solvepde` is recommended. Unfortunately, I did not have the time to refactor the code thoroughly to allow this function to be used, but might be recommended.

# Future steps:
I've done a fair bit of code reorg and commenting that did not exist originally but as such have had to use some ugly hacks and less than optimal `pdenonlin`. A proper reorg using a object-oriented method (like in Python or Julia) might be recommended to make the problem more general and
development easier. The code also currently breaks for grounded emitters and negative potential collectors (if the box is ground), and I'm not 100% sure why. Maybe its easier to define the surrounding box and emitter at a high potential and leave the collector at ground, but the characteristics will not only start at the emitter in this case, and thus the code will have some difficulty. Allowing all surfaces to generate characteristic lines might help with this problem, and may significantly speed up the code.

The code also bakes in some assumptions about symmetry that I did not have the time to fully remove.

The code also spends a lot of time trying to find the right rho0 for a target e crit and potential across the electrodes. It might be easier to define a new driver script that just sweeps through the rho0 for a given geometry, and figures out the potential across the electrodes to meet the peek's e field requirement. And then the solution can (if needed) be interpolated at the require potential differences. This way one sweep through rho0s should give most of the solutions needed, and will be significant speed up to the solution.

Also a better way to log the solutions might be warranted.

Enjoy. Feel free to contact me at devansh [dot] agrawal [at] icloud.com for questions.  
