# EAD_Numerical_Simulations

## Notes

This is a repository containing the code to be able to simulate the 2D behaviour of a electroaerodynamic thruster.

Since it was first implemented in a relatively old version of MATLAB, it uses a relatively poorly written MATLAB FEM solver. There is documentation on MATLAB's website, but it takes a while to find and understand.

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

## Code flow:

`driver script.m` - define a problem or a set of problems

For each problem:

-- TRY:

-- -- `convergent_solver()` on the problem

-- -- -- generate geometry

-- -- -- use fzero to find the required rho0 at the emitter to produce the required E field

-- -- -- solve the problem with that rho0 using `get_results_using_rho0()` and return the results.

-- Except:

-- -- note the error, continue
