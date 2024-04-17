# NUTM_timestep

On the Numerical Unified Transform Method, uncommonly known as Fokas' method, to approximate non-linear PDEs on a bounded spatial domain [0,L].

At the moment I'm trying to:
1. Verify if my code works for the first integral
2. Complete special f.jl such that x can be a complex number.
3. Write out how to compute the second and third integral. The problem is $\Delta(k)$ in the denominator. One can take the domination expential out and leave with something that is approximately 1 and hopefully the pole are not reached in the deformation. 
