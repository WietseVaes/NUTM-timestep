# NUTM_timestep

On the Numerical Unified Transform Method, uncommonly known as Fokas' method, to approximate non-linear PDEs on a bounded spatial domain [0,L].

At the moment I'm trying to:
1. Complete 'special f.jl' such that x can be a complex number:
    * Find better CP_ent and CP_ext points
    * Deal with large critical points
    * Deal with small critical points
    * Check how critical points change as rotations on x are done.
    * Deal with an infinite amount of poles
4. Write out how to compute the second and third integral. The problem is $\Delta(k)$ in the denominator. One can take the domination expential out and leave with something that is approximately 1 and hopefully the pole are not reached in the deformation. 
