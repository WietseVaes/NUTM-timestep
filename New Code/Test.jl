using Plots, OperatorApproximation
include("Cheb_NUTM.jl")
using .Cheb_NUTM



w =  [1];
L = 2;
t = 1e-5;

f = x -> sin.(pi*x)

gd = UltraMappedInterval(0.0,L,0.0);
sp = Ultraspherical(0.0,gd);

fL = t -> t.^0 - 1;
fR = t -> t.^0 - 1;
ff = BasisExpansion(f,sp);
x = 0.75;
@time f_res = Heat_eq(w,x,t,ff,fL,fR,L)
realsol = (x,t) -> exp(-π^2*t) .* f.(x);
f_res - realsol(x,t)

