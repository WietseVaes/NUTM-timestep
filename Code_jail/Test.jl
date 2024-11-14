using Plots, OperatorApproximation
include("Master_func.jl")
using .Master_func

w =  [1];
L = 2;
t = 1e-5;

f = x -> sin.(pi*x)

gd = UltraMappedInterval(0.0,L,0.0);
sp = Ultraspherical(0.0,gd);

fL = t -> t.^0 - 1;
fR = t -> t.^0 - 1;
ff = BasisExpansion(f,sp);
x = 0.25
@time f_res = Heat_eq(w,x,t,ff,fL,fR,L)
realsol = (x,t) -> exp(-π^2*t) .* f.(x);
f_res - realsol(x,t)
-realsol(x,t)







a = -1.0; b = 1.0;
gd = UltraMappedInterval(a,b,1.0);
sp = Ultraspherical(0.0,gd);
f = x -> sin(x)
D = Derivative();
gv = GridValues(gd);
E = Conversion(gv);
M = Multiplication(x -> 1im*100)
Op = E * D -  M * E
u = (Op*sp)\(x -> f.(x))
