using Plots, Revise
#using OperatorApproximation
include("../Code_jail/Master_func.jl")
using .Master_func


w =  [1];
L = 2;
t = 1e-5;
f = x-> sin(π*x)
realsol = (x,t) -> exp(-π^2*t) .* f.(x);
fL = t -> exp.(- pi^2 .* t) .* f(0);
fR = t -> exp.(- pi^2 .* t) .* f.(L);

x = .6;

Rsol = realsol.(x,t);

@time Heat_eq(w,x,t,f,fL,fR,L) - Rsol


start = 1-1im;stop = -1-1im;
g =  z-> (1im*z)^2;
w = [0,0,0,0,1];
x = 2;
t = 1e-5;
inp = Input_sf(w,x,t,start,stop);
integrand, DD =  Integrand(inp, g);
DomainPlot(DD,x,t)
PathPlot(DD,integrand)




using Plots
#using OperatorApproximation
include("../Code_jail/Master_func.jl")
using .Master_func

f = x -> sin(π*x)
Nx = 300;
L = 2
x = LinRange(0,L,Nx);
x[1] = []; x[end] = [];
t = 1e-5;
w = [1];
realsol = (x,t) -> exp(-π^2*t) .* f.(x);
fL = t -> exp.(- pi^2 .* t) .* f(0);
fR = t -> exp.(- pi^2 .* t) .* f.(L);

Rsol = realsol.(x,t);

t1 = time();
Res = Heat_eq(w, t, f, fL, fR, L);
plot(x, abs.(Res.(x)-Rsol), label = "", yaxis =:log)
plot(x, real.(Res.(x)), label = "real")
plot!(x, imag.(Res.(x)), label = "imag")
maximum(abs.(Res.(x)-Rsol))
dt = time()-t1