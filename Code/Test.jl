using Plots
include("../Code_jail/Master_func.jl")
using .Master_func

f = x -> sin(x);

Nx = 1000;
w =  [1];
L = 2*pi;
t = 1e-2;
x = collect(L/Nx:L/(Nx):L-L/Nx);
#xx = 0:L/(Nx):L;


Nx = length(x);


realsol = (x,t) -> exp(-t)*sin(x);
Rsol = realsol.(x,t);
RRsol = realsol.(x,2*t);

t1 = time();
Res = Heat_eq(w, t, f, L);
plot(x, abs.(Res.(x)-Rsol), label = "", yaxis =:log)
maximum(abs.(Res.(x)-Rsol))
Ress = Heat_eq(w, t, Res, L);
dt = time()-t1

plot(x, abs.(Ress.(x)-RRsol), label = "", yaxis =:log)
maximum(abs.(Ress.(x)-RRsol))


Nx = 1000;
w =  [1];
L = 3;
t = 1e-2;
x = collect(L/Nx:L/(Nx):L-L/Nx);
f = x -> -x*(x-L);
Sol = Heat_eq(w, 1, f, L);
plot(x,real.(Sol.(x)))
plot(x,imag.(Sol.(x)))