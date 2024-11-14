using Plots, OperatorApproximation
include("../Code_jail/Master_func.jl")
using .Master_func

w =  [1];
L = 2;
t = 1e-5;

ff = x->x.^0;

ll = 25-1
ff = x-> cos.(ll*acos.(x-1))*sqrt(2);

gd = UltraMappedInterval(0.0,L,0.0);
sp = Ultraspherical(0.0,gd);

fL = t -> t.^0 - 1;
fR = t -> t.^0 - 1;
fff = BasisExpansion(ff,sp,ll+1)
@time f_res = Heat_eq(w,t,fff,fL,fR,L);

plot(f_res)


f = x -> exp.(-(x-1).^2*40);
f = x -> sin.(pi*x)
F = BasisExpansion(f,sp,25);


w =  [1];
L = 2;
t = 1e-5;
f = x-> sin(π*x);
realsol = (x,t) -> exp(-π^2*t) .* f.(x);
fL = t -> exp.(- pi^2 .* t) .* f(0);
fR = t -> exp.(- pi^2 .* t) .* f.(L);

x = .75;

Rsol = realsol.(x,t);

@time Heat_eq(w,x,t,f,fL,fR,L) - Rsol

x = collect(0.1:0.01:1.9);
res = Heat_eq(w,x,t,f,fL,fR,L)
plot(x,abs.(res-realsol.(x,t)))

start = -1;stop = 1;
g =  z-> (1im*z)^2;
w = [0,1/2,-1/2,1im];
x = 1;
t = 1e-5;
inp = Input_sf(w,x,t,start,stop);
integrand, DD =  Integrand(inp, g);
DomainPlot(DD,x,t)
PathPlot(DD,integrand)




using Plots
#using OperatorApproximation
include("../Code_jail/Master_func.jl")
using .Master_func

f = x -> sin(π*x);
L = 2;
x = 0.1:0.01:1.9;
dt = 1e-5;
w = [1];
realsol = (x,t) -> exp(-π^2*t) .* f.(x);
fL = t -> exp.(- pi^2 .* t) .* f(0);
fR = t -> exp.(- pi^2 .* t) .* f.(L);
fL2 = t -> exp.(- pi^2 .* (t+dt)) .* f(0);
fR2 = t -> exp.(- pi^2 .* (t+dt)) .* f.(L);

Rsol = realsol.(x,dt);
Rsol2 = realsol.(x,2*dt);

t1 = time();
Res = Heat_eq(w, dt, f, fL, fR, L);
Res2 = Heat_eq(w, dt,x-> Res.(x), fL2, fR2, L);
plot(x, abs.(Res.(x)-Rsol), label = "", yaxis =:log);
plot(x, abs.(Res2.(x)-Rsol2), label = "", yaxis =:log);
maximum(abs.(Res.(x)-Rsol));
dt = time()-t1;


