using Plots
#using OperatorApproximation
include("../Code_jail/Master_func.jl")
using .Master_func

L = 2;
s = curv(x->x, 0, L, x->1,400);
f = x -> sin.(pi*x); f_ft = k -> -π*(1im*sin.(2*k)-cos(2*k)+1)/(k^2-π^2);
f = x -> sin.(pi*x).^2; f_ft = k -> -(2 * π^2 * (sin(2 * k) + 1im * cos(2 * k) - 1im)) / (k * (k^2 - 4 * π^2));
R = 2;
B = 1;
kkk = ξ -> (B .* ξ .- 2*R .+ B)./(ξ .- 1);
ff_ft = k -> Levin_ft(f,kkk.(k),s);

dk = .001;
kk = -1:dk:1-dk;
@time ff_ft(kk);
plot(kk,imag.(f_ft.(kkk.(kk))-ff_ft(kk)));plot!(kk,real.(f_ft.(kkk.(kk))-ff_ft(kk)))


f = x -> cos.(π*x);

Nx = 1000;
w =  [1];
L = 2;
t = 1e-5;
x = collect(L/Nx:L/(Nx):L-L/Nx);
#xx = 0:L/(Nx):L;


Nx = length(x);

realsol = (x,t) -> exp(-π^2*t) .* f.(x);
fL = t -> exp.(- pi^2 .* t) .* f(0);
fR = t -> exp.(- pi^2 .* t) .* f.(L);

Rsol = realsol.(x,t);
RRsol = realsol.(x,2*t);


t1 = time();
Res = Heat_eq(w, t, f, fL, fR, L);
plot(x, abs.(Res.(x)-Rsol), label = "", yaxis =:log)
plot(x, real.(Res.(x)), label = "real")
plot!(x, imag.(Res.(x)), label = "imag")
maximum(abs.(Res.(x)-Rsol))
dt = time()-t1

using OperatorApproximation, Plots, SpecialFunctions

R = 30;
sp = Ultraspherical(0.0,UltraMappedInterval(-R,R,2.0)); 
sp2 = Ultraspherical(2.0,UltraMappedInterval(-R,R,2.0));
M = Multiplication(x -> x);
D = Derivative();
Op = D^2 - Conversion(sp2)*M
lbdry = FixedGridValues([-R],ChebyshevMappedInterval(-R,R)) |> Conversion;
rbdry = FixedGridValues([R],ChebyshevMappedInterval(-R,R)) |> Conversion;
setbasis(sp)
u = (lbdry ⊘ rbdry ⊘ Op)\[[airyai(-R)];[airyai(R)]; [0;0]]
u = u[1]
u(-R)-airyai(-R)
plot(u;dx=0.0001)


using Plots
x = 1;
L = 2;
N = 3;
t = 1e-5
g = s -> (pi*L/2)^N*sin(π*L/2 .* (s .+ 1)+N*π/2);

I = (s,k) ->  exp.(1im .* k .* (x .- L/2 .* (s .+ 1)) .- k.^2*t)./(1im.*k).^N .* g.(s)
s = -1:0.1:1
k = -1:.0001:1
Inte = complex.(zeros(length(s),length(k))) 
for i1 in eachindex(k)
    Inte[:,i1] = map(s -> I(s,k[i1]),s);
end
plot(layout = grid(2, 2),size=(1000,800));
surface!(k,s,real.(Inte), camera = (45, 45), xlabel = "k", ylabel = "s",title="real", subplot = 1)
heatmap!(k,s,real.(Inte), xlabel = "k", ylabel = "s",title="real", subplot = 2)

surface!(k,s,imag.(Inte), camera = (45, 45), xlabel = "k", ylabel = "s",title="imag", subplot = 3)
heatmap!(k,s,imag.(Inte), xlabel = "k", ylabel = "s",title="imag", subplot = 4)

k = -100:.01:100
plot(k,abs.(real.(map(k -> I(s[1],k),k))), yaxis =:log)