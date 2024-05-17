include.(("Misc.jl","Spf.jl","Cheb.jl","Myquad.jl"))

function fun_transform(f,L, Nmax = 100)
    ψ = x -> (x+1)*L/2;
    F_inb = UltraFun(0,x -> f.(ψ(x)),Nmax)
    N = standardChop(F_inb.c)
    if N == Nmax
        @warn "Initial condition is badly approximated: accuracy may be effected."
    end
    UltraFun(0,F_inb.c[1:N])
end

function First_Int(w,x,t,f,N,L)
    
    Res = complex.(zeros(length(x),length(t)));
    if isa(x,Number)
       Res = 0+0im; 
    end
    SF = (m,x,t) -> SpecialFunction(w, x, t, -m, 500); #SF = ∫exp(ixk-w(k)t)*(ik)^m
    for i1 = 0:(N-1)
        Res += (2/L)^i1*(f(-1)*SF(i1+1,x,t)-f(1)*SF(i1+1,x .- L,t));
        f = Diff(f);
    end
    s = curv(x->x,-1,1,x->1,500)
    if isa(x,Number)
        Res += (2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x .- (L/2*(z+1)),t),s)
    else
        for i1 = 1:length(x)
            for i2 = 1:length(t)
                Res[i1,i2] += (2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x[i1] .- (L/2*(z+1)),t[i2]),s)
            end
        end
    end
    return Res ./ (2*π)
end

function Second_Int(w,x,t,f,N,L)
    Res = complex.(zeros(length(x),length(t)));
    if isa(x,Number)
       Res = 0+0im; 
    end
    g = (m,k) -> 1/((exp.(2im*k*L) .- 1).*(1im*k).^m)
    ms = m -> [];
    poles = []
    Dg = [k->0; k->2im*L*exp.(2im*k*L).*(1im*k); k->2im*L*exp.(2im*k*L).*(1im*k); k->2im*L*exp.(2im*k*L).*(1im*k); k->2im*L*exp.(2im*k*L).*(1im*k); k->2im*L*exp.(2im*k*L).*(1im*k); k->2im*L*exp.(2im*k*L).*(1im*k); k->2im*L*exp.(2im*k*L).*(1im*k); k->2im*L*exp.(2im*k*L).*(1im*k)]
    start = -1+1im;stop = 1+1im;
    SF = (m,x,t) -> SpecialFunction(w, x, t,start,stop,500,k -> g(m,k), Dg, poles, ms(m)); #SF = ∫exp(ixk-w(k)t)*(ik)^m
    for i1 = 0:(N-1)
        Res += (2/L)^i1*(f(-1)*SF(i1+1,x .+ 2*L,t)-f(1)*SF(i1+1,x .+ L,t));
        Res += (-2/L)^i1*(f(-1)*SF(i1+1,x,t)-f(1)*SF(i1+1,x .+ L,t));
        f = Diff(f);
    end
    s = curv(x->x,-1,1,x->1,200)
    if isa(x,Number)
        Res += (2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x .+ 2*L .- (L/2*(z+1)),t),s)
        Res += (-2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x .+ (L/2*(z+1)),t),s)
    else
        for i1 = 1:length(x)
            for i2 = 1:length(t)
                Res[i1,i2] += (2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x[i1] .+ 2*L .- (L/2*(z+1)),t[i2]),s)
                Res[i1,i2] += (-2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x[i1] .+ (L/2*(z+1)),t[i2]),s)
            end
        end
    end
    return Res ./ (2*π)
end

function Third_Int(w,x,t,f,N,L)
    Res = complex.(zeros(length(x),length(t)));
    if isa(x,Number)
       Res = 0+0im; 
    end
    g = (m,k) -> 1/((1 .- exp.(-2im*k*L)).*(1im*k).^m)
    #ms = m -> -1 .* [m+1;1;1;1;1;1;1;1;1];
    #poles = [0;-π/L;π/L;-2*π/L;2*π/L;-3*π/L;3*π/L;-4*π/L;4*π/L]
    ms = m -> [];
    poles = []
    Dg = [k->0; k->2im*L*exp.(-2im*k*L).*(1im*k); k->2im*L*exp.(-2im*k*L).*(1im*k); k->2im*L*exp.(-2im*k*L).*(1im*k); k->2im*L*exp.(-2im*k*L).*(1im*k); k->2im*L*exp.(-2im*k*L).*(1im*k); k->2im*L*exp.(-2im*k*L).*(1im*k); k->2im*L*exp.(-2im*k*L).*(1im*k); k->2im*L*exp.(-2im*k*L).*(1im*k)]
    start = 1-1im;stop = -1-1im;
    SF = (m,x,t) -> SpecialFunction(w, x, t,start,stop,500,k -> g(m,k), Dg, poles, ms(m)); #SF = ∫exp(ixk-w(k)t)*(ik)^m
    for i1 = 0:(N-1)
        Res += (2/L)^i1*(f(-1)*SF(i1+1,x .- 2*L, t)-f(1)*SF(i1+1,x .- 3 * L,t));
        Res += (-2/L)^i1*(f(-1)*SF(i1+1,x .- 2*L, t)-f(1)*SF(i1+1,x .- L,t));
        f = Diff(f);
    end
    s = curv(x->x,-1,1,x->1,200)
    if isa(x,Number)
        Res += (2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x .- 2*L .- (L/2*(z+1)),t),s)
        Res += (-2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x .- 2*L .+ (L/2*(z+1)),t),s)
    else
        for i1 = 1:length(x)
            for i2 = 1:length(t)
                Res[i1,i2] += (2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x[i1] .- 2*L .- (L/2*(z+1)),t[i2]),s)
                Res[i1,i2] += (-2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x[i1] .- 2*L .+ (L/2*(z+1)),t[i2]),s)
            end
        end
    end
    return Res ./ (2*π)
end

function Heat_eq(w::Vector,x::Vector,t::Vector,f_init::Function,L::Number)

    f = fun_transform(f_init,L);

    Res = First_Int(w,x,t,f,1,L) - Second_Int(w,x,t,f,1,L) - Third_Int(w,x,t,f,1,L)

end

function Heat_eq(w::Vector,x::Number,t::Number,f_init::Function,L::Number)

    Res = Heat_eq(w::Vector,x::Vector,t::Vector,f_init::Function,L::Number);

    Res[1]

end

function Heat_eq(w::Vector,t::Number,f_init::Function,L::Number)

    f = fun_transform(f_init,L);

    u = x-> First_Int(w,x,t,f,1,L) - Second_Int(w,x,t,f,1,L) - Third_Int(w,x,t,f,1,L);

    res = fun_transform(u,L,50);
    ψ = x -> 2/L*x .- 1;
    x -> res(ψ(x))

end