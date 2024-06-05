include.(("Misc.jl","Spf.jl","Cheb.jl","Myquad.jl"))

function fun_transform(f,L, Nmax = 100)
    ψ = x -> (x .+ 1.) .* L/2;
    F_inb = UltraFun(0,x -> f(ψ(x)),Nmax)
    N = standardChop(F_inb.c)
    if N == Nmax + 1
        @warn "Function is badly approximated: accuracy may be effected."
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
        for i1 in eachindex(x)
            for i2 in eachindex(t)
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
        for i1 in eachindex(x)
            for i2 in eachindex(t)
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
        for i1 in eachindex(x)
            for i2 in eachindex(t)
                Res[i1,i2] += (2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x[i1] .- 2*L .- (L/2*(z+1)),t[i2]),s)
                Res[i1,i2] += (-2/L)^(N-1) .* Clen_Curt(z -> f.(z).*SF.(N,x[i1] .- 2*L .+ (L/2*(z+1)),t[i2]),s)
            end
        end
    end
    return Res ./ (2*π)
end

function bd_cond(w,x,t,fL,N,L)
    T = 2*maximum(t)
    f = fun_transform(fL,T);
    Res = complex.(zeros(length(x),length(t)));
    if isa(x,Number)
       Res = 0+0im; 
    end
    g1 = (m,k) -> 1im/((exp.(2im*k*L) .- 1).*(k).^(m-1))
    g2 = (m,k) -> 1im/((1 .- exp.(-2im*k*L)).*(k).^(m-1))
    ms = m -> [];
    poles = []
    Dg = []
    start1 = -1+1im;stop1 = 1+1im;
    start2 = -1+1im;stop2 = 1+1im;
    SF1 = (m,x,t) -> SpecialFunction(w, x, t,start1,stop1,500,k -> g1(m,k), Dg, poles, ms(m)); 
    SF2 = (m,x,t) -> SpecialFunction(w, x, t,start2,stop2,500,k -> g2(m,k), Dg, poles, ms(m)); 
    for i1 = 0:(N-1)
        Res += (-2/T)^(i1) * f(-1) .* SF1(2*(i1+1),x,t)
        Res += (-2/T)^(i1) * f(-1) .* SF2(2*(i1+1),x .- 2*L,t)
        f = Diff(f);
    end
    return -1 .* Res ./ (π)
end

function Heat_eq(w::Vector,x::Vector,t::Vector,f_init::Function,fL::Function,fR::Function,L::Number)

    f = fun_transform(f_init,L);

    Res = First_Int(w,x,t,f,1,L) - Second_Int(w,x,t,f,1,L) - Third_Int(w,x,t,f,1,L)

    wtilde = [w[i1+1] * (-1)^(i1) for i1 = 0:(length(w)-1)]
    Res += bd_cond(w,x,t,fL,4,L) + bd_cond(wtilde,L .- x,t,fR,4,L)

end

function Heat_eq(w::Vector,x::Number,t::Number,f_init::Function,fL::Function,fR::Function,L::Number)

    Res = Heat_eq(w,[x],[t],f_init,fL::Function,fR::Function,L);

    Res[1]

end

function Heat_eq(w::Vector,x::Vector,t::Number,f_init::Function,fL::Function,fR::Function,L::Number)

    Res = Heat_eq(w,x,[t],f_init,fL::Function,fR::Function,L);

    Res[:,1]

end

function Heat_eq(w::Vector,t::Number,f_init::Function,fL::Function,fR::Function,L::Number)

    f = fun_transform(f_init,L);

    u = x -> Heat_eq(w,x,t,f_init,fL::Function,fR::Function,L);

    res = fun_transform(u,L,30);
    res.c|>display
    ψ = x -> 2/L .* x .- 1;
    x -> res.(ψ(x))

end
