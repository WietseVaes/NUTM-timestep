
include.(("Misc.jl","Spf.jl","Myquad.jl"))

using OperatorApproximation


function First_Int(w,x,t,f,N,L)
    
    Res = complex.(zeros(length(x),length(t)));

    D = Derivative(1);
    SF = (m,x,t) -> SpecialFunction(w, x, t, -m, 500); #SF = ∫exp(ixk-w(k)t)*(ik)^m
    t1 = time()
    for i1 = 0:(N-1)
        Res += (f(0)*SF(i1+1,x,t)-f(L)*SF(i1+1,x .- L,t));
        f = D*f;
    end
    SP_time = time()-t1;
    "Special function time"|>display
    SP_time|>display

    s = curv(x->x,0,L,x->1,500)

    for i1 in eachindex(x)
        for i2 in eachindex(t)
            
            Res[i1,i2] += int_after_IBP_trunc_2(w,x[i1],t[i2],f,k -> 1 ./ (1im*k) .^ N,x->x,L);

            #Res[i1,i2] += Clen_Curt(z -> f.(z).*SF.(N,x[i1] .- z,t[i2]),s)

        end
    end


    return Res ./ (2*π)
end

ω(w) = (z) -> sum([w[i] .* z.^(i + 1) for i in eachindex(w)])
Dω(w) = (z) -> sum([(i + 1) .* w[i] .* z.^(i) for i in eachindex(w)])

η = u -> exp.(2*exp.(-1 ./ u) ./ (u-1))  

ϕ(l,L) = begin
    H = L/4
    if 0<= l && l <= L/2-H
        return 1
    elseif l > L/2-H && L/2+H > l
        return η((l-(L/2-H))/(2*H))
    else
        return 0
    end
end

function int_after_IBP_map(w,x,t,f,g,ψ,L)
    s_ft = curv(x->x,0,L,x->1,200)
    f_ft = k -> Levin_ft(f,k,s_ft);
    # Middle piece
    R = 1
    s_out = curv(t -> R*exp.(1im*(pi-t)),0,π,t -> -1im*R*exp.(1im*(pi-t)),500)
    Center_int = Clen_Curt(z -> f_ft.(z) .* exp.(1im*z*ψ(x)-ω(w).(z)*t).*g.(z),s_out)

    #R -> infty
    B = R/2;
    ξ = u -> (B*u .- (2*R-B)) ./ (u-1);
    Dξinv = u -> (u-1) .^ 2/(2*(R-B));


    sp = Ultraspherical(0.0,UltraMappedInterval(-1,1,1.0)); 

    f_ftp = k -> Levin_ft(s -> f(s).*map(z -> ϕ(z,L),s),k,s_ft);
    f_ftm = k -> exp.(1im*k*L) .* Levin_ft(s -> f(s).*(1 .- map(z -> ϕ(z,L),s)),k,s_ft);

    RHS = (u,f_ft) -> g.(ξ.(u)) .* f_ft(ξ.(u)) .* exp.(-real.(ω(w)(ξ.(u)) .* t))
    sp = Ultraspherical(0.0,UltraMappedInterval(-1,1,1.0)); 
    sp1 = Ultraspherical(1.0,UltraMappedInterval(-1,1,1.0));
    Fp = BasisExpansion(u -> RHS(u,f_ftp),sp1,1000);
    Fm = BasisExpansion(u -> RHS(u,f_ftm),sp1,1000);
    Mp_fun = u -> 1im*(ψ(x) - imag.(Dω(w)(ξ.(u))) .*t );
    Mm_fun = u -> 1im*(ψ(x) - L - imag.(Dω(w)(ξ.(u))) .*t );
    Mp = Multiplication(Mp_fun);
    Mm = Multiplication(Mm_fun);

    M_Dξinv = Multiplication(BasisExpansion(Dξinv,sp1));
    D = Derivative();
    Opp = M_Dξinv*D + Conversion(sp1)*Mp
    Opm = M_Dξinv*D + Conversion(sp1)*Mm
    rbdry = FixedGridValues([1],ChebyshevMappedInterval(-1,1)) |> Conversion;
    Qp = \(((rbdry ⊘ Opp)*sp),[[0.0]; Fp]);
    Qm = \(((rbdry ⊘ Opm)*sp),[[0.0]; Fm]);
    Right_int = -Qp(-1)*exp(1im*(ψ(x)*R-imag.(ω(w)(R))*t)) - Qm(-1)*exp(1im*((ψ(x) - L)*R -imag.(ω(w)(R))*t)) 

    #Left integral
    RHS = (u,f_ft) -> g.(-ξ.(u)) .* f_ft(-ξ.(u)) .* exp.(-real.(ω(w)(-ξ.(u)) .* t))

    Fp = BasisExpansion(u -> RHS(u,f_ftp),sp1,1000);
    Fm = BasisExpansion(u -> RHS(u,f_ftm),sp1,1000);

    Mp_fun = u -> 1im*(ψ(x) - imag.(Dω(w)(-ξ(u))) .*t );
    Mm_fun = u -> 1im*(ψ(x) - L - imag.(Dω(w)(-ξ(u))) .*t );

    Mp = Multiplication(Mp_fun);
    Mm = Multiplication(Mm_fun);
    
    M_Dξinv = Multiplication(BasisExpansion(u -> -Dξinv(u),sp1));

    D = Derivative();
    Opp = M_Dξinv*D + Conversion(sp1)*Mp
    Opm = M_Dξinv*D + Conversion(sp1)*Mm
    rbdry = FixedGridValues([1],ChebyshevMappedInterval(-1,1)) |> Conversion;
    Qp = \(((rbdry ⊘ Opp)*sp),[[0.0]; Fp]);
    Qm = \(((rbdry ⊘ Opm)*sp),[[0.0]; Fm]);
    Left_int = Qp(-1)*exp(1im*(-R*ψ(x)-imag.(ω(w)(-R))*t)) + Qm(-1)*exp(1im*(-R*(ψ(x) - L) -imag.(ω(w)(-R))*t))

    Left_int + Center_int + Right_int
end

function int_after_IBP_trunc(w,x,t,f,g,ψ,L)


    s_ft = curv(x->x,0,L,x->1,200)
    f_ft = k -> Levin_ft(f,k,s_ft);
    # Middle piece
    R = 1
    s_out = curv(t -> R*exp.(1im*(pi-t)),0,π,t -> -1im*R*exp.(1im*(pi-t)),500)
    Center_int = Clen_Curt(z -> f_ft.(z) .* exp.(1im*z*ψ(x)-ω(w).(z)*t).*g.(z),s_out)
    InfR = 275.
    #integra = k -> exp.(1im*ψ.(k) .* x + ω(w).(k) .* t) .* g(k) .* f_ft(k) 

    f_ftp = k -> Levin_ft(s -> f(s).*map(z -> ϕ(z,L),s),k,s_ft);
    f_ftm = k -> exp.(1im*k*L) .* Levin_ft(s -> f(s).*(1 .- map(z -> ϕ(z,L),s)),k,s_ft);

    RHS = (u,f_ft) -> g.(u) .* f_ft(u) .* exp.(-real.(ω(w)(u) .* t))
    sp = Ultraspherical(0.0,UltraMappedInterval(R,InfR,1.0)); 
    sp1 = Ultraspherical(1.0,UltraMappedInterval(R,InfR,1.0));
    Fp = BasisExpansion(u -> RHS(u,f_ftp),sp1);
    Fm = BasisExpansion(u -> RHS(u,f_ftm),sp1);

    Mp_fun = u -> 1im*(ψ(x) - imag.(Dω(w)(u)) .*t );
    Mm_fun = u -> 1im*(ψ(x) - L - imag.(Dω(w)(u)) .*t );

    Mp = Multiplication(Mp_fun);
    Mm = Multiplication(Mm_fun);

    D = Derivative();
    Opp = D + Conversion(sp1)*Mp
    Opm = D + Conversion(sp1)*Mm
    rbdry = FixedGridValues([InfR],ChebyshevMappedInterval(R,InfR)) |> Conversion;
    Qp = \(((rbdry ⊘ Opp)*sp),[[0.0]; Fp],2*length(Fp.c));
    Qm = \(((rbdry ⊘ Opm)*sp),[[0.0]; Fm],2*length(Fm.c));
    Right_int = -Qp(R)*exp(1im*(ψ(x)*R-imag.(ω(w)(R))*t)) - Qm(R)*exp(1im*((ψ(x)-L)*R - imag.(ω(w)(R))*t)) 

    #Left integral
    sp = Ultraspherical(0.0,UltraMappedInterval(-InfR,-R,1.0)); 
    sp1 = Ultraspherical(1.0,UltraMappedInterval(-InfR,-R,1.0));

    Fp = BasisExpansion(u -> RHS(u,f_ftp),sp1);
    Fm = BasisExpansion(u -> RHS(u,f_ftm),sp1);

    Opp = D + Conversion(sp1)*Mp
    Opm = D + Conversion(sp1)*Mm
    rbdry = FixedGridValues([-InfR],ChebyshevMappedInterval(-InfR,-R)) |> Conversion;
    Qp = \(((rbdry ⊘ Opp)*sp),[[0.0]; Fp],2*length(Fp.c));
    Qm = \(((rbdry ⊘ Opm)*sp),[[0.0]; Fm],2*length(Fm.c));
    Left_int = Qp(-R)*exp(1im*(-R*ψ(x)-imag.(ω(w)(-R))*t)) + Qm(-R)*exp(1im*(-(ψ(x) - L)*R - imag.(ω(w)(-R))*t))

    Left_int + Center_int + Right_int
end

function int_after_IBP_map_2(w,x,t,f,g,ψ,L)
    s_ft = curv(x->x,0,L,x->1,200)
    f_ft = k -> Levin_ft(f,k,s_ft);
    # Middle piece
    R = 1
    s_out = curv(t -> R*exp.(1im*(pi-t)),0,π,t -> -1im*R*exp.(1im*(pi-t)),500)
    Center_int = Clen_Curt(z -> f_ft.(z) .* exp.(1im*z*ψ(x)-ω(w).(z)*t).*g.(z),s_out)
    "center"|>display
    #R -> infty
    B = R/2;
    ξ = u -> (B*u .- (2*R-B)) ./ (u-1);
    Dξinv = u -> (u-1) .^ 2/(2*(R-B));


    sp = Ultraspherical(0.0,UltraMappedInterval(-1,1,1.0)); 

    RHS = (u,f_ft) -> g.(ξ.(u)) .* f_ft(ξ.(u)) .* exp.(-real.(ω(w)(ξ.(u)) .* t))
    sp = Ultraspherical(0.0,UltraMappedInterval(-1,1,1.0)); 
    sp1 = Ultraspherical(1.0,UltraMappedInterval(-1,1,1.0));
    F = BasisExpansion(u -> RHS(u,f_ft),sp1,1000);
    M_fun = u -> 1im*(ψ(x) - imag.(Dω(w)(ξ.(u))) .*t );
    M = Multiplication(M_fun);

    M_Dξinv = Multiplication(BasisExpansion(Dξinv,sp1));
    D = Derivative();
    Op = M_Dξinv*D + Conversion(sp1)*M
    rbdry = FixedGridValues([1],ChebyshevMappedInterval(-1,1)) |> Conversion;
    Q = \(((rbdry ⊘ Op)*sp),[[0.0]; F]);
    Right_int = -Q(-1)*exp(1im*(ψ(x)*R-imag.(ω(w)(R))*t));
    "Right"|>display
    #Left integral
    RHS = (u,f_ft) -> g.(-ξ.(u)) .* f_ft(-ξ.(u)) .* exp.(-real.(ω(w)(-ξ.(u)) .* t))
    F = BasisExpansion(u -> RHS(u,f_ft),sp1,1000);

    Op = -M_Dξinv*D + Conversion(sp1)*M;
    rbdry = FixedGridValues([1],ChebyshevMappedInterval(-1,1)) |> Conversion;
    Q = \(((rbdry ⊘ Op)*sp),[[0.0]; F]);
    Left_int = Q(-1)*exp(1im*(-R*ψ(x)-imag.(ω(w)(-R))*t));
    "Left"|>display
    Left_int + Center_int + Right_int
end

function int_after_IBP_trunc_2(w,x,t,f,g,ψ,L)


    s_ft = curv(x->x,0,L,x->1,200)
    f_ft = k -> Levin_ft(f,k,s_ft);
    # Middle piece
    R = 1
    s_out = curv(t -> R*exp.(1im*(pi-t)),0,π,t -> -1im*R*exp.(1im*(pi-t)),500)
    Center_int = Clen_Curt(z -> f_ft.(z) .* exp.(1im*z*ψ(x)-ω(w).(z)*t).*g.(z),s_out)
    InfR = 275.

    RHS = (u,f_ft) -> g.(u) .* f_ft(u) .* exp.(-real.(ω(w)(u) .* t))
    sp = Ultraspherical(0.0,UltraMappedInterval(R,InfR,1.0)); 
    sp1 = Ultraspherical(1.0,UltraMappedInterval(R,InfR,1.0));
    F = BasisExpansion(u -> RHS(u,f_ft),sp1);

    M_fun = u -> 1im*(ψ(x)-imag.(Dω(w)(u)) .*t );

    M = Multiplication(M_fun);

    D = Derivative();
    Opp = D + Conversion(sp1)*M
    rbdry = FixedGridValues([InfR],ChebyshevMappedInterval(R,InfR)) |> Conversion;
    Q = \(((rbdry ⊘ Opp)*sp),[[0.0]; F],2*length(F.c));
    Right_int = -Q(R)*exp(1im*(R*ψ(x)-imag.(ω(w)(R))*t))

    #Left integral
    sp = Ultraspherical(0.0,UltraMappedInterval(-InfR,-R,1.0)); 
    sp1 = Ultraspherical(1.0,UltraMappedInterval(-InfR,-R,1.0));

    F = BasisExpansion(u -> RHS(u,f_ft),sp1);

    Op = D + Conversion(sp1)*M
    rbdry = FixedGridValues([-InfR],ChebyshevMappedInterval(-InfR,-R)) |> Conversion;
    Q = \(((rbdry ⊘ Op)*sp),[[0.0]; F],2*length(F.c));
    Left_int = Q(-R)*exp(1im*(-R*ψ(x)-imag.(ω(w)(-R))*t))

    Left_int + Center_int + Right_int
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
    
    D = Derivative(1);
    SF = (m,x,t) -> SpecialFunction(w, x, t,start,stop,500,k -> g(m,k), Dg, poles, ms(m)); #SF = ∫exp(ixk-w(k)t)*(ik)^m
    
    for i1 = 0:(N-1)
        Res += (f(0)*SF(i1+1,x .+ 2*L,t)-f(L)*SF(i1+1,x .+ L,t));
        Res += (f(0)*SF(i1+1,x,t)-f(L)*SF(i1+1,x .+ L,t));
        f = D*f;
    end

    s_ft = curv(x->x,0,L,x->1,200)
    f_ft = k -> Levin_ft(f,k,s_ft);

    SF_end = (x,t) -> SpecialFunction(w, x, t,start,stop,500, k -> g(N,k)*f_ft, [], [], []);


    s = curv(x->x,0,L,x->1,200)
    for i1 in eachindex(x)
        for i2 in eachindex(t)
            #Res[i1,i2] += int_after_IBP_trunc(w,x[i1],t[i2],f,k -> g(N,k), x->x .+ 2*L,x->1,L);
            #Res[i1,i2] += int_after_IBP_trunc(w,x[i1],t[i2],f,k -> g(N,k), x->x, x->1, L, -1);

            Res[i1,i2] += Clen_Curt(z -> f.(z).*SF.(N,x[i1] .+ 2*L .- z,t[i2]),s)
            Res[i1,i2] += Clen_Curt(z -> f.(z).*SF.(N,x[i1] .+ z,t[i2]),s)
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

    D = Derivative(1);
    SF = (m,x,t) -> SpecialFunction(w, x, t,start,stop,500,k -> g(m,k), Dg, poles, ms(m)); #SF = ∫exp(ixk-w(k)t)*(ik)^m

    for i1 = 0:(N-1)
        Res += (f(0)*SF(i1+1,x .- 2*L, t)-f(L)*SF(i1+1,x .- 3 * L,t));
        Res += (f(0)*SF(i1+1,x .- 2*L, t)-f(L)*SF(i1+1,x .- L,t));
        f = D*f;
    end

    s = curv(x->x,0,L,x->1,200)


    for i1 in eachindex(x)
        for i2 in eachindex(t)
            Res[i1,i2] += Clen_Curt(z -> f.(z).*SF.(N,x[i1] .- 2*L .- z,t[i2]),s)
            Res[i1,i2] += Clen_Curt(z -> f.(z).*SF.(N,x[i1] .- 2*L .+ z,t[i2]),s)
        end
    end

    return Res ./ (2*π)
end

function bd_cond(w,x,t,fL,N,L)
    T = 2*maximum(t)
    gd = UltraMappedInterval(0.0,T,0.0);
    sp = Ultraspherical(0.0,gd);    
    f =  BasisExpansion(fL,sp,N);

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
    D = Derivative(1);
    for i1 = 0:(N-1)
        Res +=  f(0) .* SF1(2*(i1+1),x,t)
        Res +=  f(0) .* SF2(2*(i1+1),x .- 2*L,t)
        f = D*f;
    end
    return -1 .* Res ./ (π)
end

function Heat_eq(w::Vector,x::Vector,t::Vector,f_init::Function,fL::Function,fR::Function,L::Number)

    gd = UltraMappedInterval(0.0,L,0.0);
    sp = Ultraspherical(0.0,gd);    
    f =  BasisExpansion(f_init,sp);

    N = 4;

    Int1 = First_Int(w,x,t,f,N,L);
    Int2 = Second_Int(w,x,t,f,N,L);
    Int3 = Third_Int(w,x,t,f,N,L);

    Res = Int1 - Int2 - Int3

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

function Heat_eq(w::Vector,x::Number,t::Number,f_init::Function,fL::Function,fR::Function,L::Number)
    Res = Heat_eq(w,[x],[t],f_init,fL,fR,L)
    Res[1,1]

end

function Heat_eq(w::Vector,t::Number,f_init::Function,fL::Function,fR::Function,L::Number)

    u = x -> Heat_eq(w,x,t,f_init,fL,fR,L);

    gd = UltraMappedInterval(0.0,L,0.0);
    sp = Ultraspherical(0.0,gd);
    #sp_grid = GridValues(gd);
    #x_grid = sp_grid.GD.grid(30);
    #u_grid = u.(x_grid)
    res =  BasisExpansion(u,sp,30);

    res.c|>display
    res

end
