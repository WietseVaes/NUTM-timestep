
include.(("Misc.jl","Spf.jl","Myquad.jl"))

using OperatorApproximation


function Common_Int(w, x, t, start, stop, ψ, α, f, g, Dg, poles, ms, N, L, PD)

    Res = complex.(zeros(length(x),length(t)));

    α₀ = α(0);
    α₁ = α(1)-α₀;
    SF = (m,x,t) -> SpecialFunction(w, x, t, start, stop, 500, k -> g(m,k), Dg, poles, ms(m))

    D = Derivative(1);    

    for i1 = 0:(N-1)
        Res += (f(0)*SF(i1+1, ψ.(x), t) - f(L)*exp(-1im*α₀*L)*SF(i1+1, ψ.(x) .- α₁*L, t));
        f = D*f;
    end

    #SP_time = time()-t1;

    s = curv(x->x,0,L,x->1,500)
    "pre quad" |> display
    for i1 in eachindex(x)
        for i2 in eachindex(t)

            # if PD[1] == -1
                Res[i1,i2] += int_after_IBP_col(w, x[i1], t[i2], f, k -> g(N,k), ψ, α, L, PD);
                #Res[i1,i2] += int_after_IBP_trunc_2(w, x[i1], t[i2], f, k -> g(N,k), ψ, α, L);
            # else
            #     Res[i1,i2] += Clen_Curt(z -> exp(-1im*α₀*z) .* f.(z).*SF.(N,ψ(x[i1]) .- α₁ * z,t[i2]),s)
            # end

            
            #Res[i1,i2] += int_after_IBP_trunc_2(w, x[i1], t[i2], f, k -> g(N,k), ψ, α, L);

            #Res[i1,i2] += int_after_IBP_col(w,x[i1],t[i2],f,k -> 1 ./ (1im * k) .^ N,x->x,L);

            #Res[i1,i2] += Clen_Curt(z -> exp(-1im*α₀*z) .* f.(z).*SF.(N,ψ(x[i1]) .- α₁ * z,t[i2]),s)

        end
    end
    "Post quad" |>display

    return Res ./ (2*π)

end


function Heat_eq(w::Vector,x::Vector,t::Vector,f_init::Function,fL::Function,fR::Function,L::Number)

    N_ft = 0;
    N_quad = 0;
    N = 4;


    gd = UltraMappedInterval(0.0,L,0.0);
    sp = Ultraspherical(0.0,gd);    
    f =  BasisExpansion(f_init,sp);

    ψ = x -> x;
    α = k -> k;
    g = (m,k) -> 1 ./ (1im*α(k)) .^ m;
    Dg = [ k-> 1im*(α(1)-α(0))];
    poles = [-α(0)/(α(1)-α(0))];
    ms = m -> -m;
    start = -1; stop  = 1;

    Int1 = Common_Int(w, x, t, -1, 1, ψ, α, f, g, Dg, poles, ms, N, L, 1 )

    ψ = x -> x .+ 2*L;
    α = k -> k;
    g = (m,k) -> 1 ./ ((exp.(2im*k*L) .- 1) .* (1im*α(k)) .^ m);
    Dg = [k -> 0];
    poles = [];
    ms = m -> [];
    start = -1+1im; stop = 1+1im;
    
    Int2_1 = Common_Int(w, x, t, start, stop, ψ, α, f, g, Dg, poles, ms, N, L, 0);

    ψ = x -> x;
    α = k -> -k;

    Int2_2 = Common_Int(w, x, t, start, stop, ψ, α, f, g, Dg, poles, ms, N, L, 0);

    ψ = x -> x .- 2*L;
    α = k -> k;
    g = (m,k) -> 1/((1 .- exp.(-2im*k*L)).*(1im*α(k)).^m)
    Dg = [k -> 0];
    poles = [];
    ms = m -> [];
    start = 1-1im; stop = -1-1im;
    
    Int3_1 = Common_Int(w, x, t, start, stop, ψ, α, f, g, Dg, poles, ms, N, L, 0);

    ψ = x -> x .- 2*L;
    α = k -> -k;

    Int3_2 = Common_Int(w, x, t, start, stop, ψ, α, f, g, Dg, poles, ms, N, L, 0);

    Res = Int1 - Int2_1 - Int2_2 - Int3_1 - Int3_2

    wtilde = [w[i1+1] * (-1)^(i1) for i1 = 0:(length(w)-1)]
    Res += bd_cond(w,x,t,fL,4,L) + bd_cond(wtilde,L .- x,t,fR,4,L)

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
    #Left integral
    RHS = (u,f_ft) -> g.(-ξ.(u)) .* f_ft(-ξ.(u)) .* exp.(-real.(ω(w)(-ξ.(u)) .* t))
    F = BasisExpansion(u -> RHS(u,f_ft),sp1,1000);

    Op = -M_Dξinv*D + Conversion(sp1)*M;
    rbdry = FixedGridValues([1],ChebyshevMappedInterval(-1,1)) |> Conversion;
    Q = \(((rbdry ⊘ Op)*sp),[[0.0]; F]);
    Left_int = Q(-1)*exp(1im*(-R*ψ(x)-imag.(ω(w)(-R))*t));
    Left_int + Center_int + Right_int
end

function int_after_IBP_trunc_2(w,x,t,f,g,ψ,α,L)

    n = length(w);

    s_ft = curv(x->x,0,L,x->1,0)
    f_ft = k -> Levin_ft(f,α(k),s_ft);
    # Middle piece
    #R = (ψ(x)/t).^(1/n)
    R = 1
    s_out = curv(t -> R*exp.(1im*(pi-t)),0,π,t -> -1im*R*exp.(1im*(pi-t)),0)
    Center_int = Clen_Curt(z -> f_ft.(z) .* exp.(1im*z*ψ(x)-ω(w).(z)*t).*g.(z),s_out)
    "Center_int"|>display
    Center_int |> display
    InfR = 400.
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
    "Right"|>display
    Right_int|>display
    #Left integral
    sp = Ultraspherical(0.0,UltraMappedInterval(-InfR,-R,1.0)); 
    sp1 = Ultraspherical(1.0,UltraMappedInterval(-InfR,-R,1.0));

    F = BasisExpansion(u -> RHS(u,f_ft),sp1);

    Op = D + Conversion(sp1)*M
    rbdry = FixedGridValues([-InfR],ChebyshevMappedInterval(-InfR,-R)) |> Conversion;
    Q = \(((rbdry ⊘ Op)*sp),[[0.0]; F],2*length(F.c));
    Left_int = Q(-R)*exp(1im*(-R*ψ(x)-imag.(ω(w)(-R))*t))
    "Left"|>display
    Left_int|>display
    Left_int + Center_int + Right_int
end

function int_after_IBP_col(w,x,t,f,g,ψ, α,L,PD)

    n = length(w);

    s_ft = curv(x->x,0,L,x->1,200)
    f_ft = k -> Levin_ft_col(f,k,s_ft);

    # Middle piece(s)
    R = 1;

    Res = 0;
    for i1 = 1:2:length(PD) # Only works for monomials
        "i1:"|>display
        i1 |>display
    #R = (ψ(x)/t).^(1/n)
        θ1 = angle(PD[i1]);
        θ2 = angle(PD[i1+1]);

        s_out = curv(t -> R*exp.(1im*(θ2 * t + (1-t)*θ1)),0,1,t -> (θ2-θ1)*1im*R*exp.(1im*(θ2 * t + (1-t)*θ1)),500)

        intC = Clen_Curt(z -> f_ft.(z) .* exp.(1im*z*ψ(x)-ω(w).(z)*t).*g.(z),s_out)
        "Center_int"|>display
        intC |> display


        InfR = 400.
        gd = UltraMappedInterval(R,InfR,1.0);
        sp = Ultraspherical(0.0,gd); 
        gv = GridValues(gd);
        E = Conversion(gv);
        bdry = FixedGridValues([InfR],ChebyshevMappedInterval(R,InfR)) |> Conversion;
        D = Derivative();

        
        M_fun = u -> 1im*(ψ(x)-imag.(Dω(w)(u)) .*t );
        RHS = (u) -> g.(u) .* f_ft(u) .* exp.(-real.(ω(w)(u) .* t))
        
        # The left an right integral turn out to be treated the same
        ML = Multiplication(u -> exp(1im*θ1)*M_fun.(u*exp(1im*θ1)));
        MR = Multiplication(u -> exp(1im*θ2)*M_fun.(u*exp(1im*θ2)));

        OppL = E*D + ML*E;
        OppR = E*D + MR*E;
        plot(real.(collect(R:10:InfR) .* exp(1im*θ1)),imag.(collect(R:10:InfR) .* exp(1im*θ1)))
        plot!(real.(collect(R:10:InfR) .* exp(1im*θ2)),imag.(collect(R:10:InfR) .* exp(1im*θ2)))|>display
        plot(collect(R:InfR),abs.(RHS(collect(R:InfR)*exp(1im*θ1))), yaxis=:log, title = "Left integrand")|>display
        "A"|>display
        plot(collect(R:InfR),abs.(g.(collect(R:InfR)*exp(1im*θ1))), yaxis=:log, title = "g")|>display
        "B"|>display
        plot(collect(R:InfR),abs.(f_ft(collect(R:InfR)*exp(1im*θ1))), yaxis=:log, title = "ft")|>display
        "C"|>display
        plot(collect(R:InfR),abs.(exp.(-real.(ω(w)(collect(R:InfR)*exp(1im*θ1)) .* t))), yaxis=:log, title = "exp")|>display
        "D"|>display

        plot(collect(R:InfR),abs.(RHS(collect(R:InfR)*exp(1im*θ2))), yaxis=:log, title = "Right integrand")|>display
        "E"|>display

        QL = \(((bdry ⊘ OppL)*sp),[[0.0]; u -> RHS(u*exp(1im*θ1))]);
        QR = \(((bdry ⊘ OppR)*sp),[[0.0]; u -> RHS(u*exp(1im*θ2))]);

        intL = -exp(1im*θ1)*QL(R)*exp(1im*(R*exp(1im*θ1)*ψ(x)-imag.(ω(w)(exp(1im*θ1)*R))*t))
        "Left"|>display
        intL|>display
        intR = -exp(1im*θ2)*QR(R)*exp(1im*(R*exp(1im*θ2)*ψ(x)-imag.(ω(w)(R*exp(1im*θ2)))*t))
        "Right"|>display
        intR|>display
        Res += intC + intR - intL

    end

    Res

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
    res =  BasisExpansion(u,sp,300);

    res.c|>display
    res

end

function Heat_eq(w::Vector,x::Vector,t::Vector,f_init::BasisExpansion,fL::Function,fR::Function,L::Number)

    gd = UltraMappedInterval(0.0,L,0.0);
    sp = Ultraspherical(0.0,gd);    
    f =  BasisExpansion(f_init,sp);

    N = 4;

    ψ = x -> x;
    α = k -> k;
    g = (m,k) -> 1 ./ (1im*α(k)) .^ m;
    Dg = [ k-> 1im*(α(1)-α(0))];
    poles = [-α(0)/(α(1)-α(0))];
    ms = m -> [-m];
    start = -1; stop  = 1;

    Int1 = Common_Int(w, x, t, -1, 1, ψ, α, f, g, Dg, poles, ms, N, L,[-1,1])

    ψ = x -> x .+ 2*L;
    α = k -> k;
    g = (m,k) -> 1 ./ ((exp.(2im*k*L) .- 1) .* (1im*α(k)) .^ m);
    Dg = [k -> 0];
    poles = [];
    ms = m -> [];
    start = -1+1im; stop = 1+1im;
    
    Int2_1 = Common_Int(w, x, t, start, stop, ψ, α, f, g, Dg, poles, ms, N, L,[exp(3im/4*π),exp(1im/4*π)]);

    ψ = x -> x;
    α = k -> -k;

    Int2_2 = Common_Int(w, x, t, start, stop, ψ, α, f, g, Dg, poles, ms, N, L,[exp(3im/4*π),exp(1im/4*π)]);

    ψ = x -> x .- 2*L;
    α = k -> k;
    g = (m,k) -> 1/((1 .- exp.(-2im*k*L)).*(1im*α(k)).^m)
    Dg = [k -> 0];
    poles = [];
    ms = m -> [];
    start = 1-1im; stop = -1-1im;
    
    Int3_1 = Common_Int(w, x, t, start, stop, ψ, α, f, g, Dg, poles, ms, N, L,[exp(-1im/4*π),exp(-3im/4*π)]);

    ψ = x -> x .- 2*L;
    α = k -> -k;

    Int3_2 = Common_Int(w, x, t, start, stop, ψ, α, f, g, Dg, poles, ms, N, L,[exp(-1im/4*π),exp(-3im/4*π)]);

    Res = Int1 - Int2_1 - Int2_2 - Int3_1 - Int3_2

    wtilde = [w[i1+1] * (-1)^(i1) for i1 = 0:(length(w)-1)]
    Res += bd_cond(w,x,t,fL,4,L) + bd_cond(wtilde,L .- x,t,fR,4,L)

end

function Heat_eq(w::Vector,x::Number,t::Number,f_init::BasisExpansion,fL::Function,fR::Function,L::Number)

    Res = Heat_eq(w,[x],[t],f_init,fL,fR,L);

    Res[1]

end

function Heat_eq(w::Vector,x::Vector,t::Number,f_init::BasisExpansion,fL::Function,fR::Function,L::Number)

    Res = Heat_eq(w,x,[t],f_init,fL,fR,L);

    Res[:,1]

end

function Heat_eq(w::Vector,x::Number,t::Number,f_init::BasisExpansion,fL::Function,fR::Function,L::Number)
    Res = Heat_eq(w,[x],[t],f_init,fL,fR,L)
    Res[1,1]

end

function Heat_eq(w::Vector,t::Number,f_init::BasisExpansion,fL::Function,fR::Function,L::Number)

    u = x -> Heat_eq(w,x,t,f_init,fL,fR,L);

    gd = UltraMappedInterval(0.0,L,0.0);
    sp = Ultraspherical(0.0,gd);
    #sp_grid = GridValues(gd);
    #x_grid = sp_grid.GD.grid(30);
    #u_grid = u.(x_grid)
    res =  BasisExpansion(u,sp,length(f_init.c)-1);

    res.c|>display
    res

end