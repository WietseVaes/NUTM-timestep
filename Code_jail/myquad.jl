using LinearAlgebra
include("..\\Code_jail\\Cheb.jl")
function standardChop(coeffs, tol = eps())
    # Set default if TOL is not provided
    if isempty(tol)
        tol = eps()
    end
    
    # Check magnitude of TOL
    if tol >= 1
        return 1
    end
    
    # Make sure COEFFS has length at least 17
    n = length(coeffs)
    if n < 17
        return n
    end
    
    # Step 1: Compute ENVELOPE
    b = abs.(coeffs)
    m = fill(b[end], n)
    for j = n-1:-1:1
        m[j] = max(b[j], m[j+1])
    end
    if abs.(m[1]) < 1e-18
        return 1
    end
    envelope = m ./ m[1]
    
    # Step 2: Find PLATEAUPOINT
    plateauPoint = 0
    j2 = 0;
    for j = 2:n
        j2 = round(Int, 1.25*j + 5)
        if j2 > n
            # No plateau: exit
            return n
        end
        e1 = envelope[j]
        e2 = envelope[j2]
        r = 3 * (1 - log(e1) / log(tol))
        plateau = (abs(e1) < 1e-18) || (e2 / e1 > r)
        if plateau
            plateauPoint = j - 1
            break
        end
    end
    
    # Step 3: Determine CUTOFF
    if envelope[plateauPoint] == 0
        return plateauPoint
    else
        j3 = sum(envelope .>= tol^(7/6))
        j2 = min(j3 + 1, j2)
        envelope[j2] = tol^(7/6)
        cc = log10.(envelope[1:j2])
        cc .= cc + range(0,stop=(-1/3)*log10(tol),length=j2)
        d = argmin(cc)
        return max(d - 1, 1)
    end

end

function adaptive_quad(f,Nmax)
    Ns = round.(Int,(1:6) ./ 6 .*Nmax)
    f_ult = [];
    N = 0;
    for NN in Ns
        f_ult = UltraFun(0,f,NN-1);
        N = standardChop(f_ult.c)
        if  N != NN
            break
        end
    end
    if mod(N,2) != 0
        N += 1 
    end
    return f_ult.c[1:min(N,length(f_ult.c))], min(N,length(f_ult.c))
end

function Gauss_grid_weights(a,b)
    
    J = SymTridiagonal(a,b[1:(end-1)]);
    
    xgrid, U = eigen(J);
    
    U = transpose(sign.(U[1,:])) .* U
    
    w = abs.(U[1,:]).^2
    
    return xgrid, w, U
    
end

#   The quads
#   ≡≡≡≡≡≡≡≡≡≡≡

function stand_int(f,s)
    
    # function f 
    
    #object s with: curve - c
    #               start value - a
    #               end value - b
    #               curve derivative - w
    
    f_curv = x-> map(f,map(s.c,x))  .* map(s.w,x); # map to real line
    
    Dtrans_func = (s.b-s.a)/2;
    trans_func = x -> (s.b + s.a)/2 + x * (s.b - s.a)/2 ;
    
    g = x -> f_curv(trans_func(x)) .* Dtrans_func; # map to [-1,1]
    return g
end

function Clen_Curt_no_adapt(f,s)
    
    # function f
    
    #object s with: curve - c
    #               start value - a
    #               end value - b
    #               curve derivative - w
    #               number of quadrature points - N
    
    N = s.N

    if mod(N,2) != 0
       N -= 1 
    end

    #fc = UltraFun(0,f,N)
    #NN = length(fc.c)
    #if sum(abs.(fc.c[NN-1:NN])) > 1e-12
        #@warn "The integrand is incorrect or not approximated well enough"
    #end
    
    f_int = stand_int(f,s)
    n = 0:N/2;
    D = 2 * cos.(2* transpose(n) .* n * pi/N)/N;
    D[1,:] = D[1,:] .* 0.5;
    d = [1; 2 ./ (1 .- (2:2:N).^2)];
    w = D * d;
    x = cos.( (0:N) * π / N );
    w = [w;w[length(w)-1:-1:1]];
    res = dot(w,f_int.(x))
    
    return res
end

function Clen_Curt(f,s)
    
    # function f
    
    #object s with: curve - c
    #               start value - a
    #               end value - b
    #               curve derivative - w
    #               number of quadrature points - N

    f_int = stand_int(f,s);

    (f_coeffs, N) = adaptive_quad(f_int,s.N)
    if s.N == N
        #findlast(abs.(f_coeffs) .> 1e-14) |>display
        @warn "Maximal amount of quadrature points insufficient: accuracy may be effected."
    end
    n = 0:N/2;
    D = 2 * cos.(2* transpose(n) .* n * pi/N)/N;
    D[1,:] = D[1,:] .* 0.5;
    d = [1; 2 ./ (1 .- (2:2:N).^2)];
    w = D * d;
    x = cos.( (0:N) * π / N );
    w = [w;w[length(w)-1:-1:1]];
    res = dot(w,f_int.(x))
    
    return res
end

function m_Filon_Clen_Curt(f_o,s)
    
    # function f
    
    #object s with: curve - c
    #               start value - a
    #               end value - b
    #               curve derivative - w
    #               number of segments - M
    #               Spacing factor - q
    #               number of quadrature points - N (per segment)
    
    M = s.M
    N = s.N 
    q = s.q

    if mod(N,2) != 0
       N -= 1 
    end
    
    f = x-> map(f_o,map(s.c,x))  .* map(s.w,x); # map to real line

    mua = s.a * ((1:M) ./ M).^q;
    mub = s.b * ((1:M) ./ M).^q;
    
    resa = 0;
    resb = 0;

    for i1 = 1:(M-1)
        sa = curv(x->x,mua[i1],mua[i1+1],x->1,N,M,q)
        sb = curv(x->x,mub[i1],mub[i1+1],x->1,N,M,q)
        resa = resa + Clen_Curt_no_adapt(f,sa);
        resb = resb + Clen_Curt_no_adapt(f,sb);
    end

    if s.a == 0
        resa = 0;
    end
    if s.b == 0
        resb = 0;
    end
    
    res = resb-resa;
    
    return res
end

function Trap(f_0,s)
    
    #Good for periodic functions
    
    N = s.N 
    
    f = x-> map(f_0,map(s.c,x))  .* map(s.w,x); # map to real line
    
    dx = (s.b-s.a)/(N);
    
    res = (f(s.a)+f(s.b))/2 * dx
    
    for i1 = 1:N-1
        res = dx * f((s.b-s.a) * i1 / N + s.a)
    end
    
    return res
end

function Laguerre_quad(f_o,s)
    
    #
    # For integrals of the the ∫_0^∞ exp(-x)*g(x)dx, here f(x) = exp(x)*g(x)
    # Problem: we need the poly function and it can't go very high
    #
    
    # xL_k = (2k+1)L_k - kL_{k-1} - (k+1)L_{k+1}, k ≥ 1
    
    N = s.N 
    
    f = x -> map(f_o,map(s.c,x)) .* exp(x)  .* map(s.w,x); # map to real line
    
    K = 1:N
    
    a = collect(2 .* K .- 1) 
    b = collect(-K)
    x, w, uu = Gauss_grid_weights(a,b)
    
    res = sum(w .* map(f,x))
    
    return res
end

function two_sided_Laguerre_quad(f_o,s)
    
    #
    # For integrals of the the ∫_0^∞ exp(-x)*g(x)dx, here f(x) = exp(-x)*g(x)
    # Problem: we need the poly function and it can't go very high
    #
    
    # xL_k = (2k+1)L_k - kL_{k-1} - (k+1)L_{k+1}, k ≥ 1
    
    N = s.N 

    f = x -> map(f_o,map(s.c,x)) .* exp(x)  .* map(s.w,x); # map to real line
    f_oppo = x -> map(f_o,map(s.c,-x)) .* exp(x)  .* map(s.w,-x); # map to real line
    
    K = 1:N
    
    a = collect(2 .* K .- 1) 
    b = collect(-K)
    x, w, uu = Gauss_grid_weights(a,b)
    
    res = sum(w .* map(f,x)) + sum(w .* map(f_oppo,x))
    
    return res
end

function Hermite_quad(f_o,s)
    
    #
    # For integrals of the the ∫_{-∞}^∞ exp(-x^2)*g(x)dx, here f(x) = exp(-x^2)*g(x)
    # Problem: Only takes like N = 5, works good, but very variable
    #
    
    # x H_n = 1/2H_{n+1}+nH_{n-1}
    # 1/2H_{n+1} = x H_n-nH_{n-1}
    
    N = s.N 

    f = x -> map(f_o,map(s.c,x)) .* exp(x .^ 2)  .* map(s.w,x); # map to real line

    a = zeros(N);
    b = sqrt.((1:N) ./ 2)
    x, w, uu = Gauss_grid_weights(a,b)
    
    w *= sqrt(pi) # multiply by int of e^(-x^2)

    res = sum(w .* map(f,x))
    
    return res

end