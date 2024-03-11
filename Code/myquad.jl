using LinearAlgebra, Polynomials, SpecialPolynomials

mutable struct curv
    c # Curve
    a # Start value
    b # End value
    w # Weights for integration
end  

#   Extra stuff
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡

function Cheb_x_w(a,b)
    
    J = Matrix(SymTridiagonal(a,b[1:(end-1)]))
    
    xgrid, U = eigen(J);
    
    U = transpose(sign.(U[1,:])) .* U
    
    w = abs.(U[1,:]).^2
    
    return xgrid, w, U
    
end

function three_term_rec(a,b,c,x)
    
    # getting function values for three term recurrence c_{n+1}T_{n+1} = (x-b_{n+1})T_n + a_{n+1} T_{n-1}
    
    N = length(a)
    fold = ones(length(x),1) # T_0(x)
    f = (x .- b[1]) .* fold ./ c[1]; # T_1(x)
    for i1 = 2:(N-1)
        fnew = ((x .- b[i1]) .* f + a[i1] .* fold)./c[i1]
        fold = f
        f = fnew
    end
    return f
end

#   The actual quads
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function stand_int(f,s)
    
    # function f 
    
    #object s with: curve - c
    #               start value - a
    #               end value - b
    #               curve derivative - w
    
    f_curv = x-> map(f,map(s.c,x))  .* map(s.w,x); # map to real line
    
    Dtrans_func = (s.b-s.a)/2;
    trans_func = x -> (s.b + s.a)/2 + x * (s.b - s.a)/2 ;
    
    g = x -> f_curv(trans_func(x)) * Dtrans_func; # map to [-1,1]
    return g
end

function Clen_Curt(f,s,N)
    
    # function f
    
    #object s with: curve - c
    #               start value - a
    #               end value - b
    #               curve derivative - w
    
    # number of quadrature points N
    
    if mod(N,2) != 0
       N -= 1 
    end
    
    f_int = stand_int(f,s)
    n = 0:N/2;
    D = 2 * cos.(2* transpose(n) .* n * pi/N)/N;
    D[1,:] = D[1,:] .* 0.5;
    d = [1; 2 ./ (1 .- (2:2:N).^2)];
    w = D * d;
    x = cos.( (0:N) * π / N );
    w = [w;w[length(w)-1:-1:1]];
    res = sum(map(f_int,x) .* w)
    
    return res
end

function m_Filon_Clen_Curt(f_o,s,N,M,q)
    
    # function f
    
    #object s with: curve - c
    #               start value - a
    #               end value - b
    #               curve derivative - w
    
    # number of segments M
    
    #Spacing factor q
    
    # number of quadrature points N per segment
    
    if mod(N,2) != 0
       N -= 1 
    end
    
    f = x-> map(f_o,map(s.c,x))  .* map(s.w,x); # map to real line

    mua = s.a * ((1:M) ./ M).^q;
    mub = s.b * ((1:M) ./ M).^q;
    
    resa = 0;
    resb = 0;

    for i1 = 1:(M-1)
        sa = curv(x->x,mua[i1],mua[i1+1],x->1)
        sb = curv(x->x,mub[i1],mub[i1+1],x->1)
        resa = resa + Clen_Curt(f,sa,N);
        resb = resb + Clen_Curt(f,sb,N);
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

function Trap(f_0,s,N)
    
    #Good for periodic functions
    
    f = x-> map(f_0,map(s.c,x))  .* map(s.w,x); # map to real line
    
    dx = (s.b-s.a)/(N);
    
    res = (f(s.a)+f(s.b))/2 * dx
    
    res = dx * sum(map(f,(s.b-s.a) .* (1:(N-1))./N .+ s.a))
    
    return res
end

function F_Laguerre_quad(f_o,s,N)
    
    #
    # For integrals of the the ∫_0^∞ exp(-x)*g(x)dx, here f(x) = exp(-x)*g(x)
    # Problem: we need the poly function and it can't go very high
    #
    
    f = x -> map(f_o,map(s.c,x)) .* exp(x)  .* map(s.w,x); # map to real line
    
    eN = zeros(N+2)
    eN[N+1] = 1 
    
    quad_points = roots(Laguerre{0}(eN))
        
    eN[N+1] = 0;eN[N+2] = 1;
    
    LNp1 = Laguerre{0}(eN)
    weights = quad_points./((N+1)^2*map(LNp1,quad_points).^2)
    
    res = sum(weights .* map(f,quad_points))
    
    return res
end

function Laguerre_quad(f_o,s,N)
    
    #
    # For integrals of the the ∫_0^∞ exp(-x)*g(x)dx, here f(x) = exp(x)*g(x)
    # Problem: we need the poly function and it can't go very high
    #
    
    # xL_k = (2k+1)L_k - kL_{k-1} - (k+1)L_{k+1}, k ≥ 1
    
    f = x -> map(f_o,map(s.c,x)) .* exp(x)  .* map(s.w,x); # map to real line
    
    K = 1:N
    
    a = collect(2 .* K .- 1) 
    b = collect(-K)
    x, w, uu = Cheb_x_w(a,b)
    
    res = sum(w .* map(f,x))
    
    return res
end

function two_sided_Laguerre_quad(f_o,s,N)
    
    #
    # For integrals of the the ∫_0^∞ exp(-x)*g(x)dx, here f(x) = exp(-x)*g(x)
    # Problem: we need the poly function and it can't go very high
    #
    
    # xL_k = (2k+1)L_k - kL_{k-1} - (k+1)L_{k+1}, k ≥ 1
    
    f = x -> map(f_o,map(s.c,x)) .* exp(x)  .* map(s.w,x); # map to real line
    f_oppo = x -> map(f_o,map(s.c,-x)) .* exp(x)  .* map(s.w,-x); # map to real line
    
    K = 1:N
    
    a = collect(2 .* K .- 1) 
    b = collect(-K)
    x, w, uu = Cheb_x_w(a,b)
    
    res = sum(w .* map(f,x)) + sum(w .* map(f_oppo,x))
    
    return res
end

function F_Hermite_quad(f_o,s,N)
    
    #
    # For integrals of the the ∫_{-∞}^∞ exp(-x^2)*g(x)dx, here f(x) = exp(-x^2)*g(x)
    # Problem: Only takes like N = 5, works good, but very variable
    #
    
    f = x -> map(f_o,map(s.c,x)) .* exp(x .^ 2)  .* map(s.w,x); # map to real line
    
    quad_points = roots(basis(Hermite,N))
    
    LNm1 = basis(Hermite,N-1)
    weights = 2^(N-1)*factorial(big(N))*sqrt(pi)./(N^2*map(LNm1,quad_points).^2)
    
    res = sum(weights .* map(f,quad_points))
    
    return res
end

function Hermite_quad(f_o,s,N)
    
    #
    # For integrals of the the ∫_{-∞}^∞ exp(-x^2)*g(x)dx, here f(x) = exp(-x^2)*g(x)
    # Problem: Only takes like N = 5, works good, but very variable
    #
    
    # x H_n = 1/2H_{n+1}+nH_{n-1}
    # 1/2H_{n+1} = x H_n-nH_{n-1}
    
    f = x -> map(f_o,map(s.c,x)) .* exp(x .^ 2)  .* map(s.w,x); # map to real line
    
    K = 1:N
    
    upperdiag = ones(N-1) ./ 2
    lowerdiag = collect(1:(N-1))
    
    DD = diagm(1 => upperdiag, -1 => lowerdiag)
    
    invvec = ones(N);
    
    for NN = 1:(N-1)
        invvec[NN+1] = invvec[NN] * sqrt(2*NN)
    end
    vvec = 1. ./ invvec
    Sim_mat = diagm(0 => vvec)
    invSim_mat = diagm(0 => invvec)
    
    J = Sim_mat * DD * invSim_mat
    
    #J - Symmetric(J) |>display
    
    x, U = eigen(J);
    
    U = transpose(sign.(U[1,:])) .* U
    
    w = abs.(U[1,:]).^2*sqrt(pi) # multiply by int of e^(-x^2)
    
    
    #M = N-0
    
    #a = -collect(0:(M-1)) 
    #b = zeros(M)
    #c = ones(M) ./ 2
    
    #Hnm1 = three_term_rec(a,b,c,x)
    
    #Hnm1 = map(basis(Hermite, N-1),x) #Slower
    
    # w = big(2^(N-1))*factorial(big(N))*sqrt(π) ./ (N^2 .* Hnm1 .^ 2) #direct use of factorial 1/thirds the amount of N we can use 
    
    #w = sqrt(π) ./ (N^2 .* Hnm1 .^ 2) / 2
    
    #for i = 1:N
    #    w *= (i*2) #loop doubles the amount of N we can use
    #end

    res = sum(w .* map(f,x))
    
    return res
end

using NBInclude
nbexport("myquad.jl", "myquad.ipynb")