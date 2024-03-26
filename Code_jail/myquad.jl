using LinearAlgebra

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
    
    g = x -> f_curv(trans_func(x)) * Dtrans_func; # map to [-1,1]
    return g
end

function Clen_Curt(f,s)
    
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
        resa = resa + Clen_Curt(f,sa);
        resb = resb + Clen_Curt(f,sb);
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