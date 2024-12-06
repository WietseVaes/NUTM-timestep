using FastGaussQuadrature, OperatorApproximation, SpecialFunctions
include.(("Cheb.jl","Misc.jl"))

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
    trans_func = x -> (s.b + s.a)/2 .+ x .* (s.b - s.a) ./ 2 ;
    
    g = x -> f_curv(trans_func(x)) .* Dtrans_func; # map to [-1,1]
    return g
end

function Clen_Curt(f,s)

    f_int = stand_int(f,s);

    gd = UltraMappedInterval(-1,1,0.0);
    sp = Ultraspherical(0.0,gd);
    
    if s.N == 0
        f_int =  BasisExpansion(f_int,sp);
    else
        f_int =  BasisExpansion(f_int,sp, s.N);
    end

    f_coeffs = f_int.c
    N = length(f_coeffs)

    w = [ k % 2 == 0 ? 2. /(1. - k^2) * sqrt(2.) : 0. for k in 0:(N-1)]
    w[1] = 2;
    w[end-1:end] ./= 2
    
    return sum(w .* f_coeffs)
end

function m_Filon_Clen_Curt(f_0,s)
    
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
    
    f = x-> map(f_0,map(s.c,x))  .* map(s.w,x); # map to real line

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

function Laguerre_quad(f_0,s)
    
    #
    # For integrals of the the ∫_0^∞ exp(-x)*g(x)dx, here f(x) = exp(x)*g(x)
    # Problem: we need the poly function and it can't go very high
    #
    
    # xL_k = (2k+1)L_k - kL_{k-1} - (k+1)L_{k+1}, k ≥ 1
    
    N = s.N 
    
    f = x -> map(f_0,map(s.c,x)) .* exp(x)  .* map(s.w,x); # map to real line
    
    x, w = gausslaguerre(N);
    
    res = sum(w .* map(f,x))
    
    return res
end

function two_sided_Laguerre_quad(f_0,s)
    
    #
    # For integrals of the the ∫_0^∞ exp(-x)*g(x)dx, here f(x) = exp(-x)*g(x)
    # Problem: we need the poly function and it can't go very high
    #
    
    # xL_k = (2k+1)L_k - kL_{k-1} - (k+1)L_{k+1}, k ≥ 1
    
    N = s.N 

    f = x -> map(f_0,map(s.c,x)) .* exp(x)  .* map(s.w,x); # map to real line
    f_0ppo = x -> map(f_0,map(s.c,-x)) .* exp(x)  .* map(s.w,-x); # map to real line
    
    x, w = gausslaguerre(N);
    
    res = sum(w .* map(f,x)) + sum(w .* map(f_0ppo,x))
    
    return res
end

function Hermite_quad(f_0,s)
    
    #
    # For integrals of the the ∫_{-∞}^∞ exp(-x^2)*g(x)dx, here f(x) = exp(-x^2)*g(x)
    # Problem: Only takes like N = 5, works good, but very variable
    #
    
    # x H_n = 1/2H_{n+1}+nH_{n-1}
    # 1/2H_{n+1} = x H_n-nH_{n-1}
    
    N = s.N 

    f = x -> map(f_0,map(s.c,x)) .* exp(x .^ 2)  .* map(s.w,x); # map to real line

    a = zeros(N);
    b = sqrt.((1:N) ./ 2)
    x, w, uu = Gauss_grid_weights(a,b)
    
    w *= sqrt(pi) # multiply by int of e^(-x^2)

    res = sum(w .* map(f,x))
    
    return res

end

function Levin_quad(f,g,Dg,s)
    a = s.a
    b = s.b
    SP = Ultraspherical(0.0,UltraMappedInterval(a,b,1.0)); 
    SP1 = Ultraspherical(1.0,UltraMappedInterval(a,b,1.0));
    F = BasisExpansion(f,SP1)
    M = Multiplication(x -> Dg(x));
    D = Derivative();
    Op = D + Conversion(SP1)*M
    lbdry = Truncation(FixedGridValues([a],ChebyshevMappedInterval(a,b)) |> Conversion, 4);
    #u = ((lbdry ⊘ Op)*SP)\[[0]; F]
    u = \((Op)*SP,F)
    u(b)*exp(g(b)) - u(a)*exp(g(a))
end