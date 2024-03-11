include("myquad.jl")
include("myCheb.jl")

function naive_ft(f,s,k,N)
    
    g = x -> exp.(-1im*k .* x).*map(f,x)
    k_amp = Clen_Curt(g,s,N)
    
    return k_amp
    
end

function IBP_ft(f,k,N,Nd)
    
    # Function f
    
    # Frequency k
    
    # Number of Chebychev polynomials N
    
    # Number of Derivatives Nd
    
    # Easily made faster.
    
    
    x = -1:2:1;
    
    g(x) = exp.(-1im .* k .* x);
    
    res = 1im ./ k .* (g(1) .* f(1) - g(-1) .* f(-1))
    
    for i1 = 1:Nd
        
        uk = Cheb_Der(f,x,i1,N);
        
        res += (-1)^i1 * (1im ./ k) .^ (i1+1) .* (g(1) .* uk[2] - g(-1) .* uk[1])
    end
    
    
   return res 
    
end

using NBInclude
nbexport("myft.jl", "myft.ipynb")