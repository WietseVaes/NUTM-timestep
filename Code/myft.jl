include("myquad.jl")
include("myCheb.jl")
include("Cheb.jl")

function naive_ft(f,s,k)
    
    g = x -> exp.(-1im * k .* x) .* f.(x)
    k_amp = Clen_Curt(g,s)
    
    return k_amp
    
end

function IBP_ft(f,k,Nd)
    
    # Ultrafun f
    
    # Frequency k
    
    # Number of Chebychev polynomials N
    
    # Number of Derivatives Nd
    
    # Easily made faster.
    
    g = x-> exp.(-1im .* k .* x);
    
    res = 1im * (g(1) * f(1) - g(-1) * f(-1)) ./ k
    
    for i1 = 1:Nd
        Df = Diff(f)
        (-1im)^(i1-1) .* (g(1) * Df(1) - g(-1) * Df(-1)) ./ ((k) .^ (i1+1))|>display
        res += (-1im)^(i1-1) .* (g(1) * Df(1) - g(-1) * Df(-1)) ./ ((k) .^ (i1+1))
        f = Df
    end
    
   return res 
    
end