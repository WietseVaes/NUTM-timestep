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
    
    g = x-> exp.(-1im .* k .* x);
    res = k .* 0
    for i1 = 0:Nd
        res += (-1. * 1im)^(i1-1) .* (g(1) * f(1) - g(-1) * f(-1)) ./ ((k) .^ (i1+1))
        f = Diff(f)
    end
    s = curv(x->x,-1,1,x->1,200)
    res += (-1im)^(Nd+1) .* Clen_Curt(x -> f.(x) .* g.(x),s) ./ ((k) .^ (Nd+1))
   return res 
    
end