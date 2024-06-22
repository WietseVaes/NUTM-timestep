include("myquad.jl")
include("Cheb.jl")
using OperatorApproximation

function naive_ft(f,s,k)
    
    g = x -> exp.(-1im * k .* x) .* f.(x)
    res = Clen_Curt(g,s)
    
    return res
    
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

function Levin_ft(f,ks,s)
    if ks isa Number
        ks = [ks];
    end
    a = s.a
    b = s.b
    SP = Ultraspherical(0.0,UltraMappedInterval(a,b,1.0)); 
    SP1 = Ultraspherical(1.0,UltraMappedInterval(a,b,1.0));
    F = BasisExpansion(f,SP1)
    D = Derivative();
    tf = true
    res = Complex.(0 .* ks)
    for i1 in eachindex(ks)
        k = ks[i1]
        Op = D - 1im*k* Conversion(SP1)
        if tf
            u = \((Op)*SP,F,length(F.c))
            res1 = u(b)*exp(-1im*k*b) - u(a)*exp(-1im*k*a)
            g = x -> exp.(-1im * k .* x) .* f.(x)
            res2 = Clen_Curt(g,s)

            if abs(res2 - res1) <= 1e-15
                tf = false
            end
            res[i1] = res2
        else
            u = \((Op)*SP,F,length(F.c))
            res[i1] = u(b)*exp(-1im*k*b) - u(a)*exp(-1im*k*a)
            #bdry = Truncation(FixedGridValues([b],ChebyshevMappedInterval(a,b)) |> Conversion, 4);
            #u = ((bdry ⊘ Op)*SP)\[[0.0]; F]
        end
    end
    return res
end