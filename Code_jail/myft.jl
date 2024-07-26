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

function Levin_ft(f,ksa,s)
    if ksa isa Number
        ks = [ksa];
    else 
        ks = ksa;
    end
    a = s.a
    b = s.b
    sp = Ultraspherical(0.0,UltraMappedInterval(a,b,1.0)); 
    sp1 = Ultraspherical(1.0,UltraMappedInterval(a,b,1.0));
    F = BasisExpansion(x -> f.(x),sp1)
    D = Derivative();
    #tf = true
    res = Complex.(0 .* ks)
    for i1 in eachindex(ks)
        k = ks[i1]
        Op = D - 1im * k* Conversion(sp1)
        if abs(k) < 30
            #"with boundary" |> display
            #bdry = Truncation(FixedGridValues([b],ChebyshevMappedInterval(a,b)) |> Conversion, 4);
            #u = ((bdry ⊘ Op)*SP)\[[0.0]; F]
            #u = \((Op)*SP,F,length(F.c))
            #res1 = u(b)*exp(-1im*k*b) - u(a)*exp(-1im*k*a)
            t1 = time();
            g = x -> exp.(-1im * k .* x) .* f.(x)
            res2 = Clen_Curt(g,s)
            dt = time() - t1;

            #if abs(res2 - res1) <= 1e-15
            #    tf = false
            #end
            res[i1] = res2
        else
            u = \((Op)*sp,F,length(F.c))
            res[i1] = u(b)*exp(-1im*k*b) - u(a)*exp(-1im*k*a)
            #bdry = Truncation(FixedGridValues([b],ChebyshevMappedInterval(a,b)) |> Conversion, 4);
            #u = ((bdry ⊘ Op)*SP)\[[0.0]; F]
            if isnan(res[i1]) || abs(res[i1])>1e2
                g = x -> exp.(-1im * k .* x) .* f.(x)
                res[i1]= Clen_Curt(g,s)
            end
        end
    end
    if ksa isa Number
        res = res[1];
    end
    return res
end