mutable struct curv
    c::Function # Curve
    a::Number # Start value  
    b::Number # End value
    w::Function  # Weights for integration
    N::Number # Amount of
    M::Number # amount of elements
    q::Number # element spacing
end  

Base.copy(s::curv) = curv(s.c, s.a,s.b,s.w,s.N,s.M,s.q)

function curv(c::Function,a::Number,b::Number,w::Function,N::Number)
    curv(c,a,b,w,N,0,0)
end

function set_precision(x::Vector)
    BigFloat.(x,256)
end
function set_precision(x::Real)
     BigFloat.(x,256)
end
function set_precision(x::Complex)
    set_precision.(real.(x)) + 1im * set_precision.(imag.(x))
end

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