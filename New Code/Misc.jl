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
