using LinearAlgebra, Plots

# Morituri te salutant

function OrderRoots(x)
    angles = angle.(x)
    sorted_indices = sortperm(angles)
    return x[sorted_indices]
end

ω(w) = (z) -> sum([w[i] .* z.^(i + 1) for i in 1:length(w)])
Dω(w) = (z) -> sum([(i + 1) .* w[i] .* z.^(i) for i in 1:length(w)])
DDω(w) = (z) -> sum([(i + 1) * i .* w[i] .* z.^(i-1) for i in 1:length(w)])

Φ(w, x, t, k) = begin
    ω_kt = ω(w)(k .* t.^(-1 / (length(w) + 1)))
    return 1im*k .- 1im .* ω_kt .* t ./ x 
end

DΦ(w, x, t, k) = begin
    Dω_kt = Dω(w)(k .* t.^(-1 / (length(w) + 1))) .* t.^(-1 / (length(w) + 1))
    return 1im .- 1im .* Dω_kt .* t ./ x 
end

DDΦ(w, x, t, k) = begin
    DDω_kt = DDω(w)(k .* t.^(-1 / (length(w) + 1))) .* t.^(-2 / (length(w) + 1))
    return - 1im .* DDω_kt .* t ./ x 
end

P(w, x, t, k) = begin
    result = exp(x * Φ(w, x, t, k))
    if result ≈ 0.0
        return 0.0
    else
        return result
    end
end

function k0(w, x, t)
    m = length(w);
    a = ones(m)
    a[1] = x/(-(m+1)*w[m])
    for i1 = 1:(m-1)
        a[i1+1] = (i1+1)*w[i1]*t^((m-i1)/(m+1))/ ((m+1)*w[m] )
    end
    A = zeros(m,m);
    for i1=1:m-1
        A[i1+1,i1] = 1.0;
    end
    A[:,m] = -1. .* a
    E = A |> eigen
    result = E.values |> filter(k -> imag(k) >= -1e-15)
    return OrderRoots(result)
end

function ArgD(w, x, t)
    arg = angle.(DDΦ(w, x, t, k0(w, x, t)))
    return [arg]
end

function mod_offset(x, m, offset)
    return offset .+ mod.(x .- offset, m)
end

function Dirs(w, x, t)
    argd = ArgD(w, x, t)
    return mod_offset(-argd[1] ./ 2 .+ π / 2, π, -π / 2)
end

function ConnectPts(R, w)
    if sign(w[end]) < 0
        return [exp(1im * (4 * m + 1) / (2 * (length(w) + 1)) * π) * R for m in 0:length(w)]
    else
        return [exp(1im * (4 * m - 1) / (2 * (length(w) + 1)) * π) * R for m in 0:length(w)]
    end
end

function SomeConnectPts(R, w)
    if sign(w[end]) < 0
        return [exp(im * (4 * m + 1) / (2 * (length(w) + 1)) * π) * R for m in 0:floor(length(w)/2.)]
    else
        return [exp(im * (4 * m - 1) / (2 * (length(w) + 1)) * π) * R for m in 0:floor(length(w)/2.) + 1]
    end
end

function Min2(X, Y)
    out = []
    for x in X
        temp = abs.(Y .- x)
        p1 = argmin(temp)
        temp[p1] = Inf
        p2 = argmin(temp)
        push!(out, (Y[p1], Y[p2]))
    end
    return out
end

function GlobalR(w, x, t)
    return (30 + abs(x)^2 + 1000 / length(w))^(1 / (length(w) + 1))
end

function Rads(w, xx, tt)
    K0 = k0(w, xx, tt)
    close = Min2(K0, ConnectPts(GlobalR(w, xx, tt), w))
    args = Dirs(w, xx, tt) 
    maxrads = []
    if length(w) > 1
        for i in 1:length(close)
            θ = angle.(close[i][1])
            ϕ = θ - args[i]
            γ = K0[i] * exp(-1im * args[i])
            s = imag(γ) / imag(exp(1im * ϕ))
            r1 = s * exp(1im * ϕ) - γ
            θ = angle.(close[i][2])
            ϕ = θ - args[i]
            γ = K0[i] * exp(-1im * args[i])
            s = imag(γ) / imag(exp(1im * ϕ))
            r2 = s * exp(1im * ϕ) - γ
            push!(maxrads, minimum(abs.([r1, r2])))
        end
    else
        maxrads = [Inf]
    end
    scalerads = [10 / sqrt(abs(DDΦ(w, xx, tt, k)) * abs(xx / 2)) for k in K0]
    
    return K0, args, close, [minimum(r) for r in zip(maxrads, scalerads)]
end

function SDPaths(w, xx, tt)
    K0, args, close, rads = Rads(w, xx, tt)
    paths = [K0[i] .+ [rads[i] * exp(1im * (args[i] + π)) 0; 0 rads[i] * exp(1im * args[i])] for i in 1:length(K0)]
    return K0, close, paths
end

function sortreal(x)
    return sort(x, by=real)
end

function SortAbsIm(x)
    return sort(x, by = x -> abs(imag(x)))
end

function LEOrder(x, y)
    if abs(real(x[1]) - real(y[1])) < eps()
        return imag(x[1]) < imag(y[1])
    else
        return real(x[1]) < real(y[1])
    end
end

function SortLE(x)
    return sort(x, lt = LEOrder)
end

function SmallXPath(w, x, t)
    p1 = sortreal(SomeConnectPts(0.5, w))
    p2 = sortreal(SomeConnectPts(GlobalR(w, x, t), w))
    if length(w) == 1
        p1 .+= 0.5
    end
    s = []
    for i in 1:length(p1) - 1
        push!(s, [p2[i], p1[i]], [p1[i], p1[i+1]], [p1[i+1], p2[i+1]])
    end
    return s
end

function FullPath(w, xx, tt)
    if abs(xx) < 0.1
        return reverse(SmallXPath(w, xx, tt))
    end
    K0, close, paths = SDPaths(w, xx, tt) # Correct
    connect = []
    flatpath = zeros(2*length(K0),2) .* 1im; # Correct
    for i1 = 1:2:2*length(K0)-1
        flatpath[i1,:] = paths[Int(floor((i1-1.)/2.)+1.)][1,:]
        flatpath[i1+1,:] = paths[Int(floor((i1-1.)/2.)+1.)][2,:]
    end
    cpts = ConnectPts(GlobalR(w, xx, tt), w) # Correct
    for i in 1:length(paths)
        ends = [paths[i][1, 1], paths[i][2, 2]]
        closer = [0., 0.] .* 1im
        if abs(ends[1] - close[i][1]) < abs(ends[2] - close[i][1])
            closer[1] = close[i][1]
        else
            closer[1] = close[i][2]
        end
        if abs(ends[2] - close[i][1]) < abs(ends[1] - close[i][1])
            closer[2] = close[i][1]
        else
            closer[2] = close[i][2]
        end
        
        push!(connect, flatpath[2*i-1,:], flatpath[2*i,:], [closer[1], ends[1]], [ends[2], closer[2]])
    end
    
    return connect
end

function plot_paths(path)
    plot()
    for i1 in 1:length(path)
        plot!([real(path[i1][1]),real(path[i1][2])],[imag(path[i1][1]),imag(path[i1][2])], arrow=true, color =:black, linewidth =2, label ="")
        scatter!([real(path[i1][1]),real(path[i1][2])],[imag(path[i1][1]),imag(path[i1][2])], color =:black, markersize =5, label ="")
    end
    plot!()|>display
end