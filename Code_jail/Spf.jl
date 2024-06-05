using Plots, Printf
include.(("Misc.jl", "Myquad.jl"))

mutable struct Input_sf
    # Dispersion relation coefficients
    w::Vector
    # abs(x)
    x::Number
    # time t
    t::Number
    # angle of x
    xθ::Number
    # Polar form of w[end]
    wr::Number
    wθ::Number
    start::Number
    stop::Number
end
function angle2pi(z)
    ang = angle.(z);
    if ang isa Array{<:Number,1}
        for i1 in eachindex(ang)
            if ang[i1] < -1e-15
                ang[i1] += 2*π
            end
        end
    else
        ang = ang < -1e-15 ? ang + 2π : ang 
    end
    return ang
end
function polar(z)
    return (abs.(z), angle2pi(z))
end
function Input_sf(ww::Vector,xx::Number,tt::Number,start::Number,stop::Number)
    xxr, xxθ = polar(xx);
    wwr, wwθ = polar(ww[end]);
    Input_sf(ww,xxr,tt,xxθ,wwr,wwθ,start,stop)
end

struct Deformation
    pp::Vector
    cc::Vector
    path::Vector
    Dpath::Vector
    tt::Vector
    meth::Vector
    w::Vector
    xθ::Number
    dir::Number
end
function dom_sectioner(path)
    a = -1; b = 1;
    distances = [abs(path[i][2]-path[i][1]) for i in eachindex(path)]; distances .*= (b-a)/sum(distances)
    t_vals = zeros(length(path)+1); t_vals[1] = a; t_vals[end] = b;
    for i1 = 1:length(path)-1
        t_vals[i1+1] = t_vals[i1]+distances[i1]
    end
    return t_vals
end
function pwlinf_maker(t,t_vals,path_nodes,meth)
    for i1 = 1:(length(t_vals)-1)
        if (t_vals[i1] <= t && t <= t_vals[i1+1])
            return (t-t_vals[i1])/(t_vals[i1+1]-t_vals[i1])*path_nodes[i1+1] + (t-t_vals[i1+1])/(t_vals[i1]-t_vals[i1+1])*path_nodes[i1]
        end
    end
    if t_vals[end] <= t
        return (t-t_vals[end-1])/(t_vals[end]-t_vals[end-1])*path_nodes[end] + (t-t_vals[end])/(t_vals[end-1]-t_vals[end])*path_nodes[end-1]
    elseif t <= t_vals[1]
        return (t-t_vals[1])/(t_vals[2]-t_vals[1])*path_nodes[2] + (t-t_vals[2])/(t_vals[1]-t_vals[2])*path_nodes[1]
    end
end
function Dpwlinf_maker(t,t_vals,path_nodes,meth)
    for i1 = 1:(length(t_vals)-1)
        if t >= t_vals[i1] && t <= t_vals[i1+1]
            return 1/(t_vals[i1+1]-t_vals[i1])*path_nodes[i1+1] + 1/(t_vals[i1]-t_vals[i1+1])*path_nodes[i1]
        end
    end
    if t_vals[end] <= t
        return 1/(t_vals[end]-t_vals[end-1])*path_nodes[end] + 1/(t_vals[end-1]-t_vals[end])*path_nodes[end-1]
    elseif t <= t_vals[1]
        return 1/(t_vals[2]-t_vals[1])*path_nodes[2] + 1/(t_vals[1]-t_vals[2])*path_nodes[1]
    end
end
function fpath_maker(path,cate)
    t_vals = dom_sectioner(path)
    meth = [];
    tt = [];
    funcs = [];
    Dfuncs = [];
    for i1 in eachindex(path)
        push!(meth, "Clenshaw-Curtis")
        push!(funcs, t -> pwlinf_maker(t,t_vals[i1:(i1+1)],path[i1],meth[i1]))
        push!(Dfuncs, t -> Dpwlinf_maker(t,t_vals[i1:(i1+1)],path[i1],meth[i1]))
        push!(tt, t_vals[i1:(i1+1)])
    end
    return funcs, Dfuncs, tt, meth
end
function Deformation(path::Vector,cate::Vector,w::Vector,xθ::Number,dir::Number)
    func, Dfunc, tt, meth = fpath_maker(path,cate)
    Deformation(path,cate,func, Dfunc, tt, meth,w,xθ,dir)
end 
function tick_maker(y)
    exponent = Int.(floor.(log10.(abs.(y))))
    if maximum(abs.(exponent)) < 3
        scaled_value = range(y[1],stop=y[2],length = 5)
    else
        RR = collect(range(y[1],stop=y[2],length = 5))
        RR2 = RR .*( 10. .^ (-maximum(exponent)))
        scaled_value = round.(RR2, digits = 2)
    end
    formatted_value = []
    for i1 in eachindex(scaled_value)
        push!(formatted_value, @sprintf("%.2f", scaled_value[i1]))
    end
    return_formatted_value = []
    for i1 in eachindex(scaled_value)
        push!(return_formatted_value,"$(formatted_value[i1]) × 10^$(maximum(exponent))")
    end
    
    # Construct the final label
    if maximum(abs.(exponent)) < 3
        return scaled_value, formatted_value
    else
        return RR, return_formatted_value
    end
end

function DomainPlot(D::Deformation,x::Number,t::Number)
    pl = plot();
    path = D.pp;
    cate = D.cc
    CP_label = "CP";
    inf_label = "inf";
    CP_ext_label = "CP_ext";
    CP_ent_label = "CP_ent";
    CP_label = "";
    inf_label = "";
    CP_ext_label = "";
    CP_ent_label = "";
    Def_label = "Deformation";
    titlestr = @sprintf("Deformation for x = %.2f and t = %.2f", x, t)
    for i1 in 1:length(path)
        plot!([real(path[i1][1]),real(path[i1][2])],[imag(path[i1][1]),imag(path[i1][2])], arrow=true, color =:black, linewidth =2, label = Def_label, title = titlestr);
        Def_label = "";
        for i2 in 1:2
            if cate[i1][i2] == "CP"
                scatter!([real(path[i1][i2])],[imag(path[i1][i2])], color =:green, markersize =5, label = CP_label);
                CP_label = "";
            elseif cate[i1][i2] == "inf"
                scatter!([real(path[i1][i2])],[imag(path[i1][i2])], color =:red, markersize =5, label = inf_label)
                inf_label = "";
            elseif cate[i1][i2] == "CP_ext"
                scatter!([real(path[i1][i2])],[imag(path[i1][i2])], color =:orange, markersize =5, label = CP_ext_label)
                CP_ext_label = "";
            elseif cate[i1][i2] == "CP_ent"
                scatter!([real(path[i1][i2])],[imag(path[i1][i2])], color =:orange, markersize =5, label = CP_ent_label)
                CP_ent_label = "";
            end 
        end
    end

    ww = D.w;
    n = length(ww) + 1;

    xr = xlims(pl);
    yr = ylims(pl);

    # Generate a grid of complex numbers
    xplt = LinRange(minimum([xr[1],yr[1]]), maximum([xr[2],yr[2]]), 500)
    yplt = LinRange(minimum([xr[1],yr[1]]), maximum([xr[2],yr[2]]), 500)
    z = [a + 1im*b for b in yplt, a in xplt]
    xθ = D.xθ
    # Create a mask for points where real((x + yi)^2) < 0
    #mask = [real(-1im*ω(ww)(zi)) >= 0 for zi in z]
    mask = [real(ww[end]*(zi)^n) < 0 for zi in z]

    mmask = [ NaN for zi in z]
    mmask[mask] .= 0
    heatmap!(xplt, yplt, mmask, c=:RdBu,  alpha =0.5, cbar=false,aspect_ratio=:equal, xlims = (minimum([xr[1],yr[1]]), maximum([xr[2],yr[2]])), ylims = (minimum([xr[1],yr[1]]), maximum([xr[2],yr[2]])))

    xθ = D.xθ;
    r = LinRange(minimum([xr[1],yr[1]]),maximum([xr[2],yr[2]]), 500);
    plot!(r*cos(xθ),r*sin(xθ),c=:blue, linewidth=3, label = "Original int", xlims = (minimum([xr[1],yr[1]]), maximum([xr[2],yr[2]])), ylims = (minimum([xr[1],yr[1]]), maximum([xr[2],yr[2]])), legend =:outerbottom);
    pl = plot!();
end
function DomainPlot(D::Deformation)
    DomainPlot(D,0,0)
end
function PathPlot(D::Deformation, integrand::Function)
    tpoint = D.tt 
    pathf = D.path 
    cate = D.cc
    funcs = [z->real.(z); z->imag.(z)]
    titles = ["Real integrand"; "Complex integrand"]
    plts = [];
    for i3 = 1:2
        plt = plot()
        gg = funcs[i3]
        for i1 in 1:length(pathf)
            ttend = tpoint[i1]
            tt = ttend[1]:0.001:ttend[2];
            plot!(tt,gg.(integrand.(pathf[i1].(tt))), color =:black, linewidth =2, label ="",title = titles[i3]);
            for i2 = 1:2
                if cate[i1][i2] == "CP"
                    scatter!([ttend[i2]],[gg(integrand(pathf[i1](ttend[i2])))], color =:green, markersize =5, label = "");
                elseif cate[i1][i2] == "inf"
                    scatter!([ttend[i2]],[gg(integrand(pathf[i1](ttend[i2])))], color =:red, markersize =5, label = "");
                elseif cate[i1][i2] == "CP_ext"
                    scatter!([ttend[i2]],[gg(integrand(pathf[i1](ttend[i2])))], color =:orange, markersize =5, label = "");
                elseif cate[i1][i2] == "CP_ent"
                    scatter!([ttend[i2]],[gg(integrand(pathf[i1](ttend[i2])))], color =:orange, markersize =5, label = "");
                end 
            end
        end
        yr = ylims(plt);
        yticks_positions, yticks_lables = tick_maker(yr);
        plt = plot!(yticks=(yticks_positions, yticks_lables))
        push!(plts, plt)
    end
    pl = plot(plts[1], plts[2], layout = grid(2, 1),size=(800,500));

end 
ω(w) = (z) -> sum([w[i] .* z.^(i + 1) for i in eachindex(w)])
DDω(w) = (z) -> sum([(i + 1) * i .* w[i] .* z.^(i-1) for i in eachindex(w)])
Φ(inp, k) = begin
    ww = inp.w; xx = inp.x; tt = inp.t;
    ω_kt = ω(ww)(k .* tt.^(-1 / (length(ww) + 1)))
    return 1im*k .- ω_kt .* tt ./ xx 
end
DDΦ(inp, k) = begin
    ww = inp.w; xx = inp.x; tt = inp.t;
    DDω_kt = DDω(ww)(k .* tt.^(-1 / (length(ww) + 1))) .* tt.^(-2 / (length(ww) + 1))
    return - DDω_kt .* tt ./ xx 
end
P(inp, k) = begin
    xx = inp.x; 
    result = exp(xx * Φ(inp, k))
    result ≈ 0.0 ? 0 : result
end
function MinN(X, Y,N)
    out = []
    for i1 in eachindex(X)
        x = X[i1]
        n = N[i1]
        p = []
        temp = abs.(Y .- x)
        for i2 = 1:n
            p = push!(p,argmin(temp))
            temp[p[end]] = Inf
        end
        push!(out, Y[p])
    end
    return out
end
function edgeθ(w)
    n = length(w)+1
    wend = w[end];
    A = real(wend); B = imag(wend);
    θ = atan(A/B)/n
    θs = [θ]
    i1 = 1
    while 0 <= θ - i1*π/n
        push!(θs,θ - i1*π/n)
        i1 += 1
    end
    i1 = 1
    while θ + i1*π/n <= 2*π
        append!(θs,θ + i1*π/n)
        i1 += 1
    end
    if (θs[1]≈ θs[end]-2*π)||(θs[end] ≈ θs[1]+2*π)
        return [θs[end-2]-2*π; θs[end-1]-2*π;θs;θs[2]+2*π;θs[3]+2*π]

    else    
        return [θs[end-2]-2*π; θs[end-1]-2*π;θs[end]-2*π;θs;θs[1]+2*π;θs[2]+2*π;θs[3]+2*π]
    end
end
function sortreal(xx)
    return sort(xx, by=real)
end
function OrderRoots(xx)
    angles = angle2pi(xx)
    sorted_indices = sortperm(angles)
    return xx[sorted_indices]
end
function kk0(inp)
    ww = inp.w; xx = inp.x; tt = inp.t; xxθ = inp.xθ

    m = length(ww);
    a = zeros(m) .* 1im
    a[1] = 1im*xx/(-(m+1)*ww[m])
    for i1 = 1:(m-1)
        a[i1+1] = (i1+1)*ww[i1]*tt^((m-i1)/(m+1))/ ((m+1)*ww[m] )
    end
    A = zeros(m,m) * 1im;
    for i1=1:m-1
        A[i1+1,i1] = 1.0;
    end
    A[:,m] = -1. .* a
    E = A |> eigen
    result = (E.values) #|> filter(k -> real(ww[end]*(k)^(m+1))<=0.01)
    return OrderRoots(result)
end
function find_indices(values, target)
    indices = Int[]
    for val in values
        push!(indices, findfirst(isequal(val), target))
    end
    return indices
end
function boundaryθ(inp)
    w = inp.w;
    n = length(w)+1;
    θs = edgeθ(w);
    θs_bis = (θs[2:end] + θs[1:end-1])/2;

    start_ind = findfirst(x -> x > 0, θs)-1
    stop_ind = findlast(x -> x <= 2*pi, θs)

    off_set = 1
    if real(w[end]*exp(1im*n*θs_bis[start_ind])) <= 0
        off_set = 0
    end

    Nneigh = 0*θs_bis[start_ind:(stop_ind-1)] .+ 2
    θneigh = sortreal.(MinN(θs_bis[start_ind:(stop_ind-1)],θs,Nneigh));
    Nneigh_bis = 0*θs_bis[start_ind:(stop_ind-1)] .+ 3
    θneigh_bis = sortreal.(MinN(θs_bis[start_ind:(stop_ind-1)],θs_bis,Nneigh_bis));

    bd = [[θneigh_bis[i1][1],θneigh[i1][1],θneigh[i1][2],θneigh_bis[i1][3]] for i1 = (1+off_set):2:length(θneigh)]
    return bd
end
function containing_bds(inp,z)
    zr, zθ = polar(z);
    w = inp.w;
    n = length(w) + 1
    if real(w[end]*exp(1im*n*zθ)) < -1e-14
        @warn "start or endpoint is in the danger zone."
    end
    θs = edgeθ(w);
    N = abs.(real(w[end]*z^(n))) < 1e-15 ? 3 : 2
    θ = sortreal.(MinN(zθ,θs,N));
    θ = real(w[end]*exp(1im*n*(θ[1][1]+zθ)/2)) >= 0 ? θ[1][1:2] : θ[1][2:3]
    θs[find_indices(θ, θs)] 
end
function next_to_indices(inp,z,BD)
    z_bd = containing_bds(inp,z)
    ind_bd = z_bd*0;
    for i1 in eachindex(BD)
        bd = BD[i1]
        if z_bd[1] in bd
            ind_bd[1] = i1
        elseif z_bd[2] in bd
            ind_bd[2] = i1
        end
    end
    if any(ind_bd .≈ 0.)
        num = sum(ind_bd)
        ind_bd = [1,length(BD)]
    end
    return Int.(ind_bd |> filter(k -> k != 0))
end
function CP_assigner(K0,bd)
    CP_bd = complex.(zeros(length(bd)))
    CP_in = Vector{Bool}(undef, length(bd))
    CP_in[:] .= false
    for i1 in eachindex(K0)
        CP = K0[i1]
        CPr, CPθ = polar(CP)
        for i2 in eachindex(bd)
            if bd[i2][2]<=CPθ && CPθ <=bd[i2][3]
                CP_bd[i2] = CP
                CP_in[i2] = true
                break
            end
            if (bd[i2][1]<=CPθ && CPθ <bd[i2][2]) || (bd[i2][3]<CPθ && CPθ <bd[i2][4])
                CP_bd[i2] = CP
                CP_in[i2] = false
                break
            end
        end
    end
    if sign(bd[1][1]) != sign(bd[1][4])
        CP = K0[end]
        CPr, CPθ = polar(CP)
        if bd[1][2]+2*pi<=CPθ
            CP_bd[1] = CP
            CP_in[1] = true
        end
        if (bd[1][1] + 2*pi <=CPθ && CPθ <bd[1][2] + 2*pi)
            CP_bd[1] = CP
            CP_in[1] = false
        end
    end
    CP_bd_filt = (CP_bd) |> filter(x -> x == 0 )
    if length(CP_bd_filt) > 1
        @warn "Not all critical points are assigned"
    end
    if !(0.0 in CP_bd)
        @warn "A critical point has been double assigned"
    end
    return CP_bd, CP_in
end
function get_direction(CP_in,CP_bd,start_bd_ind,stop_bd_ind)
    dir = 1 # Counter clockwise
    # This can be done better
    if length(CP_in) == 2
        if CP_bd[1]!= 0.0+0im
            return 0, 1, 1
        else
            return 0, 2, 2
        end
    end

    ind = start_bd_ind[1]
    if mod(ind+(dir-1),length(CP_in))+1 == start_bd_ind[2]
        ind = start_bd_ind[2]
    end

    if ind in stop_bd_ind && CP_in[ind]
        return dir, starting, ind 
    end

    while ~(ind in stop_bd_ind) && CP_in[ind]
        ind = mod(ind+(dir-1),length(CP_in))+1
        if  ind in stop_bd_ind
            return dir, starting, ind
        end 
    end

    dir = -1
    # This can be done better
    ind = start_bd_ind[1]
    ind = start_bd_ind[1]
    if mod(ind+(dir-1),length(CP_in))+1 == start_bd_ind[2]
        ind = start_bd_ind[2]
    end
    starting = ind

    if ind in stop_bd_ind && CP_in[ind]
        return dir, starting, ind 
    end

    while ~(ind in stop_bd_ind) && CP_in[ind]
        ind = mod(ind+(dir-1),length(CP_in))+1
        if  ind in stop_bd_ind
            return dir, starting, ind 
        end 
    end
    @warn "No path found, we need to go through an area with no CP in it."
    dir = 1
    # This can be done better
    ind = start_bd_ind[1]
    if mod(ind+(dir-1),length(CP_in))+1 == start_bd_ind[2]
        ind = start_bd_ind[2]
    end
    starting = ind
    if ind in stop_bd_ind && CP_bd[ind] != 0.0
        return dir, starting, ind 
    end

    while ~(ind in stop_bd_ind) && CP_bd[ind] != 0.0
        ind = mod(ind+(dir-1),length(CP_in))+1
        if  ind in stop_bd_ind
            return dir, starting, ind
        end 
    end

    dir = -1
    # This can be done better
    ind = start_bd_ind[1]
    if mod(ind+(dir-1),length(CP_in))+1 == start_bd_ind[2]
        ind = start_bd_ind[2]
    end
    starting = ind
    if ind in stop_bd_ind && CP_bd[ind] != 0.0
        return dir, starting, ind 
    end

    while ~(ind in stop_bd_ind) && CP_bd[ind] != 0.0
        ind = mod(ind+(dir-1),length(CP_in))+1
        if  ind in stop_bd_ind
            return dir, starting, ind 
        end 
    end
    @warn "No path found, start or end is inbetween inescapable danger zones"
    return 0
end
function GlobalR(inp)
    w = inp.w; x = inp.x; t = inp.t;
    return (30 + abs(x)^2 + 1000 / length(w))^(1 / (length(w) + 1))
end
function mod_offset(x, m, offset)
    return offset .+ mod.(x .- offset, m)
end
function get_CP_ent_ext(inp::Input_sf,bd::Vector,K0::Vector,i1::Number)
    x = inp.x; t = inp.t; w = inp.w;
    αs = angle.(map(k-> DDΦ(inp,k),K0)); 
    θ_ent_ext = mod_offset.(-αs ./ 2 .+ π/2, π, -π / 2);
    maxR=[]
    if length(w) > 1
        for i in 1:length(bd)
            θ = bd[i][1]
            ϕ = θ - θ_ent_ext[i]
            γ = K0[i] * exp(-1im * θ_ent_ext[i])
            s = imag(γ) / imag(exp(1im * ϕ))
            r1 = s * exp(1im * ϕ) - γ
            θ = bd[i][4]
            ϕ = θ - θ_ent_ext[i]
            γ = K0[i] * exp(-1im * θ_ent_ext[i])
            s = imag(γ) / imag(exp(1im * ϕ))
            r2 = s * exp(1im * ϕ) - γ
            push!(maxR, minimum(abs.([r1, r2])))
        end
    else
        maxR = [Inf,Inf]
    end
    scaleR = 20*sqrt.(abs.(1 ./ (αs .* x)))

    R = [minimum([scaleR[i1];maxR[i1]]) for i1 in eachindex(scaleR)];
    
    CP_ent = K0 + R.*exp.(1im .* θ_ent_ext)
    CP_ext = K0 - R.*exp.(1im .* θ_ent_ext)
    return CP_ent, CP_ext
end
function get_CP_ent_ext(inp::Input_sf,bd::Vector,K0::Vector)
    w = inp.w; x = inp.x; t = inp.t;
    θs = edgeθ(w)
    θent  = zeros(length(bd))
    θext  = zeros(length(bd))
    for i1 = 1:length(bd)
        θext[i1] = bd[i1][1]
        θent[i1] = bd[i1][4]
    end
    Rinf = GlobalR(inp)

    αs = angle.(map(k-> DDΦ(inp,k),K0))
    θmid = -αs/2

    R,θ = polar(K0)
    Δθent = θent-θmid;
    Δθext = θmid-θext;
    CPent = cos.(θ-θmid).*R./cos.(Δθent).*exp.(1im*θent)
    CPext = cos.(θ-θmid).*R./cos.(Δθext).*exp.(1im*θext)
    CP_ent = [];
    CP_ext = [];
    for i1 in eachindex(CPent)
        push!(CP_ent, [CPent[i1] + Rinf*exp(1im*θent[i1]),CPent[i1] + K0[i1]] ./ 2)
        push!(CP_ext, [CPext[i1] + K0[i1],CPext[i1] + Rinf*exp(1im*θext[i1])] ./ 2)
    end
    return CP_ent, CP_ext
end
function get_CP_ent_ext_small(inp::Input_sf,bd::Vector)
    w = inp.w; x = inp.x; t = inp.t;
    θs = edgeθ(w)
    θent  = zeros(length(bd))
    θext  = zeros(length(bd))
    for i1 = 1:length(bd)
        θext[i1] = bd[i1][1]
        θent[i1] = bd[i1][4]
    end
    if length(w) == 1;
        Rsmall = 2.5
        CP_ent = .5*1im .+ Rsmall*exp.(1im*θent)
        CP_ext = .5*1im .+ Rsmall*exp.(1im*θext)
    else
        Rsmall = .5
        CP_ent = Rsmall*exp.(1im*θent)
        CP_ext = Rsmall*exp.(1im*θext)
    end
    return CP_ent, CP_ext
end
function Danger_zone_maker(inp::Input_sf)
    w = inp.w; x = inp.x; xθ = inp.xθ; t = inp.t; start = inp.start; stop = inp.stop;
    K0 = kk0(inp)
    bd = boundaryθ(inp)
    start_bd_ind = next_to_indices(inp,start,bd) 
    stop_bd_ind = next_to_indices(inp,stop,bd)
    if ((start_bd_ind[1] == stop_bd_ind[1] && start_bd_ind[2] == stop_bd_ind[2]) || (start_bd_ind[1] == stop_bd_ind[2] && start_bd_ind[2] == stop_bd_ind[1])) && length(w) !=1
        @warn "Stop and startpoint are in the same analytic domain, do simple integration over a connecting line."
    end
    Danger_zone_m = []

    # CP_in
    (CP_bd, CP_in) = CP_assigner(K0,bd)
    #dir

    dir, start_ind, stop_ind = get_direction(CP_in,CP_bd,start_bd_ind,stop_bd_ind)
    if length(w) == 1
        dir = -Int(sign(real(stop)-real(start)));
        if dir == 0
            dir = -Int(sign(imag(start)));
        end
    end
    (CP_ent, CP_ext) = get_CP_ent_ext(inp,bd,CP_bd,2)#Using bisection to get two intermediate steps
    #(CP_ent, CP_ext) = get_CP_ent_ext(inp,bd,CP_bd) #Using variance of Gaussian to estimate distance

    ## Defor_points
    bdinf = []
    Rinf = GlobalR(inp);
    if length(w) == 1
        push!(bdinf, Rinf*exp.(1im*bd[1][[4,1]]) .+ CP_bd[1])
        push!(bdinf, Rinf*exp.(1im*bd[2][[4,1]]) .+ CP_bd[2])
    else
        for i1 in eachindex(bd)
            push!(bdinf, Rinf*exp.(1im*bd[i1][[4,1]]))
        end
    end
    CP_ent_ext = [[MinN([bdinf[i1][1]],[CP_ent[i1];CP_ext[i1]],1)[1];MinN([bdinf[i1][1]],[CP_ent[i1];CP_ext[i1]],2)[1][2]] for i1 in eachindex(bdinf)]
    Defor_points = [[bdinf[i1][1];CP_ent_ext[i1][1];CP_bd[i1];CP_ent_ext[i1][2];bdinf[i1][2]] for i1 in eachindex(bdinf)]
    Defor_points = dir == 1 ? reverse.(Defor_points) : Defor_points;
    Defor_zones = []
    for i1 in eachindex(Defor_points)
        push!(Defor_zones, Defor_points[i1])
    end
    return Defor_zones, dir, start_ind, stop_ind
end
function Danger_zone_small(inp::Input_sf)
    w = inp.w; x = inp.x; xθ = inp.xθ; t = inp.t; start = inp.start; stop = inp.stop;
    K0 = kk0(inp)
    bd = boundaryθ(inp)
    start_bd_ind = next_to_indices(inp,start,bd) 
    stop_bd_ind = next_to_indices(inp,stop,bd)
    if ((start_bd_ind[1] == stop_bd_ind[1] && start_bd_ind[2] == stop_bd_ind[2]) || (start_bd_ind[1] == stop_bd_ind[2] && start_bd_ind[2] == stop_bd_ind[1])) && length(w) !=1
        @warn "Stop and startpoint are in the same analytic domain, do simple integration over a connecting line."
    end

    # CP_in
    (CP_bd, CP_in) = CP_assigner(K0,bd)
    #dir

    dir, start_ind, stop_ind = get_direction(CP_in,CP_bd,start_bd_ind,stop_bd_ind)
    if length(w) == 1
        dir = -Int(sign(real(stop)-real(start)));
        if dir == 0
            dir = -Int(sign(imag(start)));
        end
    end
    (CP_ent, CP_ext) = get_CP_ent_ext_small(inp,bd)#Using bisection to get two intermediate steps

    ## Defor_points
    bdinf = []
    Rinf = GlobalR(inp);
    if length(w) == 1
        push!(bdinf, reverse(Rinf*exp.(1im*bd[1][[1,4]]) .+ CP_bd[1]/abs(CP_bd[1])*.5))
        push!(bdinf, reverse(Rinf*exp.(1im*bd[2][[1,4]]) .+ CP_bd[2]/abs(CP_bd[2])*.5))
    else
        for i1 in eachindex(bd)
            push!(bdinf, reverse(Rinf*exp.(1im*bd[i1][[1,4]])))
        end
    end
    Defor_points = [[bdinf[i1][1];CP_ent[i1];CP_ext[i1];bdinf[i1][2]] for i1 in eachindex(bdinf)]
    Defor_points = dir == 1 ? reverse.(Defor_points) : Defor_points;
    Defor_zones = []
    for i1 in eachindex(Defor_points)
        push!(Defor_zones, Defor_points[i1])
    end
    return Defor_zones, dir, start_ind, stop_ind
end
function FullPath(inp)
    w = inp.w; x = inp.x; t = inp.t; xθ = inp.xθ; wθ = inp.wθ
    if abs(x) < 0.1-eps() || isempty(kk0(inp))
        Defor_zones, dir, start_ind, stop_ind = Danger_zone_small(inp)
        small = true 
    else 
        Defor_zones, dir, start_ind, stop_ind = Danger_zone_maker(inp); 
        small = false
    end
    i1 = start_ind
    connect = []
    category = []
    while true
        Defor_p = Defor_zones[i1]
        NN = length(Defor_p);
        NCP = Int(ceil(NN/2));

        push!(connect,[Defor_p[1],Defor_p[2]])
        push!(category, ["inf","CP_ent"])
        if small
            for i2 = 2:(NN-2)
                push!(connect,[Defor_p[i2],Defor_p[i2+1]])
                push!(category, ["CP_ent","CP_ext"])
            end
        else
            for i2 = 2:(NCP-2)
                push!(connect,[Defor_p[i2],Defor_p[i2+1]])
                push!(category, ["CP_ent","CP_ent"])
            end

            push!(connect,[Defor_p[NCP-1],Defor_p[NCP]])
            push!(category, ["CP_ent","CP"])
            push!(connect,[Defor_p[NCP],Defor_p[NCP+1]])
            push!(category, ["CP","CP_ext"])

            for i2 = (NCP+1):(NN-2)
                push!(connect,[Defor_p[i2],Defor_p[i2+1]])
                push!(category, ["CP_ext","CP_ext"])
            end
        end
        push!(connect,[Defor_p[NN-1],Defor_p[NN]])
        push!(category, ["CP_ext","inf"])
        if i1 == stop_ind
            break
        end
        i1 = Int(mod(i1+(dir-1),length(Defor_zones)) + 1)
    end

    return Deformation(connect,category,w,xθ,dir)
end 
function My_Integrate(int_f,Defor,N)
    res = 0im;
    for i1 = 1:length(Defor.tt)
        s = curv(Defor.path[i1],Defor.tt[i1][1],Defor.tt[i1][end],Defor.Dpath[i1],N)
        if Defor.meth[i1] == "Legendre"
            f = stand_int(int_f,s)
            x, w = gausslegendre(N);
            res += dot(w,f.(x))
        elseif Defor.meth[i1] == "Clenshaw-Curtis"
            res += Clen_Curt(int_f,s)
        end
    end
    return res
end
function Integrand(inp, g)
    DD = FullPath(inp)
    integrand = z -> g(z) * P(inp, z)
    return integrand, DD
end
function Residue(g,f,z)
    F = z -> g.(z).*f.(z)
    s = curv( t -> z .+ exp.(1im*t), 0, 2*π, t -> 1im*exp.(1im*t),500)
    Clen_Curt(F,s)/(2*π*1im)
end

function SpecialFunction(w::Vector, x::Vector, t::Vector, start::Number, stop::Number, N::Number, gg::Function, Dgg::Vector, poles::Vector, m::Vector)
    n = length(w) + 1
    if real(start^n*w[end]) < 0
        @warn "Deformation cannot be made: exponential growth at start"
    end
    if real(stop^n*w[end]) < 0
        @warn "Deformation cannot be made: exponential growth at stop"
    end
 
    base_inp = Input_sf(w,x[1] * t[1]^(-1/n)+ eps(),t[1],start,stop);
    DDD, init_dir, start_init_ind, stop_init_ind = Danger_zone_maker(base_inp)
    Res = Complex.(zeros(length(x),length(t)))
    #Rotating
    elapsed_time_path = 0;
    elapsed_time_int = 0;
    for i1 = 1:length(x)
        xθ = angle(x[i1])
        w_i1 = [ w[j] * (exp(-1im*xθ))^(j + 1) for j in 1:length(w)];
        start_i1 = start * exp(1im*xθ);
        stop_i1 = stop * exp(1im*xθ);
    #Get integrand
        for i2 = 1:length(t)
            g = z -> gg(z * exp(-1im*xθ) * t[i2]^(-1/n));
            inp = Input_sf(w_i1,x[i1] * t[i2]^(-1/n)+ eps(),t[i2],start_i1,stop_i1);
            t1 = time();

            integrand, DD = Integrand(inp, g);
            elapsed_time_path += time() - t1
            #DomainPlot(DD) |> display
            #PathPlot(DD,integrand)
            t1 = time();

            vals = My_Integrate(integrand, DD,N)

            for i3 = 1:length(m)
                if m[i3] != -1
                    vals -=  (init_dir * DD.dir < 0 && m[i3] < 0) ? 2 * π * 1im * Residue(z -> g(z), z-> P(inp, z), poles[i3]) : 0
                else 
                    vals -=  (init_dir * DD.dir < 0) ? 2 * π * 1im ./ (Dgg[i3](poles[i3] * exp(-1im*xθ) * t[i2]^(-1/n))*exp(-1im*xθ) * t[i2]^(-1/n)) : 0
                end
            end
            elapsed_time_int += time() - t1
            Res[i1,i2] = (t[i2])^(-1/n) * exp(-1im*xθ) * vals
        end
    end
    #println("Elapsed pathing time: ", elapsed_time_path, " seconds");
    #println("Elapsed integrating time: ", elapsed_time_int, " seconds");
    return Res
end
function SpecialFunction(w::Vector, x::Number, t::Number, start::Number, stop::Number, N::Number, gg::Function, Dgg::Vector, poles::Vector, m::Vector)
    R = SpecialFunction(w, [x], [t],start, stop, N, gg, Dgg, poles, m)
    R[1]
end
function SpecialFunction(w::Vector, x::Vector, t::Number, start::Number, stop::Number, N::Number, gg::Function, Dgg::Vector, poles::Vector, m::Vector)
    R = SpecialFunction(w, x, [t],start, stop, N, gg, Dgg, poles, m)
    R[:,1]
end
function SpecialFunction(w::Vector, xx::Vector, tt::Vector, m::Number, N::Number)
    g = z -> (1im*z).^m
    Dg = z -> (1. * 1im).^(-m)
    SpecialFunction(w, xx, tt,-1,1, N, g, [Dg], [0], [m])
end
function SpecialFunction(w::Vector, xx::Number, tt::Number, m::Number, N::Number)
    R = SpecialFunction(w, [xx], [tt], m, N)
    R[1]
end
function SpecialFunction(w::Vector, xx::Vector, tt::Number, m::Number, N::Number)
    R = SpecialFunction(w, xx, [tt], m, N)
    R[:,1]
end 
function SpecialFunction(w::Vector, xx::Number, tt::Number, N::Number, gg::Function, poles::Vector, m::Vector)
    R = SpecialFunction(w, [xx], [tt], N, gg, poles, m)
    R[1]
end
function SpecialFunction(w::Vector, xx::Vector, tt::Number, N::Number, gg::Function, poles::Vector, m::Vector)
    R = SpecialFunction(w, xx, [tt], N, gg, poles, m)
    R[:,1]
end 

function Path_gif(w::Vector, x::Vector, t::Vector, start::Number, stop::Number, gg::Function, filename::String)
    anim = Animation()
    n = length(w)+1
    for i1 = 1:length(x)
        xθ = angle(x[i1])
        w_i1 = [ w[j] * (exp(-1im*xθ))^(j + 1) for j in 1:length(w)];
        start_i1 = start * exp(1im*xθ);
        stop_i1 = stop * exp(1im*xθ);
        for i2 = 1:length(t)
            g = z -> gg(z * exp(-1im*xθ) * t[i2]^(-1/n));
            inp = Input_sf(w_i1,x[i1] * t[i2]^(-1/n)+ eps(),t[i2],start_i1,stop_i1);
            integrand, DD = Integrand(inp, g);
            p1 = DomainPlot(DD,x[i1],t[i2]);
            p2 = PathPlot(DD,integrand);
            l = @layout [grid(1,2)]
            plot(p1,p2; layout = l);
            frame(anim)
        end
    end
    gif(anim, filename)
end

function Path_gif(w::Vector, x::Vector, t::Number, start::Number, stop::Number, gg::Function, filename::String)
    Path_gif(w, x, [t], start, stop, gg,filename)
end
