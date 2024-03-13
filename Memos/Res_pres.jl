using  Plots, LinearAlgebra
mutable struct curv
    c # Curve
    a # Start value
    b # End value
    w # Weights for integration
end 

function stand_int(f,s)
    
    f_curv = x-> map(f,map(s.c,x))  .* map(s.w,x); # map to real line
    
    Dtrans_func = (s.b-s.a)/2;
    trans_func = x -> (s.b + s.a)/2 + x * (s.b - s.a)/2 ;
    
    g = x -> f_curv(trans_func(x)) * Dtrans_func; # map to [-1,1]
    return g
end

function grid_weight_find(a,b)
    
    J = Matrix(SymTridiagonal(a,b[1:(end-1)]))
    
    xgrid, U = eigen(J);
    
    U = transpose(sign.(U[1,:])) .* U
    
    w = abs.(U[1,:]).^2
    
    return xgrid, w, U
    
end

function int_err_plot(meth,f,true_res,s, NN = 1000-1)
    res = zeros(NN)
    
    counter = NN

    for i1 = 1:NN
        res[i1] = meth(f,s,i1+1) - true_res
        if isnan(res[i1])
            counter = i1-1
            break
        end
    end
    
    plot_res = res[1:counter];

    p = plot(layout = grid(1, 2),size=(1000,400));
    plot!(s.a : .1 : s.b, map(f,s.c.(s.a : .1 : s.b)), subplot = 1, label = "", title ="Function")
    plot!((1:length(plot_res )) .+ 1,abs.(plot_res ) .+ 10^(-17),yaxis =:log,subplot = 2, label = "", title ="Integral error", show = true)
    p
end