using LinearAlgebra

include("../New Code/myquad.jl")
include("../New Code/Misc.jl")

#   Chebyshev discrete Polynomials at x and Coefficients
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function Cheb_poly(N,x,λ)
    
    # Number of cheb polynomials N
    
    # Evaluation points x
    
    # Ultrashperical constant λ
    
    # p_{j+1} = ((x-a_j)p_j - b_{j-1}p_{j-1})/b_j
    # β = α = λ-1/2;
    # d_j = 2*(j+1)+α+β
    # a_{j} = (β^2-α^2)/(d_{j}(d_{j}-2))
    # b_{j} = (2\sqrt{j+1}*\sqrt{(j+1+α)(j+1+β)}*\sqrt{j+1+α+β})/(d_{j}*\sqrt{d_{j}^2-1})
    
    α = λ - 1/2;
    β = λ - 1/2;
    
    j = 1:N;
    dj = 2*j .+ α .+ β;
    a = (β^2-α^2) ./ (dj .* (dj .- 2));
    b = 2*sqrt.(j .* (j .+ α) .* (j .+ β) .* (j .+ α .+ β)) ./ (dj .* sqrt.(dj.^2 .- 1));
    
    T = zeros(length(x),N)
    T[:,1] = ones(length(x),1);
    
    if λ == 0
        T[:,2] = sqrt(2) .* x;
        i1 = 3;
        T[:,i1] = 2 .* x .* T[:,i1-1] - sqrt(2) .* T[:,i1-2]
        for i1 = 4:N
            T[:,i1] = 2 .* x .* T[:,i1-1] - T[:,i1-2]
        end
    else
        T[:,2] = (x .- a[1]) .* T[:,1] ./ b[1];
        for i1 = 2:(N-1)
            T[:,i1+1] = ((x .- a[i1]) .* T[:,i1] - b[i1-1] .* T[:,i1-1])./b[i1]
            
        end
    end
    return T
end

function Gauss_grid_weights(a,b)
    
    J = Matrix(SymTridiagonal(a,b[1:(end-1)]))
    
    xgrid, U = eigen(J);
    
    U = transpose(sign.(U[1,:])) .* U
    
    w = abs.(U[1,:]).^2
    
    return xgrid, w, U
    
end

function Ultra_spherical_coeff(f, N = 100, λ = 0)
    
    # Function f
    
    # Number of Chebychev polynomials N
    
    # Ultrashperical constant λ
    
    α = λ-1/2; 
    β = λ-1/2;

    j = 1:N;
    dj = 2*j .+ α .+ β;
    a = (β^2-α^2) ./ (dj .* (dj .- 2));
    b = 2*sqrt.(j .* (j .+ α) .* (j .+ β) .* (j .+ α .+ β)) ./ (dj .* sqrt.(dj.^2 .- 1));
    
    if (α == -1/2) && (β == -1/2)
        b0 = 1/sqrt(2)
        b[1] = b0;
    elseif (α == 0) && (β == 0)
        a[1] = 0;
    end
    
    xgrid, w, U = Gauss_grid_weights(a,b)

    D = diagm(sqrt.(w));

    c = (U*D)*map(f,xgrid)

    return c
    
end

#   Chebyshev Derivatives
#   ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

function Dk(λ,N)
    
    # Extra constant λ
    
    # Number of Chebychev polynomials N
    
    j = 1:N-1;
    d = (λ) -> j .* sqrt.(2*(λ .+ 1).*(j .+ 2*λ)./(2*j .* λ .+ j))
    D = I;
    for i1 = 0:(λ-1)
        D = spdiagm(1 =>  d(i1)) * D
    end
    return D
end

function Cheb_Der(f, x, k, N = 100)
    
    # Function f
    
    # Evaluation points x
    
    # Derivative k
    
    # Number of Chebychev polynomials N
    
    
    c = Ultra_spherical_coeff(f, N, 0);
    
    D = Dk(k,N);
    
    Pk = Cheb_poly(N,x,k);
    
    uk = Pk*D*c;
end

#   FT through three-term-recurrence
#   ==================================
# 
#   Possible problems: Three-term-recurrence unstable Possible solutions:
#   BigFloat (Complex matrix)

function Der_sinc(x,N,N_quad)
    
    # f^(n) = -int_0^1 t^n cos(tx+n*pi/2)dt
    
    s = curv(x -> x,0,1, x ->  1,N_quad)
    
    n = length(x);
    Dsinc = zeros(n,N+1)
    
    for i1 = 0:N
        
        g = (t,x) -> t.^i1 .* cos.(t .* x .+ i1 * pi/2)
        
        for i2 = 1:n
            Dsinc[i2,i1+1] = Clen_Curt(t -> g(t,x[i2]),s)
        end
    end
    return Dsinc
    
end

function Add_Cheb_Der(D_ft_Cheb,l,ai,bim,bi)
    
    # Dimensions of D_ft_Cheb: -1.    : frequencies
    #                           -2. i1 : level of derivative
    #                           -3.  l : Chebyshev polynomial starting at l=1 (not 0)
    
    
    N = size(D_ft_Cheb,2)
    if l == 2
        for i1 = 1:(N-(l-1))
           D_ft_Cheb[:,i1,l] =  (1im .* D_ft_Cheb[:,i1+1,l-1] - ai .* D_ft_Cheb[:,i1,l-1]) ./ bi
        end
    else 
        for i1 = 1:(N-(l-1))
           D_ft_Cheb[:,i1,l] =  (1im .* D_ft_Cheb[:,i1+1,l-1] - ai .* D_ft_Cheb[:,i1,l-1] - bim .* D_ft_Cheb[:,i1,l-2]) ./ bi
        end
    end
    return D_ft_Cheb
end

function Cheb_ft(N::Number,k::Number,λ::Number,N_quad::Number)
    
    #Number of Chebyshev polynomials N
    
    #Frequencies k
    
    # f_{j+1} = (i D_k f_j-a_j f_j - b_{j-1} f_{j-1})/b_j
    
    α = set_precision(λ - 1/2);
    β = set_precision(λ - 1/2);

    a = [set_precision((β^2-α^2) / set_precision((2*j .+ α .+ β) .* (2*j .+ α .+ β .- 2))) for j in 1:N];
    b = [set_precision(2*sqrt.(j .* (j .+ α) .* (j .+ β) .* (j .+ α .+ β)) / ((2*j .+ α .+ β) .* sqrt.((2*j .+ α .+ β).^2 .- 1))) for j in 1:N];
    
    if λ == 0
        a = 0 .* a
        b = 0 .* b .+ 0.5
    end
    
    Dsinc = map(set_precision,complex(2 .* Der_sinc(k,N,N_quad)))
    nx = size(Dsinc,1); ny = size(Dsinc,2)
    
    T = map(set_precision,complex(zeros(length(k),N+1,N)));
    T[:,:,1] = Dsinc;

    if λ == 0
        T = Add_Cheb_Der(T,2,0,0,set_precision(1/sqrt(2)));
        T = Add_Cheb_Der(T,3,0,set_precision(1/sqrt(2)),1/2);
        for i1 = 4:N
            T = Add_Cheb_Der(T,i1,0,set_precision(1/2),set_precision(1/2));
        end
    else
        T = Add_Cheb_Der(T,2,a[1],0,b[1]);
        for i1 = 3:N
            T = Add_Cheb_Der(T,i1,a[i1-1],b[i1-2],b[i1-1]);
        end
    end
    return T[:,1,:]
end

function Cheb_ft(N::Number,k::Number,λ::Number)
    Cheb_ft(N,k,λ,0)
end
