using SparseArrays, LinearAlgebra

include("myquad.jl")

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

#   Trials
#   ≡≡≡≡≡≡≡≡

#   FT through Fokas paper
#   ========================
# 
#   Possible problems: binomial is to big to compute

function α_cheb(m)
    a = zeros(m+1,1);
    
    if m>0
        a[1] = (-1)^m;

        a[2] = (-1)^(m+1) * m^2;
        n = 3:(m+1)
        for n = 3:(m+1)
            for k = 1:(m-n+2)
                j = k:(n+k-3)
                a[n] += binomial(BigInt(n+k-3), BigInt(k-1)) * prod(m .- j)
            end
            a[n] *= (-1)^(m + n - 1) * 2^(n - 2) * m
        end
    else
        a[1] = (-1)^m;
    end
    return a
    
end

function Fokas_Cheb_ft(f, N, λ)
   
    λ = complex(λ)
    
    c = Ultra_spherical_coeff(f, N, λ0)
    res = λ .* 0 
    
    for m = 0:(N-1)
        res_tmp = λ .* 0
        for i1 = 1:length(λ)
            if λ[i1]!=0
                n = 1:(m+1)
                res_tmp[i1] = sum(α_cheb(m) .* (exp.(1im .* λ[i1])./((1im .* λ[i1]).^n) .+ (-1).^(n .+ m).*exp.(-1im .* λ[i1])./((1im .* λ[i1]).^n)))
            else
                if m != 1
                    res_tmp[i1] = ((-1)^(m+1)-1)/(m^2-1);
                else 
                    res_tmp[i1] = 0;
                end
            end
        end
        res += c[m+1] .* res_tmp;
    end
        
    return res
end

#   FT through three-term-recurrence
#   ==================================
# 
#   Possible problems: Three-term-recurrence unstable Possible solutions:
#   BigFloat (Complex matrix)

function Der_sinc(x,N)
    
    # f^(n) = -int_0^1 t^n cos(tx+n*pi/2)dt
    
    s = curv(x -> x,0,1, x ->  1,100)
    
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

function Cheb_ft(N,k,λ)
    
    #Number of Chebyshev polynomials N
    
    #Frequencies k
    
    # f_{j+1} = (i D_k f_j-a_j f_j - b_{j-1} f_{j-1})/b_j
    
    α = sp(λ - 1/2);
    β = sp(λ - 1/2);

    a = [sp((β^2-α^2) / sp((2*j .+ α .+ β) .* (2*j .+ α .+ β .- 2))) for j in 1:N];
    b = [sp(2*sqrt.(j .* (j .+ α) .* (j .+ β) .* (j .+ α .+ β)) / ((2*j .+ α .+ β) .* sqrt.((2*j .+ α .+ β).^2 .- 1))) for j in 1:N];
    
    if λ == 0
        a = 0 .* a
        b = 0 .* b .+ 0.5
    end
    
    Dsinc = map(sp,complex(2 .* Der_sinc(k,N)))
    nx = size(Dsinc,1); ny = size(Dsinc,2)
    
    T = map(sp,complex(zeros(length(k),N+1,N)));
    T[:,:,1] = Dsinc;

    if λ == 0
        T = Add_Cheb_Der(T,2,0,0,sp(1/sqrt(2)));
        T = Add_Cheb_Der(T,3,0,sp(1/sqrt(2)),1/2);
        for i1 = 4:N
            T = Add_Cheb_Der(T,i1,0,sp(1/2),sp(1/2));
        end
    else
        T = Add_Cheb_Der(T,2,a[1],0,b[1]);
        for i1 = 3:N
            T = Add_Cheb_Der(T,i1,a[i1-1],b[i1-2],b[i1-1]);
        end
    end
    return T
end

#   Failed
#   ≡≡≡≡≡≡≡≡
# 
#   Fdersinc: unstable
#   ––––––––––––––––––––

function F_der_sinc(x,N)
    
    # g(x) = xf(x) = sin(x)
    # f^(n) = (g^(n)-\sum_{j=0}^{n-1}f^(j))/x
    
    n = length(x);
    Dsinc = zeros(n,N+1);
    
    Dsinc[:,1] = [k == 0 ? 1 : sin.(k) ./ k for k in x]
    
    sum_prev = Dsinc[:,1];
    
    for i1 = 1:N
            
        Dsin = sin.(pi*i1/2 .+ x)
        
        atzero = mod(i1,2) ==0 ? (-1)^(i1/2)/(i1+1) : 0
        
        Dsinc[:,i1+1] = [x[i2] == 0 ? atzero : (Dsin[i2] - i1*Dsinc[i2,i1]) ./ x[i2] for i2 in 1:length(x)]
    end
    return Dsinc
end
