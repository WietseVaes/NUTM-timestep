module Master_func

using Plots, LinearAlgebra, FastGaussQuadrature, SparseArrays, Printf, SpecialFunctions, OperatorApproximation, Revise
import Plots: plot, +, -
import Base: +, -, /, *

include.(("Myquad.jl","Cheb.jl","Res_pres.jl","Misc.jl","Spf.jl","First_NUTM.jl","myCheb.jl","Myft.jl"))

export UltraFun, Diff, plot
export DomainPlot, PathPlot, SpecialFunction, Path_gif, Integrand, Input_sf
export Clen_Curt_no_adapt, Clen_Curt, m_Filon_Clen_Curt, Trap, Laguerre_quad, two_sided_Laguerre_quad, Hermite_quad, Levin_quad
export curv, sp, standardChop
export Heat_eq, int_after_IBP_map,int_after_IBP_trunc
export Cheb_ft
export Levin_ft

end