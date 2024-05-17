module Master_func

using Plots, LinearAlgebra, FastGaussQuadrature, SparseArrays, Printf, SpecialFunctions
import Plots: plot, +, -
import Base: +, -, /, *

export UltraFun, Diff, plot
export DomainPlot, PathPlot, SpecialFunction
export Clen_Curt_no_adapt, Clen_Curt, m_Filon_Clen_Curt, Trap, Laguerre_quad, two_sided_Laguerre_quad, Hermite_quad
export curv, sp, standardChop
export Heat_eq

include.(("Myquad.jl","Cheb.jl","Res_pres.jl","Misc.jl","Spf.jl","First_NUTM.jl"))

end