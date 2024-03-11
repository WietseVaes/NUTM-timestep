include("myquad.jl")
include("myCheb.jl")
include("myft.jl")

using NBInclude
nbexport("Master_func.jl", "Master_func.ipynb")