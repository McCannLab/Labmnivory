# Here is the code to install packages. 
# The list below is an example, none of the packages are actually required.
using Pkg
Pkg.add(
    ["Parameters", "RecursiveArrayTools", "Statistics", "DifferentialEquations",
    "PyPlot", "ForwardDiff", "LinearAlgebra", "QuadGK", "NLsolve"]
)