# Here is the code to install packages. 
# See readme 
using Pkg
Pkg.add(
    ["Parameters", "RecursiveArrayTools", "Statistics", "DifferentialEquations",
    "PyPlot", "ForwardDiff", "LinearAlgebra", "QuadGK", "NLsolve", "FileIO", "JLD2"]
)