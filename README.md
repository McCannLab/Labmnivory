# Labmnivory
[![Reproduce the analysis](https://github.com/McCannLab/Labmnivory/actions/workflows/reproduce.yaml/badge.svg)](https://github.com/McCannLab/Labmnivory/actions/workflows/reproduce.yaml)


## Installation


Simulation are implemented in [Julia](https://julialang.org/), a recent version of Julia is thus  required to reproduce this analysis (note that the 
code is only tested for Julia v1.6.0) along with the following packages (note that they can be installed using `src/install_packages.jl` in this repository, see below):

|Package              | Links                                                  |
|:--------------------|:-------------------------------------------------------|
|DifferentialEquations| [repo](https://github.com/SciML/DifferentialEquations.jl)|
|ForwardDiff          | [repo](https://github.com/JuliaDiff/ForwardDiff.jl)    |
|LinearAlgebra        | [doc](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)|
|NLsolve              | [repo](https://github.com/JuliaNLSolvers/NLsolve.jl)   |
|Parameters           | [repo](https://github.com/mauro3/Parameters.jl)        |
|PyPlot               | [repo](https://github.com/JuliaPy/PyPlot.jl)           |
|QuadGK               | [repo](https://github.com/JuliaMath/QuadGK.jl)         |
|RecursiveArrayTool   | [repo](https://github.com/SciML/RecursiveArrayTools.jl)|
|Statistics           | [repo](https://docs.julialang.org/en/v1/stdlib/Statistics/)|


Then, clone this repository using [git](https://git-scm.com/):

```sh
git clone https://github.com/McCannLab/Labmnivory.git
```

or downloaded the [zip file](https://github.com/McCannLab/Labmnivory/archive/refs/heads/master.zip). If you haven't yet installed the packages, you can actually 
use `src/install_packages.jl` in this repository. To do so, in a terminal, set your working directory at the root of the freshly cloned/downloaded repository, then use the following command line&nbsp;:

```julia
julia src/install_packages.jl
```

You should now be ready to reproduce the analysis.


## Repo structure 

- `src`: source code 
- `figs`: folder where figures are exported



## How to reproduce the analysis 

In a terminal, assuming the working directory is set to be the the root of the freshly cloned/downloaded repository, use the following command line&nbsp;:

```julia
julia src/pipeline.jl
```

