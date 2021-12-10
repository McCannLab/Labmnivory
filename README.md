# Labmnivory (research compendium)
[![Reproduce the analysis](https://github.com/McCannLab/Labmnivory/actions/workflows/reproduce.yaml/badge.svg)](https://github.com/McCannLab/Labmnivory/actions/workflows/reproduce.yaml)


This repository includes all the scripts to reproduce the theoretical investigations carried out in our paper "On the Dynamic Nature of Omnivory in a Changing World" (a study codenamed 'Labmnivory').


## Installation

Simulations are implemented in [Julia](https://julialang.org/), a recent version of Julia is required to reproduce the analysis (note that the code was only tested for Julia v1.7.0) along with the following packages (note that they can be installed using `src/install_packages.jl` in this repository, see below):

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
|FileIO               | [repo](https://github.com/JuliaIO/FileIO.jl)           |
|JLD2                 | [repo](https://github.com/JuliaIO/JLD2.jl)             |


Then, clone this repository using [git](https://git-scm.com/):

```sh
git clone https://github.com/McCannLab/Labmnivory.git
```

or downloaded the [zip file](https://github.com/McCannLab/Labmnivory/archive/refs/heads/master.zip). If you haven't yet installed the packages, you can actually use `src/install_packages.jl` in this repository. To do so, in a terminal, set your working directory at the root of the freshly cloned/downloaded repository, then use the following command line&nbsp;:

```julia
julia src/install_packages.jl
```

You should now be ready to reproduce the analysis.


## Structure of this repository

## General structure

- `src`: include source code.
- `fig/`: where figures are exported to.
- `fig/final`: final figures (after edition with [Inkscape](https://inkscape.org/)), available as SVG (`/svg/`) and EPS (`/eps/`).


## Scripts (content of `/scr/`)

- `install_packages.jl`: contains code to install packages required;
- `basic_omnivory_module.jl`: includes all basic building blocs including the basic set of ODEs (`model!()`) and the functions to compute the different metrics (e.g. `overshoot()`);
- `pulse.jl`: functions to run pulse analysis;
- `press.jl`: functions to run press analysis;
- `fig_press.jl`: script to create pulse figure panels (Figs. 3 (a-c), 4 (a-c), S3 (a-c), S4 (a-c));
- `fig_pulse.jl`: script to create press figure panels (Figs. 3 (d-f), 4 (d-f), S3 (d-f), S4 (d-f));
- `fig_S1.jl`: script to create figure S1;
- `fig_S2.jl`: script to create figure S2;
- `pipeline.jl`: orchestrates the different scripts to reproduce the analysis.


## How to reproduce the analysis locally

In a terminal, assuming the working directory is the root of the freshly cloned/downloaded repository, use the following command line to reproduce the analysis.

```julia
julia src/pipeline.jl
```



## Questions?

Scripts have been commented for ease of comprehension but if you still have some questions regarding the numerical implementation, do not hesitate to create an issue.

