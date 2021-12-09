# install package 
printstyled("1. Installing packages\n", color = :blue)
include("install_packages.jl")

# Code to reproduce Figure 3 

## Figure 3 
printstyled("2. Runing pulse analysis and creating Fig. 3\n", color = :blue)
include("fig_pulse.jl")

## Figure 4
printstyled("3. Runing press analysis and creating Fig. 4\n", color = :blue)
include("fig_press.jl")

## NB: Part A and Part B were combined together and some additions to them 
## have been made to them (e.g. background color) using Inkscape. 


# Figure S1
# Simply use the other set of paramters

# Figure S2 takes ~2hours
# printstyled("Runing simulations for FigS1\n", color = :blue)
# include("fig_S1.jl")

# Figue S3 takes ~2hours
# printstyled("Runing simulations for FigS1\n", color = :blue)
# include("fig_S2.jl")