# Orchestrates scripts

# Install packages
printstyled("1. Installing packages\n", color = :blue)
include("install_packages.jl")


# Main 
## Pulse 
printstyled("2. Runing pulse analysis\n", color = :blue)
include("fig_pulse.jl")
printstyled("Done! fig_pulse.jl created\n\n", color = :green)
## Press 
printstyled("3. Runing press analysis\n", color = :blue)
include("fig_press.jl")
printstyled("Done! fig_pulse.jl created\n", color = :green)


# Supplementary Information 

# The line below is used to run the SI version of the figures: when the 1st
# argument of 'ARGS' is true, then the model parameters used are the ones for 
# the SI figure. Note that we actually use `ARGS` because this is specifically
# the vector that captures external arguments. Thus using it allows us to use 
# arguments in the command line. So here, to run the SI version of fig_pulse, 
# I'd use `julia fig_pulse.jl true`. To mimic this we simply add "true" to 
# ARGS here rather than via an external argument
push!(ARGS, "true")

## Pulse
printstyled("4. Runing pulse analysis for SI\n", color = :blue)
include("fig_pulse.jl")
printstyled("Done! fig_pulse_si.jl created\n\n", color = :green)

## Press 
printstyled("5. Runing press analysis for SI\n", color = :blue)
include("fig_press.jl")
printstyled("Done! fig_press_si.jl created\n\n", color = :green)

# Figures S2 and S3 takes ~2hours each and basically call the same functions 
# as previous analyses so we skipped them on GitHub Actions
printstyled("NB: Fig. S1 and S2 are skipped because they are time/energy-consuming\n\n", color = :red)

# Figure S2 takes ~2hours
# printstyled("Runing simulations for FigS1\n", color = :blue)
# include("fig_S1.jl")

# Figue S3 takes ~2hours
# printstyled("Runing simulations for FigS1\n", color = :blue)
# include("fig_S2.jl")