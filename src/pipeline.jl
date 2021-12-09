# Install packagess
printstyled("1. Installing packages\n", color = :blue)
include("install_packages.jl")

## Pulse main 
printstyled("2. Runing pulse analysis\n", color = :blue)
include("fig_pulse.jl")
printstyled("Done! fig_pulse.jl created\n\n", color = :green)

## Press main 
printstyled("3. Runing press analysis\n", color = :blue)
include("fig_press.jl")
printstyled("Done! fig_pulse.jl created\n", color = :green)

## Pulse main 
printstyled("4. Runing pulse analysis for SI\n", color = :blue)
include("fig_pulse.jl", ARGS = ["true"])
printstyled("Done! fig_pulse_si.jl created\n\n", color = :green)

## Press main 
printstyled("5. Runing press analysis for SI\n", color = :blue)
include("fig_press.jl", ARGS = ["true"])
printstyled("Done! fig_press_si.jl created\n\n", color = :green)

# Figures S2 and S3 takes ~2hours each and basically call the same functions 
# as previous analyses so we skipped them on GitHub Actions
printstyled("Fig S1 and S2 are skipped\n\n", color = :red)

# Figure S2 takes ~2hours
# printstyled("Runing simulations for FigS1\n", color = :blue)
# include("fig_S1.jl")

# Figue S3 takes ~2hours
# printstyled("Runing simulations for FigS1\n", color = :blue)
# include("fig_S2.jl")