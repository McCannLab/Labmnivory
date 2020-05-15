"""
...
# Conpute the average over time of four biomass ratios
- `sol`: solution of the ODE problem.
- `burn::Int64 `: number of time steps to be thrown out (burn-in).
...
"""
function top_heavy(sol, burn::Int = 0)
    burn += 1 # could be time rather than number of time steps
    hcat(
        mean(sol[3, burn:end] ./ (sol[1, burn:end] .+ sol[2, burn:end])),
        mean(sol[3, burn:end] ./ sol[1, burn:end]),
        mean(sol[3, burn:end] ./ sol[2, burn:end]),
        mean(sol[1, burn:end] ./ sol[2, burn:end])
    )
end
"""


