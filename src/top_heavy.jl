using Statistics

"""
...
# Conpute the average over time of four biomass ratios
- `sol`: solution of the ODE problem.
- `burn::Int `: number of time steps to be thrown out (burn-in).
...
"""
function top_heavy(sol::Array{T}, burn::Int = 0)::Array{T} where T
    burn += 1 # could be time rather than number of time steps
    return [
        mean(sol[3, burn:end] ./ (sol[1, burn:end] .+ sol[2, burn:end])),
        mean(sol[3, burn:end] ./ sol[1, burn:end]),
        mean(sol[3, burn:end] ./ sol[2, burn:end]),
        mean(sol[1, burn:end] ./ sol[2, burn:end])
    ]
end
"""

