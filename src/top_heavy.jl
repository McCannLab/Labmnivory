using Statistics: mean

"""
...
# Conpute the average over time of four biomass ratios
- `sol`: solution of the ODE problem.
- `burn::Int `: number of time steps to be thrown out (burn-in).
...
"""
function top_heavy(sol, burn::Int = 0)
    burn += 1 # could be time rather than number of time steps
    if burn > size(sol)[2]
        error("parameters: `burn` greater than number of rows in `sol`!")
    end
    return [
        mean(sol[3, burn:end] ./ (sol[1, burn:end] .+ sol[2, burn:end])),
        mean(sol[3, burn:end] ./ sol[1, burn:end]),
        mean(sol[3, burn:end] ./ sol[2, burn:end]),
        mean(sol[1, burn:end] ./ sol[2, burn:end])
    ]
end

# top_heavy(sol, 100)
