using Statistics: cov, mean, std
using StatsBase: crosscor

# NB `sol` can be subset (e.g. sol[:, [1 3]]) to obtain the desired pair of species
function asynchrony(sol, wind = 500)
    ws1 = floor(Int, wind / 2)
    ws2 = wind - ws1
    out = zeros(size(sol, 2) - wind)
    for i in (ws1 + 1):(size(sol)[2] - ws2)
        out[i - ws1] = cov(sol[1, (i - ws1):(i + ws2 - 1)], sol[2, (i - ws1):(i + ws2 - 1)])
    end
    return out
end

function asynchrony_cc(sol, wind = 500)
    ws1 = floor(Int, wind / 2)
    ws2 = wind - ws1
    out = zeros(size(sol, 2) - wind)
    for i in (ws1 + 1):(size(sol)[2] - ws2)
        out[i - ws1] = crosscor(sol[1, (i - ws1):(i + ws2 - 1)], sol[2, (i - ws1):(i + ws2 - 1)], [0])[1]
    end
    return out
end

# Coeficient of variation with moving window
function cv(val, wind = 500)
    ws1 = floor(Int, wind / 2)
    ws2 = wind - ws1
    out = zeros(length(val) - wind)
    for i in (ws1 + 1):(length(val) - ws2)
        tmp = val[(i - ws1):(i + ws2 - 1)]
        out[i - ws1] = std(tmp) / mean(tmp)
    end
    return out
end


function auc(val, l_equil, dt, wind = 500)
    ws1 = floor(Int, wind / 2)
    ws2 = wind - ws1
    out = zeros(length(val) - wind)
    for i in (ws1 + 1):(length(val) - ws2)
        out[i - ws1] = dt * sum(abs.(val[(i - ws1):(i + ws2 - 1)] .- l_equil[(i - ws1):(i + ws2 - 1)]))
    end
    return out
end



# Returns the return time defined as the time (here an index of the vector) it
# takes so the sum of lags (abs(f(t) - f(t-1))) for R, C and P is below a
# threshold passed as a parameter.
# NB: sol should be the solution evaluated when the perturbation starts until
# the end of the simulation.
function return_time_sum(sol, dt, threshold = 1e-12)
    tmp = sum(abs.(sol[2:end] .- sol[1:end-1]), dims = 1)
    return dt * findfirst(tmp .< threshold)[2]
end
# use only one component
function return_time(sol, id, dt, threshold = 1e-12)
    tmp = [abs.(sol[id, k] .- sol[id, k+1]) for k in 1:length(sol) - 1]
    return dt * findfirst(tmp .< threshold)
end


# Returns the coeficient of variation (cv) for R, C and P.
# NB: sol should be the solution evaluated when the perturbation starts until
# the end of the simulation.
function global_cv(sol)
    return std(sol, dims = 2) # ./ mean(sol, dims = 2)
end

# Global asynchrony between 2 species
function global_asyn(sol)
    # return crosscor(sol[1, :], sol[2, :], [0])
    return cov(sol[1, :], sol[2, :])
end


# Global area under the curve
function global_auc(sol, sol_equil, dt)
    return dt * sum(abs.(sol - sol_equil), dims = 2)
end