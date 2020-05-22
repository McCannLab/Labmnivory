include("basic_omnivory_module.jl")
using Distributions
using PyPlot

# # Feral Pig Masting Events
# I want to try to model the situation where we have a relatively low
# productivity environment, but we have period high productivity pulses,
# like a oak mastings etc.
#
# Objectives are that we can see a response of omnivory in the face of this
# variable productivity situation, and test that we are measuring omnivory
# in a way that is enlightening.
t_end = 2000.0
t_span = (0.0, t_end)
#NOTE: the step time needs to be odd for the way I am doing this
masting_times = 10:101:t_end
collect(masting_times)
masting_event(u, t, integrator) = t ∈ masting_times
function forcing!(integrator)
    if integrator.t % 2 == 0
        @show "masting!"
        integrator.p.K = 4.0
    else
        @show "regular"
        integrator.p.K = 2.0
    end
    return
end
cb = DiscreteCallback(masting_event, forcing!)

u0 = [1.0, 0.5, 0.1]
par = ModelPar(a_CP = 0.3, ω = 0.1, A = 0.0)
prob = ODEProblem(model!, u0, t_span, par, tstops = masting_times)

sol = solve(prob, reltol = 1e-8, abstol = 1e-8, callback = cb)

let
    t_grid = range(0.0, t_end, length = 10000)
    sol_grid = sol(t_grid)
    plot(sol_grid.t, sol_grid[1, :], label = "R")
    plot(sol_grid.t, sol_grid[2, :], label = "C")
    plot(sol_grid.t, sol_grid[3, :], label = "P")
    legend()
end
