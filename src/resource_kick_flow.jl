include("basic_omnivory_module.jl")
using DifferentialEquations
using PyPlot

first_mast = 100.0
mast_freq = 100.0
mast_length = 2.0
t_end = mast_freq * 2
t_span = (0.0, t_end)

mast_starts = first_mast:mast_freq:t_end
mast_ends = mast_starts .+ mast_length
mast_event_times = sort(union(mast_starts, mast_ends))

masting_event(u, t, integrator) = t ∈ mast_event_times
function forcing!(integrator)
    if integrator.t ∈ mast_starts
        integrator.u[1] += 2.0
    end
    return
end
cb = DiscreteCallback(masting_event, forcing!)

let
    u0 = [1.0, 1.5, 1.5]
    t_grid = range(0.0, t_end, length = 10000)
    t_start = 75.0

    # Omnivory
    par_omn = ModelPar(a_CP = 0.25, ω = 0.1, pref = adapt_pref)

    prob_omn = ODEProblem(model!, u0, t_span, deepcopy(par_omn), tstops = mast_event_times)
    sol_omn = solve(prob_omn, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_omn_grid = sol_omn(t_grid)

    fig = figure()
    plot(sol_omn_grid.t, sol_omn_grid.u)

    return fig
end
