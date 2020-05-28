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
first_mast = 100.0
mast_freq = 100.0
mast_length = 10.0
t_end = mast_freq * 5
t_span = (0.0, t_end)

mast_starts = first_mast:mast_freq:t_end
mast_ends = mast_starts .+ mast_length
mast_event_times = sort(union(mast_starts, mast_ends))

masting_event(u, t, integrator) = t ∈ mast_event_times
function forcing!(integrator)
    if integrator.t ∈ mast_starts
        integrator.p.K = 1.5 * integrator.p.K_base
    elseif integrator.t ∈ mast_ends
        integrator.p.K = integrator.p.K_base
    end
    return
end
cb = DiscreteCallback(masting_event, forcing!)

let
    u0 = [1.0, 0.5, 0.1]
    t_grid = range(0.0, t_end, length = 10000)

    # Chain
    par_chain = ModelPar(a_CP = 0.3, ω = 0.0, K = 3.0)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain))
    sol_chain = solve(prob_chain, reltol = 1e-8, abstol = 1e-8)
    sol_chain_grid = sol_chain(t_grid)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = mast_event_times)
    sol_chain_mast = solve(prob_chain, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_chain_mast_grid = sol_chain_mast(t_grid)

    # Omnivory
    par_omn = ModelPar(a_CP = 0.3, ω = 0.05, K = 3.0, pref = adapt_pref)

    prob_omn = ODEProblem(model!, u0, t_span, deepcopy(par_omn), tstops = mast_event_times)
    sol_omn = solve(prob_omn, reltol = 1e-8, abstol = 1e-8)
    sol_omn_grid = sol_omn(t_grid)

    prob_omn = ODEProblem(model!, u0, t_span, deepcopy(par_omn), tstops = mast_event_times)
    sol_omn_mast = solve(prob_omn, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_omn_mast_grid = sol_omn_mast(t_grid)

    # Layout
    fig = figure(figsize = (10, 8))
    R_col = "#1f77b4"
    C_col = "#ff7f0e"
    P_col = "#2ca02c"
    subplot(3, 1, 1)
    plot(sol_chain_grid.t, sol_chain_grid[1, :], color = R_col, alpha = 0.5)
    plot(sol_chain_grid.t, sol_chain_grid[2, :], color = C_col, alpha = 0.5)
    plot(sol_chain_grid.t, sol_chain_grid[3, :], color = P_col, alpha = 0.5)
    plot(sol_chain_mast_grid.t, sol_chain_mast_grid[1, :], color = R_col, label = "R")
    plot(sol_chain_mast_grid.t, sol_chain_mast_grid[2, :], color = C_col, label = "C")
    plot(sol_chain_mast_grid.t, sol_chain_mast_grid[3, :], color = P_col, label = "P")
    legend()
    title("Chain")

    subplot(3, 1, 2)
    plot(sol_omn_grid.t, sol_omn_grid[1, :], color = R_col, alpha = 0.5)
    plot(sol_omn_grid.t, sol_omn_grid[2, :], color = C_col, alpha = 0.5)
    plot(sol_omn_grid.t, sol_omn_grid[3, :], color = P_col, alpha = 0.5)
    plot(sol_omn_mast_grid.t, sol_omn_mast_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_mast_grid.t, sol_omn_mast_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_mast_grid.t, sol_omn_mast_grid[3, :], color = P_col, label = "P")
    legend()
    title("Omnivory")

    subplot(3, 1, 3)
    plot(sol_omn_grid.t,  [adapt_pref(u, par_omn, 0.0) for u in sol_omn_grid], color = "#000000")
    plot(sol_omn_mast_grid.t,  [adapt_pref(u, par_omn, 0.0) for u in sol_omn_mast_grid], color = "#13b8dd")
    title("Preference")

    tight_layout()
    return fig
end
