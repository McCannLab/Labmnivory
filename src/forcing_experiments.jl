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
t_end = 55.0 * 10
t_span = (0.0, t_end)
#NOTE: the step time needs to be odd for the way I am doing this
masting_times = (0:(55):t_end)[2:end]
collect(masting_times)
masting_event(u, t, integrator) = t ∈ masting_times
function forcing!(integrator)
    if integrator.t % 2 == 0
        integrator.p.K = 4.0
    else
        integrator.p.K = 3.0
    end
    return
end
cb = DiscreteCallback(masting_event, forcing!)

let
    u0 = [1.0, 0.5, 0.1]
    t_grid = range(0.0, t_end, length = 10000)

    # Chain
    par_chain = ModelPar(a_CP = 0.3, ω = 0.0, A = 0.0, K = 3.0)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain))
    sol_chain = solve(prob_chain, reltol = 1e-8, abstol = 1e-8)
    sol_chain_grid = sol_chain(t_grid)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = masting_times)
    sol_chain_mast = solve(prob_chain, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_chain_mast_grid = sol_chain_mast(t_grid)

    # Omnivory
    par_omn = ModelPar(a_CP = 0.3, ω = 0.05, A = 0.0, K = 3.0)

    prob_omn = ODEProblem(model!, u0, t_span, deepcopy(par_omn), tstops = masting_times)
    sol_omn = solve(prob_omn, reltol = 1e-8, abstol = 1e-8)
    sol_omn_grid = sol_omn(t_grid)

    prob_omn = ODEProblem(model!, u0, t_span, deepcopy(par_omn), tstops = masting_times)
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
    plot(sol_omn_grid.t,  [pref(u, par_omn) for u in sol_omn_grid], color = "#000000")
    plot(sol_omn_mast_grid.t,  [pref(u, par_omn) for u in sol_omn_mast_grid], color = "#ff2233")
    title("Preference")

    tight_layout()
    return fig
end
