include("basic_omnivory_module.jl")
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
mast_length = 2.0
t_end = mast_freq * 2
t_span = (0.0, t_end)

mast_starts = first_mast:mast_freq:t_end
mast_ends = mast_starts .+ mast_length
mast_event_times = sort(union(mast_starts, mast_ends))

mast_strength = 2.0

masting_event(u, t, integrator) = t ∈ mast_event_times
function forcing!(integrator)
    if integrator.t ∈ mast_starts
        integrator.p.K = mast_strength * integrator.p.K_base
    elseif integrator.t ∈ mast_ends
        integrator.p.K = integrator.p.K_base
    end
    return
end
cb = DiscreteCallback(masting_event, forcing!)

let
    u0 = [1.0, 1.5, 1.5]
    t_grid = range(0.0, t_end, length = 10000)
    t_start = 75.0

    # Chain
    par_chain = ModelPar(a_CP = 0.25, ω = 0.0)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain))
    sol_chain = solve(prob_chain, reltol = 1e-8, abstol = 1e-8)
    sol_chain_grid = sol_chain(t_grid)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = mast_event_times)
    sol_chain_mast = solve(prob_chain, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_chain_mast_grid = sol_chain_mast(t_grid)

    # Omnivory
    par_omn = ModelPar(a_CP = 0.25, ω = 0.1, pref = adapt_pref)

    prob_omn = ODEProblem(model!, u0, t_span, deepcopy(par_omn), tstops = mast_event_times)
    sol_omn = solve(prob_omn, reltol = 1e-8, abstol = 1e-8)
    sol_omn_grid = sol_omn(t_grid)

    prob_omn = ODEProblem(model!, u0, t_span, deepcopy(par_omn), tstops = mast_event_times)
    sol_omn_mast = solve(prob_omn, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_omn_mast_grid = sol_omn_mast(t_grid)

    ## Fixed preference omnivory
    par_fomn = ModelPar(a_CP = 0.25, ω = 0.1, pref = fixed_pref)

    prob_fomn = ODEProblem(model!, u0, t_span, deepcopy(par_fomn), tstops = mast_event_times)
    sol_fomn = solve(prob_fomn, reltol = 1e-8, abstol = 1e-8)
    sol_fomn_grid = sol_fomn(t_grid)

    prob_fomn = ODEProblem(model!, u0, t_span, deepcopy(par_fomn), tstops = mast_event_times)
    sol_fomn_mast = solve(prob_fomn, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_fomn_mast_grid = sol_fomn_mast(t_grid)

    # Layout
    fig = figure(figsize = (8, 12))
    R_col = "#1f77b4"
    C_col = "#ff7f0e"
    P_col = "#2ca02c"
    y_max = 5
    subplot(3, 2, 1)
    plot(sol_chain_grid.t, sol_chain_grid[1, :], color = R_col, alpha = 0.5)
    plot(sol_chain_grid.t, sol_chain_grid[2, :], color = C_col, alpha = 0.5)
    plot(sol_chain_grid.t, sol_chain_grid[3, :], color = P_col, alpha = 0.5)
    plot(sol_chain_mast_grid.t, sol_chain_mast_grid[1, :], color = R_col, label = "R")
    plot(sol_chain_mast_grid.t, sol_chain_mast_grid[2, :], color = C_col, label = "C")
    plot(sol_chain_mast_grid.t, sol_chain_mast_grid[3, :], color = P_col, label = "P")
    legend()
    title("Food Chain")
    xlim(t_start, t_end)
    ylim(0, y_max)

    subplot(3, 2, 2)
    plot(sol_chain_grid.t, [degree_omnivory(u, par_chain) for u in sol_chain_grid], color = "black")
    title("Degree of Omnivory [Chain]")
    xlim(t_start, t_end)

    subplot(3, 2, 3)
    plot(sol_fomn_grid.t, sol_fomn_grid[1, :], color = R_col, alpha = 0.5)
    plot(sol_fomn_grid.t, sol_fomn_grid[2, :], color = C_col, alpha = 0.5)
    plot(sol_fomn_grid.t, sol_fomn_grid[3, :], color = P_col, alpha = 0.5)
    plot(sol_fomn_mast_grid.t, sol_fomn_mast_grid[1, :], color = R_col, label = "R")
    plot(sol_fomn_mast_grid.t, sol_fomn_mast_grid[2, :], color = C_col, label = "C")
    plot(sol_fomn_mast_grid.t, sol_fomn_mast_grid[3, :], color = P_col, label = "P")
    legend()
    title("Passive Omnivory")
    xlim(t_start, t_end)
    ylim(0, y_max)

    subplot(3, 2, 4)
    plot(sol_fomn_mast_grid.t, [degree_omnivory(u, par_fomn) for u in sol_fomn_mast_grid], color = "black")
    plot(sol_fomn_grid.t, [degree_omnivory(u, par_fomn) for u in sol_fomn_grid], color = "black", alpha = 0.3)
    title("Degree of Omnivory [Passive]")
    xlim(t_start, t_end)

    subplot(3, 2, 5)
    plot(sol_omn_grid.t, sol_omn_grid[1, :], color = R_col, alpha = 0.5)
    plot(sol_omn_grid.t, sol_omn_grid[2, :], color = C_col, alpha = 0.5)
    plot(sol_omn_grid.t, sol_omn_grid[3, :], color = P_col, alpha = 0.5)
    plot(sol_omn_mast_grid.t, sol_omn_mast_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_mast_grid.t, sol_omn_mast_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_mast_grid.t, sol_omn_mast_grid[3, :], color = P_col, label = "P")
    legend()
    title("Adaptive Omnivory")
    xlim(t_start, t_end)
    ylim(0, y_max)

    subplot(3, 2, 6)
    plot(sol_omn_mast_grid.t, [degree_omnivory(u, par_omn) for u in sol_omn_mast_grid], color = "black")
    plot(sol_omn_grid.t, [degree_omnivory(u, par_omn) for u in sol_omn_grid], color = "black", alpha = 0.3)
    title("Degree of Omnivory [Adaptive]")
    xlim(t_start, t_end)

    tight_layout()
    return fig
end
