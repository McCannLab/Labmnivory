include("basic_omnivory_module.jl")
using NLsolve
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

    # The global basic level of "Omnivory" we are looking at:
    Ω = 0.1

    # Chain
    par_chain = ModelPar(a_CP = 0.25, Ω = 0.0)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain))
    sol_chain = solve(prob_chain, reltol = 1e-8, abstol = 1e-8)
    sol_chain_grid = sol_chain(t_grid)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = mast_event_times)
    sol_chain_mast = solve(prob_chain, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_chain_mast_grid = sol_chain_mast(t_grid)

    ## Fixed preference omnivory
    par_omn_fixed = ModelPar(a_CP = 0.25, Ω = Ω, pref = fixed_pref)

    prob_omn_fixed = ODEProblem(model!, u0, t_span, deepcopy(par_omn_fixed), tstops = mast_event_times)
    sol_omn_fixed = solve(prob_omn_fixed, reltol = 1e-8, abstol = 1e-8)
    sol_omn_fixed_grid = sol_omn_fixed(t_grid)

    prob_omn_fixed_mast = ODEProblem(model!, u0, t_span, deepcopy(par_omn_fixed), tstops = mast_event_times)
    sol_omn_fixed_mast = solve(prob_omn_fixed_mast, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_omn_fixed_mast_grid = sol_omn_fixed_mast(t_grid)

    # Omnivory
    ## Solve for ω so that at equlibrium Ω_fixed = Ω_adapt
    eq = nlsolve((du, u) -> model!(du, u, par_omn_fixed, 0.0), sol_omn_fixed[end]).zero
    ## Solving for ω we have `ω = Ω * C^* / (Ω * C^* + (1 - Ω) * R^*)`
    ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
    par_omn = ModelPar(a_CP = 0.25, Ω = Ω, ω = ω, pref = adapt_pref)

    prob_omn = ODEProblem(model!, u0, t_span, deepcopy(par_omn), tstops = mast_event_times)
    sol_omn = solve(prob_omn, reltol = 1e-8, abstol = 1e-8)
    sol_omn_grid = sol_omn(t_grid)

    prob_omn = ODEProblem(model!, u0, t_span, deepcopy(par_omn), tstops = mast_event_times)
    sol_omn_mast = solve(prob_omn, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_omn_mast_grid = sol_omn_mast(t_grid)

    # Layout
    fig = figure(figsize = (8, 9))
    R_col = "#1f77b4"
    C_col = "#ff7f0e"
    P_col = "#2ca02c"
    ## Time series axis limits
    y_max = 5

    ## Food Chain
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
    ylim(0, 1)

    ## Omnivory [Fixed]
    subplot(3, 2, 3)
    plot(sol_omn_fixed_grid.t, sol_omn_fixed_grid[1, :], color = R_col, alpha = 0.5)
    plot(sol_omn_fixed_grid.t, sol_omn_fixed_grid[2, :], color = C_col, alpha = 0.5)
    plot(sol_omn_fixed_grid.t, sol_omn_fixed_grid[3, :], color = P_col, alpha = 0.5)
    plot(sol_omn_fixed_mast_grid.t, sol_omn_fixed_mast_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_fixed_mast_grid.t, sol_omn_fixed_mast_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_fixed_mast_grid.t, sol_omn_fixed_mast_grid[3, :], color = P_col, label = "P")
    legend()
    title("Omnivory [Passive]")
    xlim(t_start, t_end)
    ylim(0, y_max)

    subplot(3, 2, 4)
    plot(sol_omn_fixed_mast_grid.t, [degree_omnivory(u, par_omn_fixed) for u in sol_omn_fixed_mast_grid], color = "black")
    plot(sol_omn_fixed_grid.t, [degree_omnivory(u, par_omn_fixed) for u in sol_omn_fixed_grid], color = "black", alpha = 0.3)
    title("Degree of Omnivory [Passive]")
    xlim(t_start, t_end)
    ylim(0, 1)

    ## Omnivory [Adaptive]
    subplot(3, 2, 5)
    plot(sol_omn_grid.t, sol_omn_grid[1, :], color = R_col, alpha = 0.5)
    plot(sol_omn_grid.t, sol_omn_grid[2, :], color = C_col, alpha = 0.5)
    plot(sol_omn_grid.t, sol_omn_grid[3, :], color = P_col, alpha = 0.5)
    plot(sol_omn_mast_grid.t, sol_omn_mast_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_mast_grid.t, sol_omn_mast_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_mast_grid.t, sol_omn_mast_grid[3, :], color = P_col, label = "P")
    legend()
    title("Omnivory [Adaptive]")
    xlim(t_start, t_end)
    ylim(0, y_max)

    subplot(3, 2, 6)
    plot(sol_omn_mast_grid.t, [degree_omnivory(u, par_omn) for u in sol_omn_mast_grid], color = "black")
    plot(sol_omn_grid.t,
    [degree_omnivory(u, par_omn) for u in sol_omn_grid], color = "black", alpha = 0.3)
    title("Degree of Omnivory [Adaptive]")
    xlim(t_start, t_end)
    ylim(0, 1)

    tight_layout()
    return fig
end
