include("basic_omnivory_module.jl")
include("top_heavy.jl")
include("asynchrony.jl")
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
mast_length = 4.0
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
    dt = 0.01
    t_grid = 0.0:dt:t_end
    t_start = 88.0
    wind = 500  # 100 = 1 time unit
    t_wind = (wind * 0.5 * dt):dt:(t_end - wind * 0.5 * dt)
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
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    labs = ["R", "C", "P"]
    ## Time series axis limits
    y_max = 5


    # FIG 3A
    ## Food Chain
    subplot(3, 1, 1)
    for i in 1:3
        plot(sol_chain_grid.t, sol_chain_grid[i, :], color = colors[i], alpha = 0.5)
        plot(sol_chain_mast_grid.t, sol_chain_mast_grid[i, :], color = colors[i], label = labs[i])
    end
    legend()
    title("Food Chain")
    xlim(t_start, t_end)
    ylim(0, y_max)
    ## Omnivory [Fixed]
    subplot(3, 1, 2)
    for i in 1:3
        plot(sol_omn_fixed_grid.t, sol_omn_fixed_grid[i, :], color = colors[i], alpha = 0.5)
        plot(sol_omn_fixed_mast_grid.t, sol_omn_fixed_mast_grid[i, :], color = colors[i], label = labs[i])
    end
    title("Omnivory [Passive]")
    xlim(t_start, t_end .- 5)
    ylim(0, y_max)
    ## Omnivory [Responsive]
    subplot(3, 1, 3)
    for i in 1:3
        plot(sol_omn_grid.t, sol_omn_grid[i, :], color = colors[i], alpha = 0.5)
        plot(sol_omn_mast_grid.t, sol_omn_mast_grid[i, :], color = colors[i], label = labs[i])
    end
    title("Omnivory [Responsive]")
    xlim(t_start, t_end .- 5)
    ylim(0, y_max)



    # FIG 3B
    subplot(5, 1, 2)
    plot(sol_chain_grid.t, [top_heaviness(u, sol_chain_mast_grid[4001]) for u in sol_chain_mast_grid], color = "black", linestyle = "--", label = "Chain")
    plot(sol_chain_grid.t, [top_heaviness(u, sol_omn_fixed_mast_grid[4001]) for u in sol_omn_fixed_mast_grid], color = "black", label = "Passive")
    plot(sol_chain_grid.t, [top_heaviness(u, sol_omn_mast_grid[4001]) for u in sol_omn_mast_grid], color = "grey", label = "Active")
    title("Top heaviness")
    xlim(t_start, t_end .- 5)
    ylim(0.5, 2)

    subplot(5, 1, 1)
    plot(sol_omn_fixed_mast_grid.t, [degree_omnivory(u, par_omn_fixed) for u in sol_omn_fixed_mast_grid], color = "black", label = "Passive")
    plot(sol_omn_mast_grid.t, [degree_omnivory(u, par_omn) for u in sol_omn_mast_grid], color = "grey", label = "[Responsive")
    title("Degree of Omnivory")
    xlim(t_start, t_end .- 5)
    ylim(0, 1)

    subplot(5, 1, 3)
    plot(t_wind, cv(sol_chain_mast_grid[3, :], wind), color = "black", linestyle = "--", label = "Chain")
    plot(t_wind, cv(sol_omn_fixed_mast_grid[3, :], wind), color = "black", label = "Passive")
    plot(t_wind, cv(sol_omn_mast_grid[3, :], wind), color = "grey", label = "Responsive")
    title("CV of C")
    xlim(t_start, t_end .- 5)
    ylim(-0.01, 0.13)
    # global_cv(sol_chain_mast_grid)[3]
    # global_cv(sol_omn_fixed_mast_grid)[3]
    # global_cv(sol_chain_mast_grid)[3]
    legend()

    subplot(5, 1, 4)
    plot(t_wind, asynchrony(sol_chain_mast_grid, wind), color = "black", linestyle = "--", label = "Chain")
    plot(t_wind, asynchrony(sol_omn_fixed_mast_grid, wind), color = "black", label = "Passive")
    plot(t_wind, asynchrony(sol_omn_mast_grid, wind), color = "grey", label = "Responsive")
    title("Asynchrony C-R")
    xlim(t_start, t_end .- 5)
    ylim(-0.3, 0.3)
    # global_asyn(sol_chain_mast_grid)
    # global_asyn(sol_omn_fixed_mast_grid)
    # global_asyn(sol_chain_mast_grid)

    subplot(5, 1, 5)
    plot(t_wind, auc(sol_chain_mast_grid[3, :], sol_chain_grid[3, :], dt, wind), color = "black", linestyle = "--", label = "Chain")
    plot(t_wind, auc(sol_omn_fixed_mast_grid[3, :], sol_omn_fixed_grid[3, :], dt, wind), color = "black", label = "Passive")
    plot(t_wind, auc(sol_omn_mast_grid[3, :], sol_omn_grid[3, :], dt, wind), color = "grey", label = "Responsive")
    title("AUC C")
    xlim(t_start, t_end .- 5)
    ylim(-0.1, 3)
    # global_auc(sol_chain_mast_grid, sol_chain_grid, dt)[3]
    # global_auc(sol_omn_fixed_mast_grid, sol_omn_fixed_grid, dt)[3]
    # global_auc(sol_omn_mast_grid, sol_omn_grid, dt)[3]

    # tight_layout()
    return fig
end
