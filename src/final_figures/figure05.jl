include("basic_omnivory_module.jl")
using DifferentialEquations
using NLsolve
using QuadGK
using PyPlot

# # Resource Pulse/Masting Events
# I want to try to model the situation where we have a relatively low
# productivity environment, but we have period high productivity pulses,
# like a oak mastings etc.
#
# Objectives are that we can see a response of omnivory in the face of this
# variable productivity situation, and test that we are measuring omnivory
# in a way that is enlightening.
global first_mast = 100.0
global mast_freq = 100.0
mast_length = 2.0
global t_end = mast_freq * 2
global t_span = (0.0, t_end)

global mast_starts = first_mast:mast_freq:t_end
global mast_ends = mast_starts .+ mast_length
global mast_event_times = sort(union(mast_starts, mast_ends))

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

    # Eigenvalue analysis
    ## Chain
    eq_chain = find_eq(sol_chain[end], par_chain)
    chain_λ1 = λ1_stability(cmat(eq_chain, par_chain))
    chain_react = ν_stability(cmat(eq_chain, par_chain))
    
    ## Passive Omnivory
    eq_omn_fixed = find_eq(sol_omn_fixed[end], par_omn_fixed)
    omn_fixed_λ1 = λ1_stability(cmat(eq_omn_fixed, par_omn_fixed))
    omn_fixed_react = ν_stability(cmat(eq_omn_fixed, par_omn_fixed))
    
    ## Adaptive Omnivory
    #TODO: this naming for the par is bad
    eq_omn = find_eq(sol_omn[end], par_omn)
    omn_λ1 = λ1_stability(cmat(eq_omn, par_omn))
    omn_react = ν_stability(cmat(eq_omn, par_omn))
    
    # Measure of Overshoot
    ## What we are asking here is what is the total time * maginitute that the state variables are above or below the equilibrium after a perturbation
    chain_overshoot(t) = abs.(sol_chain_mast(t) .- eq_chain)
    omn_fixed_overshoot(t) = abs.(sol_omn_fixed_mast(t) .- eq_omn_fixed)
    omn_overshoot(t) = abs.(sol_omn_mast(t) .- eq_omn_mast)

    # let
    #     t_check = range(first_mast, first_mast + mast_freq, length = 1000)
    #     figure()
    #     subplot(3, 1, 1)
    #     plot(t_check, chain_overshoot.(t_check))
    #     subplot(3, 1, 2)
    #     plot(t_check, omn_fixed_overshoot.(t_check))
    #     subplot(3, 1, 3)
    #     plot(t_check, omn_mast_overshoot.(t_check))
    # end


    #NOTE: I am not yet dealing with looking for the peak of the Resource to start integration -- this is the naive approach of just using the start of the masting event
    chain_OS = [quadgk(t -> chain_overshoot(t)[1], first_mast, first_mast + mast_freq)[1],
                quadgk(t -> chain_overshoot(t)[2], first_mast, first_mast + mast_freq)[1],
                quadgk(t -> chain_overshoot(t)[3], first_mast, first_mast + mast_freq)[1]]
    
    omn_fixed_OS = [quadgk(t -> omn_fixed_overshoot(t)[1], first_mast, first_mast + mast_freq)[1],
                    quadgk(t -> omn_fixed_overshoot(t)[2], first_mast, first_mast + mast_freq)[1],
                    quadgk(t -> omn_fixed_overshoot(t)[3], first_mast, first_mast + mast_freq)[1]]
    
    omn_OS = [quadgk(t -> omn_overshoot(t)[1], first_mast, first_mast + mast_freq)[1],
                   quadgk(t -> omn_overshoot(t)[2], first_mast, first_mast + mast_freq)[1],
                   quadgk(t -> omn_overshoot(t)[3], first_mast, first_mast + mast_freq)[1]]
    
    
    # Layout
    fig = figure(figsize = (8, 9))
    R_col = "#1f77b4"
    C_col = "#ff7f0e"
    P_col = "#2ca02c"
    ## Time series axis limits
    y_max = 5

    ## First Column: Visualizations of Time Series
    ### Food Chain
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

    ### Omnivory [Fixed]
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
    
    ### Omnivory [Adaptive]
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
    
    ## Second Column: Dynamic Stability / Overshoot Metrics
    subplot(3, 2, 2)
    x_loc = [0, 1, 2]
    x_labels = ["FC", "P", "A"]
    plt.bar(x_loc, 1 ./ abs.([chain_λ1, omn_fixed_λ1, omn_λ1]))
    plt.xticks(x_loc, x_labels)
    ylabel("Local Return Time")
    
    subplot(3, 2, 4)
    g1 = [chain_OS[1], omn_fixed_OS[1], omn_OS[1]]
    g2 = [chain_OS[2], omn_fixed_OS[2], omn_OS[2]]
    g3 = [chain_OS[3], omn_fixed_OS[3], omn_OS[3]]

    # Set position of bar on X axis
    ## set width of bar
    bar_width = 0.25
    r1 = 1:length(g1)
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]
 
    # Make the plot
    plt.bar(r1, g1, width = bar_width, edgecolor = "white", label = "R")
    plt.bar(r2, g2, width = bar_width, edgecolor = "white", label = "C")
    plt.bar(r3, g3, width = bar_width, edgecolor = "white", label = "P")
    
    #TODO: Fix this up so we get meaninful x-axis labels
    # Add xticks on the middle of the group bars
    #plt.xlabel("State Variable", fontweight = "bold")
    #plt.xticks([r + bar_width for r in 0:(length(bars1))], ["R", "C", "P"])
    ylabel("Degree of Overshoot")
    
    # Create legend & Show graphic
    plt.legend()

    subplot(3, 2, 6)
    
    g1mm = [
        maximum(sol_chain_grid[1, :]) - minimum(sol_chain_grid[1, :]),
        maximum(sol_omn_fixed_grid[1, :]) - minimum(sol_omn_fixed_grid[1, :]),
        maximum(sol_omn_grid[1, :]) - minimum(sol_omn_grid[1, :])
    ]
    g2mm = [
        maximum(sol_chain_grid[2, :]) - minimum(sol_chain_grid[2, :]),
        maximum(sol_omn_fixed_grid[2, :]) - minimum(sol_omn_fixed_grid[2, :]),
        maximum(sol_omn_grid[2, :]) - minimum(sol_omn_grid[2, :])
    ]
    g3mm = [
        maximum(sol_chain_grid[3, :]) - minimum(sol_chain_grid[3, :]),
        maximum(sol_omn_fixed_grid[3, :]) - minimum(sol_omn_fixed_grid[3, :]),
        maximum(sol_omn_grid[3, :]) - minimum(sol_omn_grid[3, :])
    ]
        
    # Make the plot
    plt.bar(r1, g1mm, width = bar_width, edgecolor = "white", label = "R")
    plt.bar(r2, g2mm, width = bar_width, edgecolor = "white", label = "C")
    plt.bar(r3, g3mm, width = bar_width, edgecolor = "white", label = "P")
    ylabel("Max - Min")

    plt.legend()
            
    tight_layout()
    return fig
end
