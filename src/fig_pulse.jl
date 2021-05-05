include("basic_omnivory_module.jl")
using DifferentialEquations, NLsolve, QuadGK, PyPlot
pygui(true)

# SINGLE PULSE EVENT
# Objectives are that we can see a response of omnivory in the face of this
# variable productivity situation, and test that we are measuring omnivory
# in a way that is enlightening.


#TODO Clean and make metrics start calculating after resource peak

global pulse_length = 2.0
global pulse_start = 200
global pulse_end = pulse_start + pulse_length
global pulse_event_times = union(pulse_start, pulse_end)
pulse_strength = 2.0

pulse_event(u, t, integrator) = t ∈ pulse_event_times

function forcing!(integrator)
    if integrator.t == pulse_start
        integrator.p.K = pulse_strength * integrator.p.K_base
    elseif integrator.t == pulse_end
        integrator.p.K = integrator.p.K_base
    end
    return
end

cb = DiscreteCallback(pulse_event, forcing!)

function find_time_hit_res_max(data)
    max_res = findmax(data[1,:])
    return data.t[max_res[2]]
end



let
    # initial values + time setup
    u0 = [1.0, 1.5, 1.5]
    t_end = 350
    t_span = (0.0, t_end)
    t_start = 175.0
    t_grid = range(t_start, t_end, length = 10000)
    # The global basic level of "Omnivory" we are looking at:
    Ω = 0.1

    # ODE
    ## FOOD CHAIN
    par_chain = ModelPar(Ω = 0.0)
    prob_chain_pulse = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = pulse_event_times)
    sol_chain_pulse = solve(prob_chain_pulse, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_chain_pulse_grid = sol_chain_pulse(t_grid)

    ## PASSIVE OMNIVORY
    par_omn_fixed = ModelPar(Ω = Ω, pref = fixed_pref)
    prob_omn_fixed_pulse = ODEProblem(model!, u0, t_span, deepcopy(par_omn_fixed), tstops = pulse_event_times)
    sol_omn_fixed_pulse = solve(prob_omn_fixed_pulse, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_omn_fixed_pulse_grid = sol_omn_fixed_pulse(t_grid)

    ## RESPONSIVE OMNIVORY
    ### Solve for ω so that at equlibrium Ω_fixed = Ω_adapt
    eq = nlsolve((du, u) -> model!(du, u, par_omn_fixed, 0.0), sol_omn_fixed_pulse_grid[1]).zero
    ### Solving for ω we have `ω = Ω * C^* / (Ω * C^* + (1 - Ω) * R^*)`
    ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
    par_omn_responsive = ModelPar(Ω = Ω, ω = ω, pref = adapt_pref)
    prob_omn_responsive_pulse = ODEProblem(model!, u0, t_span, deepcopy(par_omn_responsive), tstops = pulse_event_times)
    sol_omn_responsive_pulse = solve(prob_omn_responsive_pulse, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_omn_responsive_pulse_grid = sol_omn_responsive_pulse(t_grid)


    # EIGENVALUE ANALYSIS
    ## FOOD CHAIN
    eq_chain = find_eq(sol_chain_pulse_grid[1], par_chain)
    chain_λ1 = λ1_stability(cmat(eq_chain, par_chain))
    # chain_react = ν_stability(cmat(eq_chain, par_chain)) CAN DELETE?

    ## PASSIVE OMNIVORY
    eq_omn_fixed = find_eq(sol_omn_fixed_pulse_grid[1], par_omn_fixed)
    omn_fixed_λ1 = λ1_stability(cmat(eq_omn_fixed, par_omn_fixed))
    # omn_fixed_react = ν_stability(cmat(eq_omn_fixed, par_omn_fixed)) CAN DELETE?

    ## RESPONSIVE OMNIVORY
    eq_omn_responsive = find_eq(sol_omn_responsive_pulse_grid[1], par_omn_responsive)
    omn_responsive_λ1 = λ1_stability(cmat(eq_omn_responsive, par_omn_responsive))
    # omn_responsive_react = ν_stability(cmat(eq_omn_responsive, par_omn_responsive)) CAN DELETE?

    #Find times for each model where Resource maximized after pulse
    res_max_times = [find_time_hit_res_max(sol_chain_pulse_grid), find_time_hit_res_max(sol_omn_fixed_pulse_grid), find_time_hit_res_max(sol_omn_responsive_pulse_grid)]

    # Measure of Overshoot
    ## What we are asking here is what is the total time * maginitute that the state variables are above or below the equilibrium after a perturbation
    chain_overshoot(t) = abs.(sol_chain_pulse(t) .- eq_chain)
    omn_fixed_overshoot(t) = abs.(sol_omn_fixed_pulse(t) .- eq_omn_fixed)
    omn_responsive_overshoot(t) = abs.(sol_omn_responsive_pulse(t) .- eq_omn_responsive)

    function overshoot(animal, t_end)
        return [quadgk(t -> chain_overshoot(t)[animal], 100.0, t_end)[1],
        quadgk(t -> omn_fixed_overshoot(t)[animal], 100.0, t_end)[1],
        quadgk(t -> omn_responsive_overshoot(t)[animal], 100.0, t_end)[1]]
    end

    resource_OS = overshoot(1, t_end)
    consumer_OS = overshoot(2, t_end)
    predator_OS = overshoot(3, t_end)

    #Calculate max-min metric
    function min_max(animal, t_end)
        return [maximum(sol_chain_pulse(range(res_max_times[animal], t_end, length = 10000))[animal, :]) - minimum(sol_chain_pulse(range(res_max_times[animal], t_end, length = 10000))[animal, :]),
        maximum(sol_omn_fixed_pulse(range(res_max_times[animal], t_end, length = 10000))[animal, :]) - minimum(sol_omn_fixed_pulse(range(res_max_times[animal], t_end, length = 10000))[animal, :]),
        maximum(sol_omn_responsive_pulse(range(res_max_times[animal], t_end, length = 10000))[animal, :]) - minimum(sol_omn_responsive_pulse(range(res_max_times[animal], t_end, length = 10000))[animal, :])
            ]
    end

    resource_mm = min_max(1, t_end)
    consumer_mm = min_max(2, t_end)
    predator_mm = min_max(3, t_end)
    #Min max
    # g1mm = [
    #     maximum(sol_chain_grid[1, :]) - minimum(sol_chain_grid[1, :]),
    #     maximum(sol_omn_fixed_grid[1, :]) - minimum(sol_omn_fixed_grid[1, :]),
    #     maximum(sol_omn_responsive_grid[1, :]) - minimum(sol_omn_responsive_grid[1, :])
    # ]
    # g2mm = [
    #     maximum(sol_chain_grid[2, :]) - minimum(sol_chain_grid[2, :]),
    #     maximum(sol_omn_fixed_grid[2, :]) - minimum(sol_omn_fixed_grid[2, :]),
    #     maximum(sol_omn_responsive_grid[2, :]) - minimum(sol_omn_responsive_grid[2, :])
    # ]
    # g3mm = [
    #     maximum(sol_chain_grid[3, :]) - minimum(sol_chain_grid[3, :]),
    #     maximum(sol_omn_fixed_grid[3, :]) - minimum(sol_omn_fixed_grid[3, :]),
    #     maximum(sol_omn_responsive_grid[3, :]) - minimum(sol_omn_responsive_grid[3, :])
    # ]

    # Layout
    fig = figure(figsize = (8, 9))
    R_col = "#1f77b4"
    C_col = "#ff7f0e"
    P_col = "#2ca02c"
    ## Time series axis limits
    y_max = 5
    
    
    # MAX DEGREE OF OMNIVORY PER PHASE 
    ## Fixed omnivory
    ### OmEq
    sol = sol_omn_fixed_pulse(180:0.01:199)
    println(
        "Fixed - Equilibrium: ", 
        maximum([degree_omnivory(sol[:,i], par_omn_fixed) for i in 1:size(sol)[2]])
    )
    ### Pulse
    sol = sol_omn_fixed_pulse(200:0.01:205)
    println(
        "Fixed - Pulse: ", 
        maximum([degree_omnivory(sol[:,i], par_omn_fixed) for i in 1:size(sol)[2]])
    )
    ### OmT
    sol = sol_omn_fixed_pulse(205:0.01:325)
    println(
        "Fixed - Transient: ", 
        maximum([degree_omnivory(sol[:,i], par_omn_fixed) for i in 1:size(sol)[2]])
    )

    ## Responsive omnivory
    ### Equilibrium
    sol = sol_omn_responsive_pulse(180:0.01:199)
    println(
        "Responsive - Equilibrium: ", 
        maximum([degree_omnivory(sol[:,i], par_omn_responsive) for i in 1:size(sol)[2]])
    )
    ### OmB 
    sol = sol_omn_responsive_pulse(200:0.01:205)
    println(
        "Responsive - Pulse: ", 
        maximum([degree_omnivory(sol[:,i], par_omn_responsive) for i in 1:size(sol)[2]])
    )
    ### OmT
    sol = sol_omn_responsive_pulse(205:0.01:310)
    println(
        "Responsive - Transient: ", 
        maximum([degree_omnivory(sol[:,i], par_omn_responsive) for i in 1:size(sol)[2]])
    )


    ## First Column: Visualizations of Time Series
    ### Food Chain
    subplot(3, 2, 1)
    hlines(eq_chain[1], t_start, t_end, color = R_col, alpha = 0.5)
    hlines(eq_chain[2], t_start, t_end, color = C_col, alpha = 0.5)
    hlines(eq_chain[3], t_start, t_end, color = P_col, alpha = 0.5)
    plot(sol_chain_pulse_grid.t, sol_chain_pulse_grid[1, :], color = R_col, label = "R")
    plot(sol_chain_pulse_grid.t, sol_chain_pulse_grid[2, :], color = C_col, label = "C")
    plot(sol_chain_pulse_grid.t, sol_chain_pulse_grid[3, :], color = P_col, label = "P")
    legend()
    title("Food Chain")
    xlim(t_start, t_end)
    ylim(0, y_max)
    ylabel("Density")

    ### Omnivory [Fixed]
    subplot(3, 2, 3)
    hlines(eq_omn_fixed[1], t_start, t_end, color = R_col, alpha = 0.5)
    hlines(eq_omn_fixed[2], t_start, t_end, color = C_col, alpha = 0.5)
    hlines(eq_omn_fixed[3], t_start, t_end, color = P_col, alpha = 0.5)
    plot(sol_omn_fixed_pulse_grid.t, sol_omn_fixed_pulse_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_fixed_pulse_grid.t, sol_omn_fixed_pulse_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_fixed_pulse_grid.t, sol_omn_fixed_pulse_grid[3, :], color = P_col, label = "P")
    legend()
    title("Passive Omnivory")
    xlim(t_start, t_end)
    ylim(0, y_max)
    ylabel("Density")

    ### Omnivory [responsive]
    subplot(3, 2, 5)
    hlines(eq_omn_responsive[1], t_start, t_end, color = R_col, alpha = 0.5)
    hlines(eq_omn_responsive[2], t_start, t_end, color = C_col, alpha = 0.5)
    hlines(eq_omn_responsive[3], t_start, t_end, color = P_col, alpha = 0.5)
    plot(sol_omn_responsive_pulse_grid.t, sol_omn_responsive_pulse_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_responsive_pulse_grid.t, sol_omn_responsive_pulse_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_responsive_pulse_grid.t, sol_omn_responsive_pulse_grid[3, :], color = P_col, label = "P")
    legend()
    title("Responsive Omnivory")
    xlim(t_start, t_end)
    ylim(0, y_max)
    ylabel("Density")
    xlabel("Time")

    ## Second Column: Dynamic Stability / Overshoot Metrics
    subplot(3, 2, 2)
    x_loc = [0, 1, 2]
    x_labels = ["FC", "P", "R"]
    plt.bar(x_loc, 1 ./ abs.([chain_λ1, omn_fixed_λ1, omn_responsive_λ1]))
    plt.xticks(x_loc, x_labels)
    plt.tick_params(axis = "x", which = "both", length = 0)
    ylabel("Local Return Time")

    subplot(3, 2, 4)
    # Set position of bar on X axis
    ## set width of bar
    bar_width = 0.25
    r1 = 1:length(resource_mm)
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]

    # Make the plot
    plt.bar(r1, resource_OS, width = bar_width, edgecolor = "white", label = "R")
    plt.bar(r2, consumer_OS, width = bar_width, edgecolor = "white", label = "C")
    plt.bar(r3, predator_OS, width = bar_width, edgecolor = "white", label = "P")
    plt.xticks(r2, x_labels)
    plt.tick_params(axis = "x", which = "both", length=0)
    # ylim(0, 11)
    ylabel("Degree of Overshoot")

    # Create legend & Show graphic
    plt.legend()

    subplot(3, 2, 6)

    # Make the plot
    plt.bar(r1, resource_mm, width = bar_width, edgecolor = "white", label = "R")
    plt.bar(r2, consumer_mm, width = bar_width, edgecolor = "white", label = "C")
    plt.bar(r3, predator_mm, width = bar_width, edgecolor = "white", label = "P")
    plt.xticks(r2, x_labels)
    plt.tick_params(axis = "x", which = "both", length=0)
    ylabel("Max - Min")

    plt.legend()

    tight_layout()
    savefig("figs/fig3A.svg")
    return fig
end
