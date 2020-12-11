include("basic_omnivory_module.jl")
using DifferentialEquations
using NLsolve
using QuadGK
using PyPlot

# # Press
press_start = 100.0
press_strength = 1.3

press_event(u, t, integrator) = t ∈ press_start

function forcing_press!(integrator)
    if integrator.t == press_start
        integrator.p.K = press_strength * integrator.p.K_base
    end
    return
end

cb_press = DiscreteCallback(press_event, forcing_press!)

let
    u0 = [1.0, 1.5, 1.5]
    t_end = 300
    t_span = (0.0, t_end)
    t_start = 75.0
    t_grid = range(t_start, t_end, length = 10000)
    t_press = range(press_start, t_end, length = 10000)

    # The global basic level of "Omnivory" we are looking at:
    Ω = 0.1

    # Chain
    par_chain = ModelPar(a_CP = 0.25, Ω = 0.0)
    par_chain_afterpress = ModelPar(K = 3.0 * press_strength, a_CP = 0.25, Ω = 0.0)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain))
    sol_chain = solve(prob_chain, reltol = 1e-8, abstol = 1e-8)
    sol_chain_grid = sol_chain(t_grid)

    prob_chain_press = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = press_start)
    sol_chain_press = solve(prob_chain_press, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_chain_press_grid = sol_chain_press(t_grid)

    ## Passive Omnivory
    par_omn_fixed = ModelPar(a_CP = 0.25, Ω = Ω, pref = fixed_pref)
    par_omn_fixed_afterpress = ModelPar(K = 3.0 *press_strength,a_CP = 0.25, Ω = Ω, pref = fixed_pref)

    prob_omn_fixed = ODEProblem(model!, u0, t_span, deepcopy(par_omn_fixed), tstops = press_start)
    sol_omn_fixed = solve(prob_omn_fixed, reltol = 1e-8, abstol = 1e-8)
    sol_omn_fixed_grid = sol_omn_fixed(t_grid)

    prob_omn_fixed_press = ODEProblem(model!, u0, t_span, deepcopy(par_omn_fixed), tstops = press_start)
    sol_omn_fixed_press = solve(prob_omn_fixed_press, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_omn_fixed_press_grid = sol_omn_fixed_press(t_grid)

    deg_omn_fixed_beforepress = round(degree_omnivory(sol_omn_fixed_press_grid.u[1], par_omn_fixed), digits = 2)
    deg_omn_fixed_afterpress = round(degree_omnivory(sol_omn_fixed_press_grid.u[end], par_omn_fixed_afterpress), digits = 2)
    max_deg_omn_fixed = round(maximum([degree_omnivory(u, par_omn_fixed) for u in sol_omn_fixed_press_grid]), digits = 2)
    # Responsive Omnivory
    ## Solve for ω so that at equlibrium Ω_fixed = Ω_adapt
    eq = nlsolve((du, u) -> model!(du, u, par_omn_fixed, 0.0), sol_omn_fixed[end]).zero

    ## Solving for ω we have `ω = Ω * C^* / (Ω * C^* + (1 - Ω) * R^*)`
    ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
    par_omn_responsive = ModelPar(a_CP = 0.25, Ω = Ω, ω = ω, pref = adapt_pref)
    par_omn_responsive_afterpress = ModelPar(K = 3.0 * press_strength, a_CP = 0.25, Ω = Ω, ω = ω, pref = adapt_pref)

    prob_omn_responsive = ODEProblem(model!, u0, t_span, deepcopy(par_omn_responsive), tstops = press_start)
    sol_omn_responsive = solve(prob_omn_responsive, reltol = 1e-8, abstol = 1e-8)
    sol_omn_responsive_grid = sol_omn_responsive(t_grid)

    prob_omn_responsive_press = ODEProblem(model!, u0, t_span, deepcopy(par_omn_responsive), tstops = press_start)
    sol_omn_responsive_press = solve(prob_omn_responsive, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_omn_responsive_press_grid = sol_omn_responsive_press(t_grid)

    deg_omn_responsive_beforepress = round(degree_omnivory(sol_omn_responsive_press_grid.u[1], par_omn_responsive), digits = 2)
    deg_omn_responsive_afterpress = round(degree_omnivory(sol_omn_responsive_press_grid.u[end], par_omn_responsive_afterpress), digits = 2)
    max_deg_omn_responsive = round(maximum([degree_omnivory(u, par_omn_responsive) for u in sol_omn_responsive_press_grid]), digits = 2)

    # println("$deg_omn_responsive_beforepress jjhjh")
    # Eigenvalue analysis
    ## Chain
    eq_chain = find_eq(sol_chain[end], par_chain_afterpress)
    chain_λ1 = λ1_stability(cmat(eq_chain, par_chain_afterpress))
    chain_react = ν_stability(cmat(eq_chain, par_chain_afterpress))

    ## Passive Omnivory
    eq_omn_fixed = find_eq(sol_omn_fixed[end], par_omn_fixed_afterpress)
    omn_fixed_λ1 = λ1_stability(cmat(eq_omn_fixed, par_omn_fixed_afterpress))
    omn_fixed_react = ν_stability(cmat(eq_omn_fixed, par_omn_fixed_afterpress))

    ## Responsive Omnivory
    eq_omn_responsive = find_eq(sol_omn_responsive[end], par_omn_responsive_afterpress)
    omn_responsive_λ1 = λ1_stability(cmat(eq_omn_responsive, par_omn_responsive_afterpress))
    omn_responsive_react = ν_stability(cmat(eq_omn_responsive, par_omn_responsive_afterpress))

    # Measure of Overshoot
    ## What we are asking here is what is the total time * maginitute that the state variables are above or below the equilibrium after a perturbation
    chain_overshoot(t) = abs.(sol_chain_press(t) .- eq_chain)
    omn_fixed_overshoot(t) = abs.(sol_omn_fixed_press(t) .- eq_omn_fixed)
    omn_responsive_overshoot(t) = abs.(sol_omn_responsive_press(t) .- eq_omn_responsive)

    # let
    #     t_check = range(press_start, t_end, length = 1000)
    #     figure()
    #     subplot(3, 1, 1)
    #     plot(t_check, chain_overshoot.(t_check))
    #     subplot(3, 1, 2)
    #     plot(t_check, omn_fixed_overshoot.(t_check))
    #     subplot(3, 1, 3)
    #     plot(t_check, omn_press_overshoot.(t_check))
    # end


    #NOTE: I am not yet dealing with looking for the peak of the Resource to start integration -- this is the naive approach of just using the start of the pressing event
    chain_OS = [quadgk(t -> chain_overshoot(t)[1], press_start, t_end)[1],
                quadgk(t -> chain_overshoot(t)[2], press_start, t_end)[1],
                quadgk(t -> chain_overshoot(t)[3], press_start, t_end)[1]]

    chain_OS_standardized = [chain_OS[1] / abs(sol_chain_grid[1,end] / eq_chain[1] ),
                             chain_OS[2] / abs(sol_chain_grid[2,end] / eq_chain[2] ),
                             chain_OS[3] / abs(sol_chain_grid[3,end] / eq_chain[3] )]

    omn_fixed_OS = [quadgk(t -> omn_fixed_overshoot(t)[1], press_start, t_end)[1],
                    quadgk(t -> omn_fixed_overshoot(t)[2], press_start, t_end)[1],
                    quadgk(t -> omn_fixed_overshoot(t)[3], press_start, t_end)[1]]

    omn_fixed_OS_standardized = [omn_fixed_OS[1] / abs(sol_omn_fixed_grid[1,end] / eq_omn_fixed[1] ),
                                 omn_fixed_OS[2] / abs(sol_omn_fixed_grid[2,end] / eq_omn_fixed[2] ),
                                 omn_fixed_OS[3] / abs(sol_omn_fixed_grid[3,end] / eq_omn_fixed[3] )]

    omn_responsive_OS = [quadgk(t -> omn_responsive_overshoot(t)[1], press_start, t_end)[1],
                   quadgk(t -> omn_responsive_overshoot(t)[2], press_start, t_end)[1],
                   quadgk(t -> omn_responsive_overshoot(t)[3], press_start, t_end)[1]]

    omn_responsive_OS_standardized = [omn_responsive_OS[1] / abs(sol_omn_responsive_grid[1,end] / eq_omn_responsive[1] ),
                                      omn_responsive_OS[2] / abs(sol_omn_responsive_grid[2,end] / eq_omn_responsive[2] ),
                                      omn_responsive_OS[3] / abs(sol_omn_responsive_grid[3,end] / eq_omn_responsive[3] )]
    #Calculate max-min metric
    resource_mm = [
        maximum(sol_chain_press(t_press)[1,:]) - minimum(sol_chain_press(t_press)[1,:]),
        maximum(sol_omn_fixed_press(t_press)[1,:]) - minimum(sol_omn_fixed_press(t_press)[1,:]),
        maximum(sol_omn_responsive_press(t_press)[1,:]) - minimum(sol_omn_responsive_press(t_press)[1,:])
    ]
    consumer_mm = [
        maximum(sol_chain_press(t_press)[2,:]) - minimum(sol_chain_press(t_press)[2,:]),
        maximum(sol_omn_fixed_press(t_press)[2,:]) - minimum(sol_omn_fixed_press(t_press)[2,:]),
        maximum(sol_omn_responsive_press(t_press)[2,:]) - minimum(sol_omn_responsive_press(t_press)[2,:])
    ]
    predator_mm = [
        maximum(sol_chain_press(t_press)[3,:]) - minimum(sol_chain_press(t_press)[3,:]),
        maximum(sol_omn_fixed_press(t_press)[3,:]) - minimum(sol_omn_fixed_press(t_press)[3,:]),
        maximum(sol_omn_responsive_press(t_press)[3,:]) - minimum(sol_omn_responsive_press(t_press)[3,:])
    ]

    # Layout
    fig = figure(figsize = (8, 9))
    R_col = "#1f77b4"
    C_col = "#ff7f0e"
    P_col = "#2ca02c"
    ## Time series axis limits
    y_max = 4

    ## First Column: Visualizations of Time Series
    ### Food Chain
    subplot(3, 2, 1)
    hlines(eq_chain[1], 100, 300, color = R_col, alpha = 0.5)
    hlines(eq_chain[2], 100, 300, color = C_col, alpha = 0.5)
    hlines(eq_chain[3], 100, 300, color = P_col, alpha = 0.5)
    plot(sol_chain_press_grid.t, sol_chain_press_grid[1, :], color = R_col, label = "R")
    plot(sol_chain_press_grid.t, sol_chain_press_grid[2, :], color = C_col, label = "C")
    plot(sol_chain_press_grid.t, sol_chain_press_grid[3, :], color = P_col, label = "P")
    legend()
    title("Food Chain")
    xlim(t_start, t_end)
    ylim(0, y_max)

    ### Omnivory [Fixed]
    subplot(3, 2, 3)
    hlines(eq_omn_fixed[1], 100, 300, color = R_col, alpha = 0.5)
    hlines(eq_omn_fixed[2], 100, 300, color = C_col, alpha = 0.5)
    hlines(eq_omn_fixed[3], 100, 300, color = P_col, alpha = 0.5)
    plot(sol_omn_fixed_press_grid.t, sol_omn_fixed_press_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_fixed_press_grid.t, sol_omn_fixed_press_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_fixed_press_grid.t, sol_omn_fixed_press_grid[3, :], color = P_col, label = "P")
    annotate("° Omn\n $deg_omn_fixed_beforepress", (75.0,3.5), xycoords = "data", fontsize = 10)
    annotate("° Omn\n $deg_omn_fixed_afterpress", (270.0,3.5), xycoords = "data", fontsize = 10)
    annotate("Max ° Omn\n $max_deg_omn_fixed", (180.0,3.5), xycoords = "data", fontsize = 10)
    legend()
    title("Omnivory [Passive]")
    xlim(t_start, t_end)
    ylim(0, y_max)

    ### Omnivory [Adaptive]
    subplot(3, 2, 5)
    hlines(eq_omn_responsive[1], 100, 300, color = R_col, alpha = 0.5)
    hlines(eq_omn_responsive[2], 100, 300, color = C_col, alpha = 0.5)
    hlines(eq_omn_responsive[3], 100, 300, color = P_col, alpha = 0.5)
    plot(sol_omn_responsive_press_grid.t, sol_omn_responsive_press_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_responsive_press_grid.t, sol_omn_responsive_press_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_responsive_press_grid.t, sol_omn_responsive_press_grid[3, :], color = P_col, label = "P")
    annotate("° Omn\n $deg_omn_responsive_beforepress", (75.0,0.1), xycoords = "data", fontsize = 10)
    annotate("° Omn\n $deg_omn_responsive_afterpress", (270.0,0.1), xycoords = "data", fontsize = 10)
    annotate("Max ° Omn\n $max_deg_omn_responsive", (180.0,0.1), xycoords = "data", fontsize = 10)
    legend()
    title("Omnivory [Adaptive]")
    xlim(t_start, t_end)
    ylim(0, y_max)

    ## Second Column: Dynamic Stability / Overshoot Metrics
    subplot(3, 2, 2)
    x_loc = [0, 1, 2]
    x_labels = ["Food\nChain", "Passive", "Responsive"]
    plt.bar(x_loc, 1 ./ abs.([chain_λ1, omn_fixed_λ1, omn_responsive_λ1]))
    plt.xticks(x_loc, x_labels)
    plt.tick_params(axis = "x", which = "both", length=0)
    ylabel("Local Return Time")

    subplot(3, 2, 4)
    # g1 = [chain_OS[1], omn_fixed_OS[1], omn_responsive_OS[1]]
    # g2 = [chain_OS[2], omn_fixed_OS[2], omn_responsive_OS[2]]
    # g3 = [chain_OS[3], omn_fixed_OS[3], omn_responsive_OS[3]]

    g1 = [chain_OS_standardized[1], omn_fixed_OS_standardized[1], omn_responsive_OS_standardized[1]]
    g2 = [chain_OS_standardized[2], omn_fixed_OS_standardized[2], omn_responsive_OS_standardized[2]]
    g3 = [chain_OS_standardized[3], omn_fixed_OS_standardized[3], omn_responsive_OS_standardized[3]]

    #Set position of bar on X axis
    # set width of bar
    bar_width = 0.25
    r1 = 1:length(g1)
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]

    # Make the plot
    plt.bar(r1, g1, width = bar_width, edgecolor = "white", label = "R")
    plt.bar(r2, g2, width = bar_width, edgecolor = "white", label = "C")
    plt.bar(r3, g3, width = bar_width, edgecolor = "white", label = "P")
    plt.xticks(r2, x_labels)
    plt.tick_params(axis = "x", which = "both", length=0)
    ylabel("Standardized\nDegree of Overshoot")

    # Create legend & Show graphic
    plt.legend()

    subplot(3, 2, 6)

    plt.bar(r1, resource_mm, width = bar_width, edgecolor = "white", label = "R")
    plt.bar(r2, consumer_mm, width = bar_width, edgecolor = "white", label = "C")
    plt.bar(r3, predator_mm, width = bar_width, edgecolor = "white", label = "P")
    plt.xticks(r2, x_labels)
    plt.tick_params(axis = "x", which = "both", length=0)
    ylabel("Max - Min")

    plt.legend()

    tight_layout()
    return fig
end
