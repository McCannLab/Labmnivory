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


function find_times_hit_equil(data)
    eq = data[1,end], data[2,end], data[3,end]
    times = zeros(3)
    for animal in 1:3
        for i in 20:length(data)
            if isapprox(data[animal,i],eq[animal], atol = 0.001)
                times[animal] = data.t[i]
                break
            end
        end
    end
    return times
end


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

    prob_chain_press = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = press_start)
    sol_chain_press = solve(prob_chain_press, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_chain_press_grid = sol_chain_press(t_grid)
    chain_press_hit_equil = find_times_hit_equil(sol_chain_press(t_press))

    ## Passive Omnivory
    par_omn_fixed = ModelPar(a_CP = 0.25, Ω = Ω, pref = fixed_pref)
    par_omn_fixed_afterpress = ModelPar(K = 3.0 *press_strength,a_CP = 0.25, Ω = Ω, pref = fixed_pref)

    prob_omn_fixed_press = ODEProblem(model!, u0, t_span, deepcopy(par_omn_fixed), tstops = press_start)
    sol_omn_fixed_press = solve(prob_omn_fixed_press, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_omn_fixed_press_grid = sol_omn_fixed_press(t_grid)
    fixed_press_hit_equil = find_times_hit_equil(sol_omn_fixed_press(t_press))

    deg_omn_fixed_beforepress = round(degree_omnivory(sol_omn_fixed_press_grid.u[1], par_omn_fixed), digits = 2)
    deg_omn_fixed_afterpress = round(degree_omnivory(sol_omn_fixed_press_grid.u[end], par_omn_fixed_afterpress), digits = 2)
    max_deg_omn_fixed = round(maximum([degree_omnivory(u, par_omn_fixed) for u in sol_omn_fixed_press_grid]), digits = 2)

    # Responsive Omnivory
    ## Solve for ω so that at equlibrium Ω_fixed = Ω_adapt
    eq = nlsolve((du, u) -> model!(du, u, par_omn_fixed, 0.0), sol_omn_fixed_press_grid[1]).zero

    ## Solving for ω we have `ω = Ω * C^* / (Ω * C^* + (1 - Ω) * R^*)`
    ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
    par_omn_responsive = ModelPar(a_CP = 0.25, Ω = Ω, ω = ω, pref = adapt_pref)
    par_omn_responsive_afterpress = ModelPar(K = 3.0 * press_strength, a_CP = 0.25, Ω = Ω, ω = ω, pref = adapt_pref)

    prob_omn_responsive_press = ODEProblem(model!, u0, t_span, deepcopy(par_omn_responsive), tstops = press_start)
    sol_omn_responsive_press = solve(prob_omn_responsive_press, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_omn_responsive_press_grid = sol_omn_responsive_press(t_grid)
    responsive_press_hit_equil = find_times_hit_equil(sol_omn_responsive_press(t_press))

    deg_omn_responsive_beforepress = round(degree_omnivory(sol_omn_responsive_press_grid.u[1], par_omn_responsive), digits = 2)
    deg_omn_responsive_afterpress = round(degree_omnivory(sol_omn_responsive_press_grid.u[end], par_omn_responsive_afterpress), digits = 2)
    max_deg_omn_responsive = round(maximum([degree_omnivory(u, par_omn_responsive) for u in sol_omn_responsive_press_grid]), digits = 2)

    # Eigenvalue analysis
    ## Chain
    eq_chain_afterpress = find_eq(sol_chain_press[end], par_chain_afterpress)
    chain_λ1 = λ1_stability(cmat(eq_chain_afterpress, par_chain_afterpress))
    # chain_react = ν_stability(cmat(eq_chain, par_chain_afterpress)) CAN DELETE?

    ## Passive Omnivory
    eq_omn_fixed_afterpress = find_eq(sol_omn_fixed_press[end], par_omn_fixed_afterpress)
    omn_fixed_λ1 = λ1_stability(cmat(eq_omn_fixed_afterpress, par_omn_fixed_afterpress))
    # omn_fixed_react = ν_stability(cmat(eq_omn_fixed, par_omn_fixed_afterpress)) CAN DELETE?

    ## Responsive Omnivory
    eq_omn_responsive_afterpress = find_eq(sol_omn_responsive_press[end], par_omn_responsive_afterpress)
    omn_responsive_λ1 = λ1_stability(cmat(eq_omn_responsive_afterpress, par_omn_responsive_afterpress))
    # omn_responsive_react = ν_stability(cmat(eq_omn_responsive, par_omn_responsive_afterpress)) CAN DELETE?

    # Measure of Overshoot
    ## What we are asking here is what is the total time * maginitute that the state variables are above or below the equilibrium after a perturbation
    chain_overshoot(t) = abs.(sol_chain_press(t) .- eq_chain_afterpress)
    omn_fixed_overshoot(t) = abs.(sol_omn_fixed_press(t) .- eq_omn_fixed_afterpress)
    omn_responsive_overshoot(t) = abs.(sol_omn_responsive_press(t) .- eq_omn_responsive_afterpress)

    function overshoot(animal, t_end)
        return [quadgk(t -> chain_overshoot(t)[animal], chain_press_hit_equil[animal], t_end)[1],
        quadgk(t -> omn_fixed_overshoot(t)[animal], fixed_press_hit_equil[animal], t_end)[1],
        quadgk(t -> omn_responsive_overshoot(t)[animal], responsive_press_hit_equil[animal], t_end)[1]]
    end

    resource_OS = overshoot(1, t_end)
    consumer_OS = overshoot(2, t_end)
    predator_OS = overshoot(3, t_end)

    #Calculate max-min metric
    function min_max(animal, t_end)
        return [maximum(sol_chain_press(range(chain_press_hit_equil[animal], t_end, length = 10000))[animal,:]) - minimum(sol_chain_press(range(chain_press_hit_equil[animal], t_end, length = 10000))[animal,:]),
        maximum(sol_omn_fixed_press(range(fixed_press_hit_equil[animal], t_end, length = 10000))[animal,:]) - minimum(sol_omn_fixed_press(range(fixed_press_hit_equil[animal], t_end, length = 10000))[animal,:]),
        maximum(sol_omn_responsive_press(range(responsive_press_hit_equil[animal], t_end, length = 10000))[animal,:]) - minimum(sol_omn_responsive_press(range(responsive_press_hit_equil[animal], t_end, length = 10000))[animal,:])
            ]
    end

    resource_mm = min_max(1, t_end)
    consumer_mm = min_max(2, t_end)
    predator_mm = min_max(3, t_end)

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
    hlines(eq_chain_afterpress[1], 100, 300, color = R_col, alpha = 0.5)
    hlines(eq_chain_afterpress[2], 100, 300, color = C_col, alpha = 0.5)
    hlines(eq_chain_afterpress[3], 100, 300, color = P_col, alpha = 0.5)
    plot(sol_chain_press_grid.t, sol_chain_press_grid[1, :], color = R_col, label = "R")
    plot(sol_chain_press_grid.t, sol_chain_press_grid[2, :], color = C_col, label = "C")
    plot(sol_chain_press_grid.t, sol_chain_press_grid[3, :], color = P_col, label = "P")
    legend()
    title("Food Chain")
    xlim(t_start, t_end)
    ylim(0, y_max)
    ylabel("Density")

    ### Omnivory [Fixed]
    subplot(3, 2, 3)
    hlines(eq_omn_fixed_afterpress[1], 100, 300, color = R_col, alpha = 0.5)
    hlines(eq_omn_fixed_afterpress[2], 100, 300, color = C_col, alpha = 0.5)
    hlines(eq_omn_fixed_afterpress[3], 100, 300, color = P_col, alpha = 0.5)
    plot(sol_omn_fixed_press_grid.t, sol_omn_fixed_press_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_fixed_press_grid.t, sol_omn_fixed_press_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_fixed_press_grid.t, sol_omn_fixed_press_grid[3, :], color = P_col, label = "P")
    # annotate("° Omn\n $deg_omn_fixed_beforepress", (75.0,3.5), xycoords = "data", fontsize = 10)
    # annotate("° Omn\n $deg_omn_fixed_afterpress", (270.0,3.5), xycoords = "data", fontsize = 10)
    # annotate("Max ° Omn\n $max_deg_omn_fixed", (180.0,3.5), xycoords = "data", fontsize = 10)
    legend()
    title("Passive Omnivory")
    xlim(t_start, t_end)
    ylim(0, y_max)
    ylabel("Density")

    ### Omnivory [Adaptive]
    subplot(3, 2, 5)
    hlines(eq_omn_responsive_afterpress[1], 100, 300, color = R_col, alpha = 0.5)
    hlines(eq_omn_responsive_afterpress[2], 100, 300, color = C_col, alpha = 0.5)
    hlines(eq_omn_responsive_afterpress[3], 100, 300, color = P_col, alpha = 0.5)
    plot(sol_omn_responsive_press_grid.t, sol_omn_responsive_press_grid[1, :], color = R_col, label = "R")
    plot(sol_omn_responsive_press_grid.t, sol_omn_responsive_press_grid[2, :], color = C_col, label = "C")
    plot(sol_omn_responsive_press_grid.t, sol_omn_responsive_press_grid[3, :], color = P_col, label = "P")
    # annotate("° Omn\n $deg_omn_responsive_beforepress", (75.0,0.1), xycoords = "data", fontsize = 10)
    # annotate("° Omn\n $deg_omn_responsive_afterpress", (270.0,0.1), xycoords = "data", fontsize = 10)
    # annotate("Max ° Omn\n $max_deg_omn_responsive", (180.0,0.1), xycoords = "data", fontsize = 10)
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
    #Set position of bar on X axis
    # set width of bar
    bar_width = 0.25
    r1 = 1:length(resource_OS)
    r2 = [x + bar_width for x in r1]
    r3 = [x + bar_width for x in r2]

    # Make the plot
    plt.bar(r1, resource_OS, width = bar_width, edgecolor = "white", label = "R")
    plt.bar(r2, consumer_OS, width = bar_width, edgecolor = "white", label = "C")
    plt.bar(r3, predator_OS, width = bar_width, edgecolor = "white", label = "P")
    plt.xticks(r2, x_labels)
    plt.tick_params(axis = "x", which = "both", length = 0)
    ylabel("Degree of Overshoot")

    # Create legend & Show graphic
    plt.legend()

    subplot(3, 2, 6)

    plt.bar(r1, resource_mm, width = bar_width, edgecolor = "white", label = "R")
    plt.bar(r2, consumer_mm, width = bar_width, edgecolor = "white", label = "C")
    plt.bar(r3, predator_mm, width = bar_width, edgecolor = "white", label = "P")
    plt.xticks(r2, x_labels)
    plt.tick_params(axis = "x", which = "both", length = 0)
    ylabel("Max - Min")

    plt.legend()

    tight_layout()
    savefig("figs/fig3b.svg")
    return fig

end
