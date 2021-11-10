function press(K, a, e, m)
    
    # Parameters
    u0 = [1.0, 1.5, 1.5]
    t_end = 10000
    t_span = (0.0, t_end)
    t_start = 220
    press_start = 300.0
    press_strength = .5
    t_grid = range(t_start, t_end, length = 10000)
    t_press = range(press_start, t_end, length = 10000)
    Ω = 0.1
    
    function find_times_hit_equil(data)
        eq = data[1, end], data[2, end], data[3, end]
        times = zeros(3)
        for animal in 1:3
            for i in 20:length(data)
                # cannot be too strict here otherwise the value of the 
                # first ht time varies a lort which will have serious 
                # impact on min and max (overshoot too) leading to major oscillations 
                if isapprox(data[animal,i], eq[animal], atol = 0.05)
                    times[animal] = data.t[i]
                    break
                end
            end
        end
        return times
    end

    
    press_event(u, t, integrator) = t ∈ press_start

    function forcing_press!(integrator)
        if integrator.t == press_start
            integrator.p.K = press_strength + integrator.p.K_base
        end
        return
    end

    cb_press = DiscreteCallback(press_event, forcing_press!)

    
    # ODE
    ## FOOD CHAIN
    par_chain = ModelPar(
        a_CP = a,
        e_CP = e,
        m_P = m,
        Ω = 0.0
    )
    par_chain.K_base = par_chain.K = K
    par_chain_afterpress = ModelPar(
        K = par_chain.K + press_strength, 
        a_CP = a,
        e_CP = e,
        m_P = m,
        Ω = 0.0
    )
    prob_chain_press = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = press_start)
    sol_chain_press = solve(prob_chain_press, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_chain_press_grid = sol_chain_press(t_grid)
    chain_press_hit_equil = find_times_hit_equil(sol_chain_press(t_press))

    ## PASSIVE OMNIVORY
    par_omn_fixed = ModelPar(
        a_CP = a,
        e_CP = e,
        m_P = m,
        Ω = Ω,
        pref = fixed_pref
    )
    par_omn_fixed.K_base = par_omn_fixed.K = K
    par_omn_fixed_afterpress = ModelPar(
        K = par_omn_fixed.K + press_strength, 
        a_CP = a,
        e_CP = e,
        m_P = m,
        Ω = Ω, 
        pref = fixed_pref)
    prob_omn_fixed_press = ODEProblem(model!, u0, t_span, deepcopy(par_omn_fixed), tstops = press_start)
    sol_omn_fixed_press = solve(prob_omn_fixed_press, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_omn_fixed_press_grid = sol_omn_fixed_press(t_grid)
    omn_fixed_press_hit_equil = find_times_hit_equil(sol_omn_fixed_press(t_press))

    ## RESPONSIVE OMNIVORY
    ### Solve for ω so that at equlibrium Ω_fixed = Ω_adapt
    eq = nlsolve((du, u) -> model!(du, u, par_omn_fixed, 0.0), sol_omn_fixed_press_grid[1]).zero
    ### Solving for ω we have `ω = Ω * C^* / (Ω * C^* + (1 - Ω) * R^*)`
    ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
    par_omn_responsive = ModelPar(
        a_CP = a,
        e_CP = e,
        m_P = m,
        Ω = Ω,
        ω = ω, 
        pref = adapt_pref
    )
    par_omn_responsive.K_base = par_omn_responsive.K = K
    par_omn_responsive_afterpress = ModelPar(
        K = par_omn_responsive.K + press_strength, 
        a_CP = a,
        e_CP = e,
        m_P = m,
        Ω = Ω, 
        ω = ω, 
        pref = adapt_pref
    )
    prob_omn_responsive_press = ODEProblem(model!, u0, t_span, deepcopy(par_omn_responsive), tstops = press_start)
    sol_omn_responsive_press = solve(prob_omn_responsive_press, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_omn_responsive_press_grid = sol_omn_responsive_press(t_grid)
    omn_responsive_press_hit_equil = find_times_hit_equil(sol_omn_responsive_press(t_press))


    # EIGENVALUE ANALYSIS
    ## FOOD CHAIN
    eq_chain_afterpress = find_eq(sol_chain_press_grid[end], par_chain_afterpress)
    chain_λ1 = λ1_stability(cmat(eq_chain_afterpress, par_chain_afterpress))
    chain_λ1_imag = λ1_stability_imag(cmat(eq_chain_afterpress, par_chain_afterpress))
    chain_react = ν_stability(cmat(eq_chain_afterpress, par_chain_afterpress))

    ## PASSIVE OMNIVORY
    eq_omn_fixed_afterpress = find_eq(sol_omn_fixed_press_grid[end], par_omn_fixed_afterpress)
    omn_fixed_λ1 = λ1_stability(cmat(eq_omn_fixed_afterpress, par_omn_fixed_afterpress))
    omn_fixed_λ1_imag = λ1_stability_imag(cmat(eq_omn_fixed_afterpress, par_omn_fixed_afterpress))
    omn_fixed_react = ν_stability(cmat(eq_omn_fixed_afterpress, par_omn_fixed_afterpress))

    ## RESPONSIVE OMNIVORY
    eq_omn_responsive_afterpress = find_eq(sol_omn_responsive_press_grid[1], par_omn_responsive_afterpress)
    omn_responsive_λ1 = λ1_stability(cmat(eq_omn_responsive_afterpress, par_omn_responsive_afterpress))
    omn_responsive_λ1_imag = λ1_stability_imag(cmat(eq_omn_responsive_afterpress, par_omn_responsive_afterpress))
    omn_responsive_react = ν_stability(cmat(eq_omn_responsive_afterpress, par_omn_responsive_afterpress))

    # Find times for each model where Resource maximized after press
    # res_max_times = [
    #     find_time_hit_res_max(sol_chain_press_grid), 
    #     find_time_hit_res_max(sol_omn_fixed_press_grid), 
    #     find_time_hit_res_max(sol_omn_responsive_press_grid)
    # ]

    # Measure of Overshoot
    chain_overshoot(t) = abs.(sol_chain_press(t) .- eq_chain_afterpress)
    omn_fixed_overshoot(t) = abs.(sol_omn_fixed_press(t) .- eq_omn_fixed_afterpress)
    omn_responsive_overshoot(t) = abs.(sol_omn_responsive_press(t) .- eq_omn_responsive_afterpress)

    function overshoot(animal, t_end)
        return [quadgk(t -> chain_overshoot(t)[animal], chain_press_hit_equil[animal], t_end)[1],
        quadgk(t -> omn_fixed_overshoot(t)[animal], omn_fixed_press_hit_equil[animal], t_end)[1],
        quadgk(t -> omn_responsive_overshoot(t)[animal], omn_responsive_press_hit_equil[animal], t_end)[1]]
    end
    
    println(chain_press_hit_equil)

    resource_OS = overshoot(1, t_end)
    consumer_OS = overshoot(2, t_end)
    predator_OS = overshoot(3, t_end)

    # Calculate max-min metric
    function min_max(animal, t_end)
        return [
        maximum(sol_chain_press(range(chain_press_hit_equil[animal], t_end, length = 10000))[animal, :]) - minimum(sol_chain_press(range(chain_press_hit_equil[animal], t_end, length = 10000))[animal, :]),
        maximum(sol_omn_fixed_press(range(omn_fixed_press_hit_equil[animal], t_end, length = 10000))[animal, :]) - minimum(sol_omn_fixed_press(range(omn_fixed_press_hit_equil[animal], t_end, length = 10000))[animal, :]),
        maximum(sol_omn_responsive_press(range(omn_responsive_press_hit_equil[animal], t_end, length = 10000))[animal, :]) - minimum(sol_omn_responsive_press(range(omn_responsive_press_hit_equil[animal], t_end, length = 10000))[animal, :])
        ]
    end

    resource_mm = min_max(1, t_end)
    consumer_mm = min_max(2, t_end)
    predator_mm = min_max(3, t_end)

    push!(
        [],
        [chain_λ1, omn_fixed_λ1, omn_responsive_λ1],
        [-1 ./ chain_λ1, -1 ./ omn_fixed_λ1, -1 ./ omn_responsive_λ1],
        [chain_λ1_imag, omn_fixed_λ1_imag, omn_responsive_λ1_imag],
        [chain_react, omn_fixed_react, omn_responsive_react],
        resource_OS,
        consumer_OS,
        predator_OS,
        resource_mm,
        consumer_mm,
        predator_mm
    )
    
end
    