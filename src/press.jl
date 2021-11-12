function press_unit(par, p_strength) 
    # Parameters
    u0 = [1.0, 1.5, 1.5]
    t_end = 5000 # way more than needed
    t_span = (0.0, t_end)
    t_start = 220
    p_start = 300.0
    len = 100000
    t_grid = range(t_start, t_end, length = len)
    t_press = range(p_start, t_end, length = len)
    Ω = par.Ω
    
    # functions
    press_event(u, t, integrator) = t ∈ p_start
    function forcing_press!(integrator)
        if integrator.t == p_start
            integrator.p.K = p_strength + integrator.p.K_base
        end
        return
    end
    cb_press = DiscreteCallback(press_event, forcing_press!)

    # run
    if par.pref == adapt_pref
        # compute w
        # NB: This may often increases computation time as we may compute
        # results for 'fixed' twice but the code is much easier to read. 
        # Here the time is not a major concern as it takes few minutes to run.
        par_fix = deepcopy(par)
        par_fix.pref = fixed_pref
        prob_fix = ODEProblem(model!, u0, t_span, deepcopy(par_fix),
            tstops = p_start)
        sol_fix = solve(prob_fix, reltol = 1e-8, abstol = 1e-8, 
            callback = cb_press)
        # sol_fix(t_grid)[1] starts at t_start
        eq = nlsolve(
            (du, u) -> model!(du, u, par_fix, 0.0), 
            sol_fix(t_grid)[1]
        ).zero
        ### Solving for ω we have `ω = Ω * C^* / (Ω * C^* + (1 - Ω) * R^*)`
        ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
        par.ω = ω
    end
    
    par_after = deepcopy(par)
    par_after.K = par.K + p_strength
    prob = ODEProblem(model!, u0, t_span, deepcopy(par), tstops = p_start)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_grid = sol(t_grid)
    t_hit_eq = find_times_hit_equil_press(sol(t_press))
    
    eq_befor = find_eq(sol_grid[1], par)
    eq_after = find_eq(sol_grid[end], par_after)
    λ1 = λ1_stability(cmat(eq_after, par_after))
    λ1_imag = λ1_stability_imag(cmat(eq_after, par_after))
    react = ν_stability(cmat(eq_after, par_after))

    # Measure of Overshoot
    predator_OS = overshoot(sol, eq_after, 3, t_hit_eq[3], t_end)
    consumer_OS = overshoot(sol, eq_after, 2, t_hit_eq[2], t_end)
    resource_OS = overshoot(sol, eq_after, 1, t_hit_eq[1], t_end)

    # Calculate max-min metric
    predator_mm = min_max(sol, 3, t_hit_eq[3], t_end)
    consumer_mm = min_max(sol, 2, t_hit_eq[2], t_end)
    resource_mm = min_max(sol, 1, t_hit_eq[1], t_end)

    push!(
        [],
        λ1,
        -1 ./ λ1,
        λ1_imag,
        react,
        predator_OS,
        consumer_OS, 
        resource_OS, 
        predator_mm, 
        consumer_mm, 
        resource_mm,  
        eq_befor, 
        eq_after, 
        par,
        sol
    )
    
end

    
function press(base_par, Ω, p_strength) 
    
    par_chain = deepcopy(base_par)
    par_chain.Ω = 0.0
    par_fixed = deepcopy(base_par)
    par_fixed.Ω = Ω
    par_fixed.pref = fixed_pref
    par_respo = deepcopy(par_fixed)
    par_respo.pref = adapt_pref
    
    push!(
        [],
        press_unit(par_chain, p_strength),
        press_unit(par_fixed, p_strength),
        press_unit(par_respo, p_strength)
    )

end










function press_old(K, a, e, m, )
    
    # Parameters
    u0 = [1.0, 1.5, 1.5]
    t_end = 2000
    t_span = (0.0, t_end)
    t_start = 220
    press_start = 300.0
    press_strength = 0.6
    t_grid = range(t_start, t_end, length = 100000)
    t_press = range(press_start, t_end, length = 100000)
    Ω = 0.1


        
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
    chain_press_hit_equil = find_times_hit_equil_press(sol_chain_press(t_press))

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
    omn_fixed_press_hit_equil = find_times_hit_equil_press(sol_omn_fixed_press(t_press))

    ## RESPONSIVE OMNIVORY
    ### Solve for ω so that at equlibrium Ω_fixed = Ω_adapt
    eq = nlsolve((du, u) -> model!(du, u, par_omn_fixed, 0.0), sol_omn_fixed_press_grid[1]).zero
    ### Solving for ω we have `ω = Ω * C^* / (Ω * C^* + (1 - Ω) * R^*)`
    ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
    par_omn_responsive = ModelPar(
        K_base = K,
        K = K,
        a_CP = a,
        e_CP = e,
        m_P = m,
        Ω = Ω,
        ω = ω, 
        pref = adapt_pref
    )
    par_omn_responsive_afterpress = deepcopy(par_omn_responsive)
    par_omn_responsive.K = par_omn_responsive.K
    
    prob_omn_responsive_press = ODEProblem(model!, u0, t_span, deepcopy(par_omn_responsive), tstops = press_start)
    sol_omn_responsive_press = solve(prob_omn_responsive_press, reltol = 1e-8, abstol = 1e-8, callback = cb_press)
    sol_omn_responsive_press_grid = sol_omn_responsive_press(t_grid)
    omn_responsive_press_hit_equil = find_times_hit_equil_press(sol_omn_responsive_press(t_press))


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


    # Measure of Overshoot
    chain_overshoot(t) = abs.(sol_chain_press(t) .- eq_chain_afterpress)
    omn_fixed_overshoot(t) = abs.(sol_omn_fixed_press(t) .- eq_omn_fixed_afterpress)
    omn_responsive_overshoot(t) = abs.(sol_omn_responsive_press(t) .- eq_omn_responsive_afterpress)

    function overshoot(spc, t_end)
        return [quadgk(t -> chain_overshoot(t)[spc], chain_press_hit_equil[spc], t_end)[1],
        quadgk(t -> omn_fixed_overshoot(t)[spc], omn_fixed_press_hit_equil[spc], t_end)[1],
        quadgk(t -> omn_responsive_overshoot(t)[spc], omn_responsive_press_hit_equil[spc], t_end)[1]]
    end
    
    println(chain_press_hit_equil)

    resource_OS = overshoot(1, t_end)
    consumer_OS = overshoot(2, t_end)
    predator_OS = overshoot(3, t_end)

    # Calculate max-min metric
    function min_max(spc, t_end)
        return [
        maximum(sol_chain_press(range(chain_press_hit_equil[spc], t_end, length = 10000))[spc, :]) - minimum(sol_chain_press(range(chain_press_hit_equil[spc], t_end, length = 10000))[spc, :]),
        maximum(sol_omn_fixed_press(range(omn_fixed_press_hit_equil[spc], t_end, length = 10000))[spc, :]) - minimum(sol_omn_fixed_press(range(omn_fixed_press_hit_equil[spc], t_end, length = 10000))[spc, :]),
        maximum(sol_omn_responsive_press(range(omn_responsive_press_hit_equil[spc], t_end, length = 10000))[spc, :]) - minimum(sol_omn_responsive_press(range(omn_responsive_press_hit_equil[spc], t_end, length = 10000))[spc, :])
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
    