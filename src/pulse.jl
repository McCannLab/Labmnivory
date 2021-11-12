using DifferentialEquations, NLsolve

function find_time_hit_res_max(res)
    max_res = findmax(res[1,:])
    return res.t[max_res[2]]
end


function pulse_unit(par, p_length, p_strength) 
    # Parameters   
    u0 = [1.0, 1.5, 1.5]
    pulse_start = 200
    pulse_end = pulse_start + p_length
    pulse_event_times = union(pulse_start, pulse_end)
    t_end = 350
    t_span = (0.0, t_end)
    t_start = 175.0
    len = 100000
    t_grid = range(t_start, t_end, length = len)
    t_after = range(pulse_end, t_end, length = len)
    Ω = par.Ω
    
    # functions
    pulse_event(u, t, integrator) = t ∈ pulse_event_times
    # NB: for pulse p_strength is a multiplicator
    function forcing_pulse!(integrator)
        if integrator.t == pulse_start
            integrator.p.K = p_strength * integrator.p.K_base
        elseif integrator.t == pulse_end
            integrator.p.K = integrator.p.K_base
        end
        return
    end
    cb_pulse = DiscreteCallback(pulse_event, forcing_pulse!)

    # run
    if par.pref == adapt_pref
        # compute w
        # NB: This may often increases computation time as we may compute
        # results for 'fixed' twice but the code is much easier to read. 
        # Here the time is not a major concern as it takes few minutes to run.
        par_fix = deepcopy(par)
        par_fix.pref = fixed_pref
        prob_fix = ODEProblem(model!, u0, t_span, deepcopy(par_fix),
            tstops = pulse_event_times)
        sol_fix = solve(prob_fix, reltol = 1e-8, abstol = 1e-8, 
            callback = cb_pulse)
        # sol_fix(t_grid)[1] starts at t_start
        eq = nlsolve(
            (du, u) -> model!(du, u, par_fix, 0.0), 
            sol_fix(t_grid)[1]
        ).zero
        ### Solving for ω we have `ω = Ω * C^* / (Ω * C^* + (1 - Ω) * R^*)`
        ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
        par.ω = ω
    end
    
    par_pulse = deepcopy(par)
    par_pulse.K = p_strength * par.K_base
    prob = ODEProblem(model!, u0, t_span, deepcopy(par), tstops = pulse_event_times)
    sol = solve(prob, reltol = 1e-8, abstol = 1e-8, callback = cb_pulse)
    sol_grid = sol(t_grid)
    t_hit_eq = find_times_hit_equil_press(sol(t_after))
    # println(t_hit_eq)
    
    # Same equilibrium before and after
    eq_before = find_eq(sol_grid[1], par)

    # evaluate on 
    λ1 = λ1_stability(cmat(eq_before, par))
    λ1_imag = λ1_stability_imag(cmat(eq_before, par))
    react = ν_stability(cmat(eq_before, par))

    predator_OS = overshoot(sol, eq_before, 3, t_hit_eq[3], t_end)
    consumer_OS = overshoot(sol, eq_before, 2, t_hit_eq[2], t_end)
    resource_OS = overshoot(sol, eq_before, 1, t_hit_eq[1], t_end)

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
        eq_before,  # for consistency in output size with `press()`
        eq_before,   
        par,
        sol
    )
    
end


function pulse(base_par, Ω, p_length, p_strength) 
    
    par_chain = deepcopy(base_par)
    par_chain.Ω = 0.0
    par_fixed = deepcopy(base_par)
    par_fixed.Ω = Ω
    par_fixed.pref = fixed_pref
    par_respo = deepcopy(par_fixed)
    par_respo.pref = adapt_pref
    
    push!(
        [],
        pulse_unit(par_chain, p_length, p_strength),
        pulse_unit(par_fixed, p_length, p_strength),
        pulse_unit(par_respo, p_length, p_strength)
    )

end


function pulse_old(K, a, e, m)
    
    # Parameters
    u0 = [1.0, 1.5, 1.5]
    t_end = 350
    t_span = (0.0, t_end)
    t_start = 175.0
    t_grid = range(pulse_end, t_end, length = 10000)
    t_after = range(pulse_end, t_end, length = 10000)
    Ω = 0.1
    
    
    
    # ODE
    ## FOOD CHAIN
    par_chain = ModelPar(Ω = 0.0)
    par_chain.K_base = par_chain.K = K
    par_chain.a_CP = a
    par_chain.e_CP = e
    par_chain.m_P = m
    prob_chain_pulse = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = pulse_event_times)
    sol_chain_pulse = solve(prob_chain_pulse, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_chain_pulse_grid = sol_chain_pulse(t_grid)
    chain_pulse_hit_equil = find_times_hit_equil(sol_chain_pulse(t_after))

    ## PASSIVE OMNIVORY
    par_omn_fixed = ModelPar(Ω = Ω, pref = fixed_pref)
    par_omn_fixed.K_base = par_omn_fixed.K = K
    par_omn_fixed.a_CP = a
    par_omn_fixed.e_CP = e
    par_omn_fixed.m_P = m
    prob_omn_fixed_pulse = ODEProblem(model!, u0, t_span, deepcopy(par_omn_fixed), tstops = pulse_event_times)
    sol_omn_fixed_pulse = solve(prob_omn_fixed_pulse, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_omn_fixed_pulse_grid = sol_omn_fixed_pulse(t_grid)
    fixed_pulse_hit_equil = find_times_hit_equil(sol_omn_fixed_pulse(t_after))

    ## RESPONSIVE OMNIVORY
    ### Solve for ω so that at equlibrium Ω_fixed = Ω_adapt
    eq = nlsolve((du, u) -> model!(du, u, par_omn_fixed, 0.0), sol_omn_fixed_pulse_grid[1]).zero
    ### Solving for ω we have `ω = Ω * C^* / (Ω * C^* + (1 - Ω) * R^*)`
    ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
    par_omn_responsive = ModelPar(Ω = Ω, ω = ω, pref = adapt_pref)
    par_omn_responsive.K_base = par_omn_responsive.K = K
    par_omn_responsive.a_CP = a
    par_omn_responsive.e_CP = e
    par_omn_responsive.m_P = m
    prob_omn_responsive_pulse = ODEProblem(model!, u0, t_span, deepcopy(par_omn_responsive), tstops = pulse_event_times)
    sol_omn_responsive_pulse = solve(prob_omn_responsive_pulse, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_omn_responsive_pulse_grid = sol_omn_responsive_pulse(t_grid)
    responsive_pulse_hit_equil = find_times_hit_equil(sol_omn_responsive_pulse(t_after))


    # EIGENVALUE ANALYSIS
    ## FOOD CHAIN
    eq_chain = find_eq(sol_chain_pulse_grid[1], par_chain)
    chain_λ1 = λ1_stability(cmat(eq_chain, par_chain))
    chain_λ1_imag = λ1_stability_imag(cmat(eq_chain, par_chain))
    chain_react = ν_stability(cmat(eq_chain, par_chain))

    ## PASSIVE OMNIVORY
    eq_omn_fixed = find_eq(sol_omn_fixed_pulse_grid[1], par_omn_fixed)
    omn_fixed_λ1 = λ1_stability(cmat(eq_omn_fixed, par_omn_fixed))
    omn_fixed_λ1_imag = λ1_stability_imag(cmat(eq_omn_fixed, par_omn_fixed))
    omn_fixed_react = ν_stability(cmat(eq_omn_fixed, par_omn_fixed))

    ## RESPONSIVE OMNIVORY
    eq_omn_responsive = find_eq(sol_omn_responsive_pulse_grid[1], par_omn_responsive)
    omn_responsive_λ1 = λ1_stability(cmat(eq_omn_responsive, par_omn_responsive))
    omn_responsive_λ1_imag = λ1_stability_imag(cmat(eq_omn_responsive, par_omn_responsive))
    omn_responsive_react = ν_stability(cmat(eq_omn_responsive, par_omn_responsive))

    # Find times for each model where Resource maximized after pulse
    res_max_times = [
        find_time_hit_res_max(sol_chain_pulse_grid), 
        find_time_hit_res_max(sol_omn_fixed_pulse_grid), 
        find_time_hit_res_max(sol_omn_responsive_pulse_grid)
    ]

    # Measure of Overshoot
    chain_overshoot(t) = abs.(sol_chain_pulse(t) .- eq_chain)
    omn_fixed_overshoot(t) = abs.(sol_omn_fixed_pulse(t) .- eq_omn_fixed)
    omn_responsive_overshoot(t) = abs.(sol_omn_responsive_pulse(t) .- eq_omn_responsive)

    function overshoot(animal, t_end)
        return [quadgk(t -> chain_overshoot(t)[animal], chain_pulse_hit_equil[animal], t_end)[1],
        quadgk(t -> omn_fixed_overshoot(t)[animal], fixed_pulse_hit_equil[animal], t_end)[1],
        quadgk(t -> omn_responsive_overshoot(t)[animal], responsive_pulse_hit_equil[animal], t_end)[1]]
    end

    resource_OS = overshoot(1, t_end)
    consumer_OS = overshoot(2, t_end)
    predator_OS = overshoot(3, t_end)

    # Calculate max-min metric
    function min_max(animal, t_end)
        return [maximum(sol_chain_pulse(range(chain_pulse_hit_equil[animal], t_end, length = 10000))[animal,:]) - minimum(sol_chain_pulse(range(chain_pulse_hit_equil[animal], t_end, length = 10000))[animal,:]),
        maximum(sol_omn_fixed_pulse(range(fixed_pulse_hit_equil[animal], t_end, length = 10000))[animal,:]) - minimum(sol_omn_fixed_pulse(range(fixed_pulse_hit_equil[animal], t_end, length = 10000))[animal,:]),
        maximum(sol_omn_responsive_pulse(range(responsive_pulse_hit_equil[animal], t_end, length = 10000))[animal,:]) - minimum(sol_omn_responsive_pulse(range(responsive_pulse_hit_equil[animal], t_end, length = 10000))[animal,:])
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
    