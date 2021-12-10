using DifferentialEquations, NLsolve

# Entire analysis for one pulse simulation (1 system)
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
            # multiplication for pulse
            integrator.p.K = p_strength * integrator.p.K_base
        elseif integrator.t == pulse_end
            integrator.p.K = integrator.p.K_base
        end
        return
    end
    cb_pulse = DiscreteCallback(pulse_event, forcing_pulse!)

    # run
    if par.pref == adapt_pref
        # compute ω
        # NB: This may often increases computation time as we may compute
        # results for passive omnivory twice but the code is much easier to
        # read. Here the time is not a major concern as it takes few minutes to
        # run.
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

# Perform the pulse for the 3 systems: 
# - Food chain: "chain";C OP (fixed)
# - Passive omnivory: "passive";
# - Active omnivory: "active";
function pulse(base_par, Ω, p_length, p_strength) 
    
    par_chain = deepcopy(base_par)
    par_chain.Ω = 0.0
    par_passive = deepcopy(base_par)
    par_passive.Ω = Ω
    par_passive.pref = fixed_pref
    par_active = deepcopy(par_passive)
    par_active.pref = adapt_pref
    
    push!(
        [],
        pulse_unit(par_chain, p_length, p_strength),
        pulse_unit(par_passive, p_length, p_strength),
        pulse_unit(par_active, p_length, p_strength)
    )

end
