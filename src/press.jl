using DifferentialEquations, NLsolve

# Full analysis for one press simulation (1 system)
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
            # addition for press
            integrator.p.K = p_strength + integrator.p.K_base
        end
        return
    end
    cb_press = DiscreteCallback(press_event, forcing_press!)

    # run fixed case first for active omnivory
    if par.pref == adapt_pref
        # compute ω
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

# Perform the pulse for the 3 systems: 
# - Food chain: "chain";C OP (fixed)
# - Passive omnivory: "passive";
# - Active omnivory: "active";
function press(base_par, Ω, p_strength) 
    
    par_chain = deepcopy(base_par)
    par_chain.Ω = 0.0
    par_passive = deepcopy(base_par)
    par_passive.Ω = Ω
    par_passive.pref = fixed_pref
    par_active = deepcopy(par_passive)
    par_active.pref = adapt_pref
    
    push!(
        [],
        press_unit(par_chain, p_strength),
        press_unit(par_passive, p_strength),
        press_unit(par_active, p_strength)
    )

end



