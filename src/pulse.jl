function pulse(K, a, e, m)
    
    u0 = [1.0, 1.5, 1.5]
    t_grid = range(0.0, t_end, length = 10000)
    t_start = 75.0

    # The global basic level of "Omnivory" we are looking at:
    Ω = 0.1
    
    par_chain = ModelPar(Ω = 0.0)
    par_chain.K_base = par_chain.K = K
    par_chain.a_CP = a
    par_chain.e_CP = e
    par_chain.m_P = m

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain))
    sol_chain = solve(prob_chain, reltol = 1e-8, abstol = 1e-8)
    sol_chain_grid = sol_chain(t_grid)

    prob_chain = ODEProblem(model!, u0, t_span, deepcopy(par_chain), tstops = mast_event_times)
    sol_chain_mast = solve(prob_chain, reltol = 1e-8, abstol = 1e-8, callback = cb)
    sol_chain_mast_grid = sol_chain_mast(t_grid)

    
    ## Fixed preference omnivory
    par_omn_fixed = ModelPar(Ω = Ω, pref = fixed_pref)
    par_chain.K_base = par_omn_fixed.K = K
    par_omn_fixed.a_CP = a
    par_omn_fixed.e_CP = e
    par_omn_fixed.m_P = m
    
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
    par_omn = ModelPar(Ω = Ω, ω = ω, pref = adapt_pref)
    par_chain.K_base = par_omn.K = K
    par_omn.a_CP = a
    par_omn.e_CP = e
    par_omn.m_P = m
    
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
    chain_λ1_imag = λ1_stability_imag(cmat(eq_chain, par_chain))
    chain_react = ν_stability(cmat(eq_chain, par_chain))
    
    ## Passive Omnivory
    eq_omn_fixed = find_eq(sol_omn_fixed[end], par_omn_fixed)
    omn_fixed_λ1 = λ1_stability(cmat(eq_omn_fixed, par_omn_fixed))
    omn_fixed_λ1_imag = λ1_stability_imag(cmat(eq_omn_fixed, par_omn_fixed))
    omn_fixed_react = ν_stability(cmat(eq_omn_fixed, par_omn_fixed))
    
    ## Adaptive Omnivory
    #TODO: this naming for the par is bad
    eq_omn = find_eq(sol_omn[end], par_omn)
    omn_λ1 = λ1_stability(cmat(eq_omn, par_omn))
    omn_λ1_imag = λ1_stability_imag(cmat(eq_omn, par_omn))
    omn_react = ν_stability(cmat(eq_omn, par_omn))
    
    # Measure of Overshoot
    ## What we are asking here is what is the total time * maginitute that the state variables are above or below the equilibrium after a perturbation
    chain_overshoot(t) = abs.(sol_chain_mast(t) .- eq_chain)
    omn_fixed_overshoot(t) = abs.(sol_omn_fixed_mast(t) .- eq_omn_fixed)
    omn_overshoot(t) = abs.(sol_omn(t) .- eq_omn)
    
    ##
    chain_OS = [quadgk(t -> chain_overshoot(t)[1], first_mast, first_mast + mast_freq)[1],
                quadgk(t -> chain_overshoot(t)[2], first_mast, first_mast + mast_freq)[1],
                quadgk(t -> chain_overshoot(t)[3], first_mast, first_mast + mast_freq)[1]
                ]
    
    omn_fixed_OS = [quadgk(t -> omn_fixed_overshoot(t)[1], first_mast, first_mast + mast_freq)[1],
                    quadgk(t -> omn_fixed_overshoot(t)[2], first_mast, first_mast + mast_freq)[1],
                    quadgk(t -> omn_fixed_overshoot(t)[3], first_mast, first_mast + mast_freq)[1]]
    
    omn_OS = [quadgk(t -> omn_overshoot(t)[1], first_mast, first_mast + mast_freq)[1],
              quadgk(t -> omn_overshoot(t)[2], first_mast, first_mast + mast_freq)[1],
              quadgk(t -> omn_overshoot(t)[3], first_mast, first_mast + mast_freq)[1]]
    
    
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
    
    push!(
        [],
        [chain_λ1, omn_fixed_λ1, omn_λ1],
        [chain_react, omn_fixed_react, omn_react],
        chain_OS,
        omn_fixed_OS,
        omn_OS,
        g1mm,
        g2mm,
        g3mm,
        [chain_λ1_imag, omn_fixed_λ1_imag, omn_λ1_imag]
    )
    
end
    