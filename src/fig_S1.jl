include("basic_omnivory_module.jl")
include("pulse.jl")
using DifferentialEquations, NLsolve, QuadGK, PyPlot
pygui(true)

# Resource Pulse/Masting Events
# We model the situation where we have a relatively low
# productivity environment, but we have period high productivity pulses,
# like a oak mastings etc.
# Objectives are that we can see a response of omnivory in the face of this
# variable productivity situation, and test that we are measuring omnivory
# in a way that is enlightening.
global pulse_length = 2.0
global pulse_start = 200
global pulse_end = pulse_start + pulse_length
global pulse_event_times = union(pulse_start, pulse_end)
global pulse_strength = 2.0


pulse_event(u, t, integrator) = t ∈ pulse_event_times

function forcing!(integrator)
    if integrator.t == pulse_start
        integrator.p.K = pulse_strength + integrator.p.K_base
    elseif integrator.t == pulse_end
        integrator.p.K = integrator.p.K_base
    end
    return
end

cb = DiscreteCallback(pulse_event, forcing!)

function find_times_hit_equil(data)
    eq = data[1, end], data[2, end], data[3, end]
    times = zeros(3)
    for animal in 1:3
        for i in eachindex(data)
            if isapprox(data[animal,i], eq[animal], atol = 0.001)
                times[animal] = data.t[i]
                break
            end
        end
    end
    return times
end

function plot_unit(rg, res, id, leg = false, xlb = "", ylb = "")
    cols = ["#000000" "#555555" "#cccccc"]
    labs = ["FC" "P" "R"]
    lty = ["solid" "dashed" "solid"]
    for j in 1:3
        plot(rg, [res[i][id][j] for i in eachindex(rg)], 
            label = labs[j], color = cols[j], linestyle = lty[j])
    end 
    if leg 
        legend()
    end 
    xlabel(xlb)
    ylabel(ylb)
end 


let

    # NB: `pulse()` returns :
    # 1: real part of first eigen value, and 
    # 2: imaginary part of first eigen value
    # 3: reactivity 
    # 4-6: degree of overshoot
    # 7-9: max-min
    
    # range of K
    println("Simulations for a range of K values")
    res_K = []
    rg_K = 1.1:0.05:4
    for k in rg_K
        res_K = push!(res_K, pulse(k, 0.5, 0.6, 0.2))
    end 
    
    # range of a_CP
    println("Simulations for a range of aCP values")
    res_aCP = []
    rg_aCP = 0.25:.0025:.55
    for a in rg_aCP
        res_aCP = push!(res_aCP, pulse(3.0, a, 0.6, 0.2))
    end 
    
    # range of e_CP
    println("Simulations for a range of eCP values")
    res_eCP = []
    rg_eCP = 0.3:.0025:0.7
    for e in rg_eCP
        res_eCP = push!(res_eCP, pulse(3.0, 0.5, e, 0.2))
    end 
    
    # range of m_P
    println("Simulations for a range of mP values")
    res_mP = []
    rg_mP = 0.15:.005:.4
    for m in rg_mP
        res_mP = push!(res_mP, pulse(3.0, 0.5, 0.6, m))
    end 
    
    
    
    println("Drawing Figure S1")
    
    fig = figure(figsize = (14, 9))
    ind = [1 2 6 9]
    
    ##----------- K
    tly = ["Re(λ1)" "Im(λ1)" "Degree of Overshoot" "Max-Min"]
    tlx = ["" "" "" "K"]
    for i in 1:4
        subplot(4, 4, (i - 1) * 4 + 1)
        plot_unit(rg_K, res_K, ind[i], false, tlx[i], tly[i])
    end 
    
    ##----------- aCP
    tlx = ["" "" "" "aCP"]
    for i in 1:4
        subplot(4, 4, (i-1)*4 + 2)
        plot_unit(rg_aCP, res_aCP, ind[i], false, tlx[i], "")
    end 

    ##----------- eCP
    tlx = ["" "" "" "eCP"]
    for i in 1:4
        subplot(4, 4, (i-1)*4 + 3)
        plot_unit(rg_eCP, res_eCP, ind[i], false, tlx[i], "")
    end 

    ##----------- mP
    tlx = ["" "" "" "mP"]
    for i in 1:4
        subplot(4, 4, (i-1)*4 + 4)
        plot_unit(rg_mP, res_mP, ind[i], true, tlx[i], "")
    end 

    tight_layout()
    savefig("figs/figS1_plus.svg")

    
end    
    