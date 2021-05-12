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
pulse_strength = 2.0

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



let

    # NB: `pulse()` returns first eigen value, degree of overshoot and 
    # Max-Min it uses the same setup as `fig_press.jl`
    
    # range of K
    println("Simulations for a range of K values")
    res_K = []
    rg_K = 1:0.05:4
    for k in rg_K
        res_K = push!(res_K, pulse(k, 0.5, 0.6, 0.2))
    end 
    
    # range of a_CP
    println("Simulations for a range of aCP values")
    res_aCP = []
    rg_aCP = 0.25:.0025:.65
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
    # colors
    col_fc = "#000000"
    col_p = "#555555"
    col_r = "#cccccc"
    
    
    ##----------- K
    subplot(4, 4, 1)
    # fig = figure(figsize = (4, 4))
    plot(rg_K, [-1/res_K[i][1][1] for (i, K) in enumerate(rg_K)], label = "FC", color = col_fc)
    plot(rg_K, [-1/res_K[i][1][2] for (i, K) in enumerate(rg_K)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_K, [-1/res_K[i][1][3] for (i, K) in enumerate(rg_K)], label = "R", color = col_r)
    # legend()
    xlabel("K")
    ylabel("Re(λ1)")
    # tight_layout()
    # savefig("figs/fig_RT_K.svg")
    
    subplot(4, 4, 2)
    # fig = figure(figsize = (4, 4))
    plot(rg_K, [res_K[i][9][1] for (i, K) in enumerate(rg_K)], label = "FC", color = col_fc)
    plot(rg_K, [res_K[i][9][2] for (i, K) in enumerate(rg_K)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_K, [res_K[i][9][3] for (i, K) in enumerate(rg_K)], label = "R", color = col_r)
    # legend()
    xlabel("K")
    ylabel("Im(λ1)")
    
    # fig = figure(figsize = (4, 4))
    subplot(4, 4, 3)
    plot(rg_K, [res_K[i][5][1] for (i, K) in enumerate(rg_K)], label = "FC", color = col_fc)
    plot(rg_K, [res_K[i][5][2] for (i, K) in enumerate(rg_K)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_K, [res_K[i][5][3] for (i, K) in enumerate(rg_K)], label = "R", color = col_r)
    # legend()
    xlabel("K")
    ylabel("Degree of Overshoot")
    # tight_layout()
    # savefig("figs/fig_OS_K.svg")
    
    # fig = figure(figsize = (4, 4))
    subplot(4, 4, 4)
    plot(rg_K, [res_K[i][8][1] for (i, K) in enumerate(rg_K)], label = "FC", color = col_fc)
    plot(rg_K, [res_K[i][8][2] for (i, K) in enumerate(rg_K)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_K, [res_K[i][8][3] for (i, K) in enumerate(rg_K)], label = "R", color = col_r)
    legend()
    xlabel("K")
    ylabel("Max-Min")
    # tight_layout()
    # savefig("figs/fig_MM_K.svg")
    
        
        
        
    ##----------- aCP
    subplot(4, 4, 5)
    # fig = figure(figsize = (4, 4))
    plot(rg_aCP, [res_aCP[i][1][1] for (i, a) in enumerate(rg_aCP)], label = "FC", color = col_fc)
    plot(rg_aCP, [res_aCP[i][1][2] for (i, a) in enumerate(rg_aCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_aCP, [res_aCP[i][1][3] for (i, a) in enumerate(rg_aCP)], label = "R", color = col_r)
    # legend()
    xlabel("aCP")
    ylabel("Re(λ1)")
    # tight_layout()
    # savefig("figs/fig_RT_aCP.svg")
    
    subplot(4, 4, 6)
    # fig = figure(figsize = (4, 4))
    plot(rg_aCP, [res_aCP[i][9][1] for (i, K) in enumerate(rg_aCP)], label = "FC", color = col_fc)
    plot(rg_aCP, [res_aCP[i][9][2] for (i, K) in enumerate(rg_aCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_aCP, [res_aCP[i][9][3] for (i, K) in enumerate(rg_aCP)], label = "R", color = col_r)
    # legend()
    xlabel("aCP")
    ylabel("Im(λ1)")
    
    subplot(4, 4, 7)
    # fig = figure(figsize = (4, 4))
    plot(rg_aCP, [res_aCP[i][5][1] for (i, a) in enumerate(rg_aCP)], label = "FC", color = col_fc)
    plot(rg_aCP, [res_aCP[i][5][2] for (i, a) in enumerate(rg_aCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_aCP, [res_aCP[i][5][3] for (i, a) in enumerate(rg_aCP)], label = "R", color = col_r)
    # legend()
    xlabel("aCP")
    ylabel("Degree of Overshoot")
    # tight_layout()
    # savefig("figs/fig_OS_aCP.svg")
    
    subplot(4, 4, 8)
    # fig = figure(figsize = (4, 4))
    plot(rg_aCP, [res_aCP[i][8][1] for (i, a) in enumerate(rg_aCP)], label = "FC", color = col_fc)
    plot(rg_aCP, [res_aCP[i][8][2] for (i, a) in enumerate(rg_aCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_aCP, [res_aCP[i][8][3] for (i, a) in enumerate(rg_aCP)], label = "R", color = col_r)
    legend()
    xlabel("aCP")
    ylabel("Max-Min")
    # tight_layout()
    # savefig("figs/fig_MM_aCP.svg")
        
        
        
    ##----------- eCP
    subplot(4, 4, 9)
    # fig = figure(figsize = (4, 4))
    plot(rg_eCP, [res_eCP[i][1][1] for (i, K) in enumerate(rg_eCP)], label = "FC", color = col_fc)
    plot(rg_eCP, [res_eCP[i][1][2] for (i, K) in enumerate(rg_eCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_eCP, [res_eCP[i][1][3] for (i, K) in enumerate(rg_eCP)], label = "R", color = col_r)
    # legend()
    xlabel("eCP")
    ylabel("Re(λ1)")
    # tight_layout()
    # savefig("figs/fig_RT_eCP.svg")
    
    
    subplot(4, 4, 10)
    # fig = figure(figsize = (4, 4))
    plot(rg_eCP, [res_eCP[i][9][1] for (i, K) in enumerate(rg_eCP)], label = "FC", color = col_fc)
    plot(rg_eCP, [res_eCP[i][9][2] for (i, K) in enumerate(rg_eCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_eCP, [res_eCP[i][9][3] for (i, K) in enumerate(rg_eCP)], label = "R", color = col_r)
    # legend()
    xlabel("eCP")
    ylabel("Im(λ1)")
    
    subplot(4, 4, 11)
    # fig = figure(figsize = (4, 4))
    plot(rg_eCP, [res_eCP[i][5][1] for (i, K) in enumerate(rg_eCP)], label = "FC", color = col_fc)
    plot(rg_eCP, [res_eCP[i][5][2] for (i, K) in enumerate(rg_eCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_eCP, [res_eCP[i][5][3] for (i, K) in enumerate(rg_eCP)], label = "R", color = col_r)
    # legend()
    xlabel("eCP")
    ylabel("Degree of Overshoot")
    # tight_layout()
    # savefig("figs/fig_OS_eCP.svg")
    
    subplot(4, 4, 12)
    # fig = figure(figsize = (4, 4))
    plot(rg_eCP, [res_eCP[i][8][1] for (i, K) in enumerate(rg_eCP)], label = "FC",color = col_fc)
    plot(rg_eCP, [res_eCP[i][8][2] for (i, K) in enumerate(rg_eCP)], label = "P",color = col_p, linestyle = "dashed")
    plot(rg_eCP, [res_eCP[i][8][3] for (i, K) in enumerate(rg_eCP)], label = "R",color = col_r)
    legend()
    xlabel("eCP")
    ylabel("Max-Min")
    # tight_layout()
    # savefig("figs/fig_MM_eCP.svg")
    
    
    ##----------- mP
    subplot(4, 4, 13)
    # fig = figure(figsize = (4, 4))
    plot(rg_mP, [res_mP[i][1][1] for (i, K) in enumerate(rg_mP)], label = "FC", color = col_fc)
    plot(rg_mP, [res_mP[i][1][2] for (i, K) in enumerate(rg_mP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_mP, [res_mP[i][1][3] for (i, K) in enumerate(rg_mP)], label = "R", color = col_r)
    # legend()
    xlabel("mP")
    ylabel("Re(λ1)")
    # tight_layout()
    # savefig("figs/fig_RT_mP.svg")
    
    subplot(4, 4, 14)
    # fig = figure(figsize = (4, 4))
    plot(rg_mP, [res_mP[i][9][1] for (i, K) in enumerate(rg_mP)], label = "FC", color = col_fc)
    plot(rg_mP, [res_mP[i][9][2] for (i, K) in enumerate(rg_mP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_mP, [res_mP[i][9][3] for (i, K) in enumerate(rg_mP)], label = "R", color = col_r)
    # legend()
    xlabel("mP")
    ylabel("Im(λ1)")
    
    subplot(4, 4, 15)
    # fig = figure(figsize = (4, 4))
    plot(rg_mP, [res_mP[i][5][1] for (i, K) in enumerate(rg_mP)], label = "FC", color = col_fc)
    plot(rg_mP, [res_mP[i][5][2] for (i, K) in enumerate(rg_mP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_mP, [res_mP[i][5][3] for (i, K) in enumerate(rg_mP)], label = "R", color = col_r)
    # legend()
    xlabel("mP")
    ylabel("Degree of Overshoot")
    # tight_layout()
    # savefig("figs/fig_OS_mP.svg")
    
    #
    subplot(4, 4, 16)
    # fig = figure(figsize = (4, 4))
    plot(rg_mP, [res_mP[i][8][1] for (i, K) in enumerate(rg_mP)], label = "FC", color = col_fc)
    plot(rg_mP, [res_mP[i][8][2] for (i, K) in enumerate(rg_mP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_mP, [res_mP[i][8][3] for (i, K) in enumerate(rg_mP)], label = "R", color = col_r)
    legend()
    xlabel("mP")
    ylabel("Max-Min")
    tight_layout()
    savefig("figs/figS1_plus.svg")
    # savefig("figs/fig_MM_mP.svg")


    
end    
    
    # Complemetary figures (no longer used)
    # ## Imag
    # fig = figure(figsize = (4, 4))
    # plot(rg_K, [res_K[i][9][1] for (i, K) in enumerate(rg_K)], label = "FC", color = col_fc)
    # plot(rg_K, [res_K[i][9][2] for (i, K) in enumerate(rg_K)], label = "P", color = col_p, linestyle = "dashed")
    # plot(rg_K, [res_K[i][9][3] for (i, K) in enumerate(rg_K)], label = "R", color = col_r)
    # legend()
    # xlabel("K")
    # ylabel("Im(λ1)")
    # tight_layout()
    # savefig("figs/fig_im_K.svg")
    # 
    # 
    # fig = figure(figsize = (4, 4))
    # plot(rg_aCP, [res_aCP[i][9][1] for (i, a) in enumerate(rg_aCP)], label = "FC", color = col_fc)
    # plot(rg_aCP, [res_aCP[i][9][2] for (i, a) in enumerate(rg_aCP)], label = "P", color = col_p, linestyle = "dashed")
    # plot(rg_aCP, [res_aCP[i][9][3] for (i, a) in enumerate(rg_aCP)], label = "R", color = col_r)
    # legend()
    # xlabel("aCP")
    # ylabel("Im(λ1)")
    # tight_layout()
    # savefig("figs/fig_im_aCP.svg")
    # 
    # # subplot(4, 3, 7)
    # fig = figure(figsize = (4, 4))
    # plot(rg_eCP, [res_eCP[i][9][1] for (i, K) in enumerate(rg_eCP)], label = "FC", color = col_fc)
    # plot(rg_eCP, [res_eCP[i][9][2] for (i, K) in enumerate(rg_eCP)], label = "P", color = col_p, linestyle = "dashed")
    # plot(rg_eCP, [res_eCP[i][9][3] for (i, K) in enumerate(rg_eCP)], label = "R", color = col_r)
    # legend()
    # xlabel("eCP")
    # ylabel("Im(λ1)")
    # tight_layout()
    # savefig("figs/fig_im_eCP.svg")
    # 
    # fig = figure(figsize = (4, 4))
    # plot(rg_mP, [res_mP[i][9][1] for (i, K) in enumerate(rg_mP)], label = "FC", color = col_fc)
    # plot(rg_mP, [res_mP[i][9][2] for (i, K) in enumerate(rg_mP)], label = "P", color = col_p, linestyle = "dashed")
    # plot(rg_mP, [res_mP[i][9][3] for (i, K) in enumerate(rg_mP)], label = "R", color = col_r)
    # legend()
    # xlabel("mP")
    # ylabel("Img()")
    # tight_layout()
    # savefig("figs/fig_im_mP.svg")
    