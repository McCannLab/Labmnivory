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
global first_mast = 100.0
global mast_freq = 100.0
mast_length = 2.0
global t_end = mast_freq * 2
global t_span = (0.0, t_end)

global mast_starts = first_mast:mast_freq:t_end
global mast_ends = mast_starts .+ mast_length
global mast_event_times = sort(union(mast_starts, mast_ends))

mast_strength = 2.0

masting_event(u, t, integrator) = t ∈ mast_event_times

function forcing!(integrator)
    if integrator.t ∈ mast_starts
        integrator.p.K = mast_strength * integrator.p.K_base
    elseif integrator.t ∈ mast_ends
        integrator.p.K = integrator.p.K_base
    end
    return
end

cb = DiscreteCallback(masting_event, forcing!)




let

    # NB: `pulse()` returns first eigen value, degree of overshoot and 
    # Max-Min
    
    # range of K
    res_K = []
    rg_K = 2.25:0.01:3.75
    for k in rg_K
        res_K = push!(res_K, pulse(k, 0.25, 0.6, 0.2))
    end 
    
    # range of a_CP
    res_aCP = []
    rg_aCP = 0.201:.001:.35
    for a in rg_aCP
        res_aCP = push!(res_aCP, pulse(3.0, a, 0.6, 0.2))
    end 
    
    # range of e_CP
    res_eCP = []
    rg_eCP = 0.51:.0025:0.8
    for e in rg_eCP
        res_eCP = push!(res_eCP, pulse(3.0, 0.25, e, 0.2))
    end 
    
    # range of m_P
    res_mP = []
    rg_mP = 0.16:.0005:.23
    for m in rg_mP
        res_mP = push!(res_mP, pulse(3.0, 0.25, 0.6, m))
    end 
    
    # colors
    col_fc = "#000000"
    col_p = "#555555"
    col_r = "#cccccc"
    
    ##----------- K
    println("Range of K")
    fig = figure(figsize = (16, 12))
    subplot(4, 3, 1)
    # fig = figure(figsize = (4, 4))
    plot(rg_K, [-1/res_K[i][1][1] for (i, K) in enumerate(rg_K)], label = "FC", color = col_fc)
    plot(rg_K, [-1/res_K[i][1][2] for (i, K) in enumerate(rg_K)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_K, [-1/res_K[i][1][3] for (i, K) in enumerate(rg_K)], label = "R", color = col_r)
    legend()
    xlabel("K")
    ylabel("-1/Re(λ1)")
    # tight_layout()
    # savefig("figs/fig_RT_K.svg")
    
    # fig = figure(figsize = (4, 4))
    subplot(4, 3, 2)
    plot(rg_K, [res_K[i][3][3] for (i, K) in enumerate(rg_K)], label = "FC", color = col_fc)
    plot(rg_K, [res_K[i][4][3] for (i, K) in enumerate(rg_K)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_K, [res_K[i][5][3] for (i, K) in enumerate(rg_K)], label = "R", color = col_r)
    # legend()
    xlabel("K")
    ylabel("Degree of Overshoot")
    # tight_layout()
    # savefig("figs/fig_OS_K.svg")
    
    # fig = figure(figsize = (4, 4))
    subplot(4, 3, 3)
    plot(rg_K, [res_K[i][6][3] for (i, K) in enumerate(rg_K)], label = "FC", color = col_fc)
    plot(rg_K, [res_K[i][7][3] for (i, K) in enumerate(rg_K)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_K, [res_K[i][8][3] for (i, K) in enumerate(rg_K)], label = "R", color = col_r)
    # legend()
    xlabel("K")
    ylabel("Max-Min")
    # tight_layout()
    # savefig("figs/fig_MM_K.svg")
    
        
    ##----------- aCP
    println("Range of aCP")
    subplot(4, 3, 4)
    # fig = figure(figsize = (4, 4))
    plot(rg_aCP, [-1/res_aCP[i][1][1] for (i, a) in enumerate(rg_aCP)], label = "FC", color = col_fc)
    plot(rg_aCP, [-1/res_aCP[i][1][2] for (i, a) in enumerate(rg_aCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_aCP, [-1/res_aCP[i][1][3] for (i, a) in enumerate(rg_aCP)], label = "R", color = col_r)
    legend()
    xlabel("aCP")
    ylabel("-1/Re(λ1)")
    # tight_layout()
    # savefig("figs/fig_RT_aCP.svg")
    
    subplot(4, 3, 5)
    # fig = figure(figsize = (4, 4))
    plot(rg_aCP, [res_aCP[i][3][3] for (i, a) in enumerate(rg_aCP)], label = "FC", color = col_fc)
    plot(rg_aCP, [res_aCP[i][4][3] for (i, a) in enumerate(rg_aCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_aCP, [res_aCP[i][5][3] for (i, a) in enumerate(rg_aCP)], label = "R", color = col_r)
    # legend()
    xlabel("aCP")
    ylabel("Degree of Overshoot")
    # tight_layout()
    # savefig("figs/fig_OS_aCP.svg")
    
    subplot(4, 3, 6)
    # fig = figure(figsize = (4, 4))
    plot(rg_aCP, [res_aCP[i][6][3] for (i, a) in enumerate(rg_aCP)], label = "FC", color = col_fc)
    plot(rg_aCP, [res_aCP[i][7][3] for (i, a) in enumerate(rg_aCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_aCP, [res_aCP[i][8][3] for (i, a) in enumerate(rg_aCP)], label = "R", color = col_r)
    # legend()
    xlabel("aCP")
    ylabel("Max-Min")
    # tight_layout()
    # savefig("figs/fig_MM_aCP.svg")
        
        
    ##----------- eCP
    println("Range of eCP")
    subplot(4, 3, 7)
    # fig = figure(figsize = (4, 4))
    plot(rg_eCP, [-1/res_eCP[i][1][1] for (i, K) in enumerate(rg_eCP)], label = "FC", color = col_fc)
    plot(rg_eCP, [-1/res_eCP[i][1][2] for (i, K) in enumerate(rg_eCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_eCP, [-1/res_eCP[i][1][3] for (i, K) in enumerate(rg_eCP)], label = "R", color = col_r)
    legend()
    xlabel("eCP")
    ylabel("-1/Re(λ1)")
    # tight_layout()
    # savefig("figs/fig_RT_eCP.svg")
    
    subplot(4, 3, 8)
    # fig = figure(figsize = (4, 4))
    plot(rg_eCP, [res_eCP[i][3][3] for (i, K) in enumerate(rg_eCP)], label = "FC", color = col_fc)
    plot(rg_eCP, [res_eCP[i][4][3] for (i, K) in enumerate(rg_eCP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_eCP, [res_eCP[i][5][3] for (i, K) in enumerate(rg_eCP)], label = "R", color = col_r)
    # legend()
    xlabel("eCP")
    ylabel("Degree of Overshoot")
    # tight_layout()
    # savefig("figs/fig_OS_eCP.svg")
    
    subplot(4, 3, 9)
    # fig = figure(figsize = (4, 4))
    plot(rg_eCP, [res_eCP[i][6][3] for (i, K) in enumerate(rg_eCP)], label = "FC",color = col_fc)
    plot(rg_eCP, [res_eCP[i][7][3] for (i, K) in enumerate(rg_eCP)], label = "P",color = col_p, linestyle = "dashed")
    plot(rg_eCP, [res_eCP[i][8][3] for (i, K) in enumerate(rg_eCP)], label = "R",color = col_r)
    # legend()
    xlabel("eCP")
    ylabel("Max-Min")
    # tight_layout()
    # savefig("figs/fig_MM_eCP.svg")
    
    
    ##----------- mP
    println("Range of mP")
    subplot(4, 3, 10)
    # fig = figure(figsize = (4, 4))
    plot(rg_mP, [-1/res_mP[i][1][1] for (i, K) in enumerate(rg_mP)], label = "FC", color = col_fc)
    plot(rg_mP, [-1/res_mP[i][1][2] for (i, K) in enumerate(rg_mP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_mP, [-1/res_mP[i][1][3] for (i, K) in enumerate(rg_mP)], label = "R", color = col_r)
    legend()
    xlabel("mP")
    ylabel("-1/Re(λ1)")
    # tight_layout()
    # savefig("figs/fig_RT_mP.svg")
    
    subplot(4, 3, 11)
    # fig = figure(figsize = (4, 4))
    plot(rg_mP, [res_mP[i][3][3] for (i, K) in enumerate(rg_mP)], label = "FC", color = col_fc)
    plot(rg_mP, [res_mP[i][4][3] for (i, K) in enumerate(rg_mP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_mP, [res_mP[i][5][3] for (i, K) in enumerate(rg_mP)], label = "R", color = col_r)
    # legend()
    xlabel("mP")
    ylabel("Degree of Overshoot")
    # tight_layout()
    # savefig("figs/fig_OS_mP.svg")
    
    #
    subplot(4, 3, 12)
    # fig = figure(figsize = (4, 4))
    plot(rg_mP, [res_mP[i][6][3] for (i, K) in enumerate(rg_mP)], label = "FC", color = col_fc)
    plot(rg_mP, [res_mP[i][7][3] for (i, K) in enumerate(rg_mP)], label = "P", color = col_p, linestyle = "dashed")
    plot(rg_mP, [res_mP[i][8][3] for (i, K) in enumerate(rg_mP)], label = "R", color = col_r)
    # legend()
    xlabel("mP")
    ylabel("Max-Min")
    tight_layout()
    savefig("figs/fig4.svg")
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
    