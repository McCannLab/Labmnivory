include("basic_omnivory_module.jl")
include("press.jl")
using DifferentialEquations, NLsolve, QuadGK, PyPlot
pygui(true)

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
    rg_K = 2:0.05:4
    for k in rg_K
        res_K = push!(res_K, press(k, 0.5, 0.6, 0.2))
    end 
    
    # range of a_CP
    println("Simulations for a range of aCP values")
    res_aCP = []
    rg_aCP = 0.3:.0025:.55
    for a in rg_aCP
        res_aCP = push!(res_aCP, press(3.0, a, 0.6, 0.2))
    end 
    
    # range of e_CP
    println("Simulations for a range of eCP values")
    res_eCP = []
    rg_eCP = 0.35:.0025:0.7
    for e in rg_eCP
        res_eCP = push!(res_eCP, press(3.0, 0.5, e, 0.2))
    end 
    
    # range of m_P
    println("Simulations for a range of mP values")
    res_mP = []
    rg_mP = 0.15:.005:.3
    for m in rg_mP
        res_mP = push!(res_mP, press(3.0, 0.5, 0.6, m))
    end 
    
    
    
    println("Drawing Figure S2")
    
    fig = figure(figsize = (11, 6))
    ind = [2 7 10]
    
    ##----------- K
    tly = ["-1/Re(λ1)" "Degree of Overshoot" "Max-Min"]
    tlx = ["" "" "K"]
    for i in 1:3
        subplot(3, 4, (i - 1) * 4 + 1)
        plot_unit(rg_K, res_K, ind[i], false, tlx[i], tly[i])
    end 
    
    ##----------- aCP
    tlx = ["" "" "aCP"]
    for i in 1:3
        subplot(3, 4, (i-1) * 4 + 2)
        plot_unit(rg_aCP, res_aCP, ind[i], false, tlx[i], "")
    end 

    ##----------- eCP
    tlx = ["" "" "eCP"]
    for i in 1:3
        subplot(3, 4, (i-1) * 4 + 3)
        plot_unit(rg_eCP, res_eCP, ind[i], false, tlx[i], "")
    end 

    ##----------- mP
    tlx = ["" "" "mP"]
    for i in 1:3
        subplot(3, 4, (i-1) * 4 + 4)
        plot_unit(rg_mP, res_mP, ind[i], true, tlx[i], "")
    end 

    tight_layout()
    savefig("figs/figS2.svg")

    
end    
    