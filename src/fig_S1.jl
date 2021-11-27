include("basic_omnivory_module.jl")
include("pulse.jl")
using PyPlot
pygui(true)



let

    # range of K
    println("Simulations for a range of K values")
    res_K = []
    rg_K = 2:0.05:4
    for k in rg_K
        par = ModelPar()
        par.K_base = par.K = k
        res_K = push!(res_K, pulse(par, 0.1, 2.0, 2.0))
    end 
    
    # range of a_CP
    println("Simulations for a range of aCP values")
    res_aCP = []
    rg_aCP = 0.3:.0025:.55
    for a in rg_aCP
        par = ModelPar()
        par.a_CP = a
        res_aCP = push!(res_aCP, pulse(par, 0.1, 2.0, 2.0))
    end 
    
    # range of e_CP
    println("Simulations for a range of eCP values")
    res_eCP = []
    rg_eCP = 0.35:.0025:0.7
    for e in rg_eCP
        par = ModelPar()
        par.e_CP = e
        res_eCP = push!(res_eCP, pulse(par, 0.1, 2.0, 2.0))
    end 
    
    # range of m_P
    println("Simulations for a range of mP values")
    res_mP = []
    rg_mP = 0.15:.005:.3
    for m in rg_mP
        par = ModelPar()
        par.m_P = m
        res_mP = push!(res_mP, pulse(par, 0.1, 2.0, 2.0))
    end 
    
    # save takes ~500Mo
    # using FileIO, JLD2
    # mkdir("res")
    # save(
    #     "res/res_S1.jld2", 
    #     "rg_K", rg_K, "res_K", res_K,
    #     "rg_aCP", rg_aCP, "res_aCP", res_aCP,
    #     "rg_eCP", rg_eCP, "res_eCP", res_eCP,
    #     "rg_mP", rg_mP, "res_mP", res_mP
    #     )
    
    
    println("Drawing Figure S1")
    
    fig = figure(figsize = (12, 8))
    ind = [2 5 8]
    
    ##----------- K
    tly = ["-1/Re(λ1)" "Degree of Overshoot" "Max-Min"]
    tlx = ["" "" "K"]
    for i in 1:3
        subplot(3, 4, (i - 1) * 4 + 1)
        plot_sa_unit(rg_K, res_K, ind[i], false, tlx[i], tly[i])
    end 
    
    ##----------- aCP
    tlx = ["" "" L"a_{CP}"]
    for i in 1:3
        subplot(3, 4, (i-1) * 4 + 2)
        plot_sa_unit(rg_aCP, res_aCP, ind[i], false, tlx[i], "")
    end 

    ##----------- eCP
    tlx = ["" "" L"e_{CP}"]
    for i in 1:3
        subplot(3, 4, (i-1) * 4 + 3)
        plot_sa_unit(rg_eCP, res_eCP, ind[i], false, tlx[i], "")
    end 

    ##----------- mP
    tlx = ["" "" L"m_P"]
    for i in 1:3
        subplot(3, 4, (i-1) * 4 + 4)
        plot_sa_unit_rev(rg_mP, res_mP, ind[i], i == 1, tlx[i], "")
        val = [0.15, 0.20, 0.25, 0.30]
        xticks(val, reverse(val))
    end 
    
    
    tight_layout()
    savefig("figs/figS1_v3.svg")

    
end    
    