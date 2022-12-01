include("basic_omnivory_module.jl")
include("press.jl")
using Hwloc
using .Threads
using PyPlot
pygui(true)

# PRESS SIMULATIONS for ranges of K, aCP, eCP and mP.

# range of K
println("Simulations for a range of K values")
rg_K = 2:0.5:4
# rg_K = 2:0.05:4
res_K_press = Vector{Any}(undef,length(rg_K))
@threads for i in eachindex(rg_K)
    par = ModelPar()
    par.K_base = par.K = rg_K[i]
    res_K_press[i] = press(par, 0.1, 0.5)
end 

# range of a_CP
println("Simulations for a range of aCP values")
rg_aCP = 0.3:0.025:.55
# rg_aCP = 0.3:0.0025:.55
res_aCP_press = Vector{Any}(undef,length(rg_aCP))
@threads for i in eachindex(rg_aCP)
    par = ModelPar()
    par.a_CP = rg_aCP[i]
    res_aCP_press[i] = press(par, 0.1, 0.5)
end 

# range of e_CP
println("Simulations for a range of eCP values")
rg_eCP = 0.35:0.05:0.7
# rg_eCP = 0.35:.025:0.7
res_eCP_press = Vector{Any}(undef,length(rg_eCP))
@threads for i in eachindex(rg_eCP)
    par = ModelPar()
    par.e_CP = rg_eCP[i]
    res_eCP_press[i] = press(par, 0.1, 0.5)
end 

# range of m_P
println("Simulations for a range of mP values")
rg_mP = 0.15:0.05:0.3
# rg_mP = 0.15:0.005:0.3
res_mP_press = Vector{Any}(undef,length(rg_mP))
@threads for i in eachindex(rg_mP)
    par = ModelPar()
    par.m_P = rg_mP[i]
    res_mP_press[i] = press(par, 0.1, 0.5)
end 

# As this takes some time to run you may want to save this (it takes ~500Mo), 
# if so uncomment the lines below.
# using FileIO, JLD2
# # mkdir("res")
# save(
#     "res/res_S2.jld2", 
#     "rg_K", rg_K, "res_K_press", res_K_press,
#     "rg_aCP", rg_aCP, "res_aCP_press", res_aCP_press,
#     "rg_eCP", rg_eCP, "res_eCP_press", res_eCP,_press
#     "rg_mP", rg_mP, "res_mP_press", res_mP_press
#     )


# FIGURE
println("Drawing Figure S2")

fig = figure(figsize = (12, 8))
ind = [2 5 8]

##----------- K
tly = ["-1/Re(λ1)" "Degree of Overshoot" "Max-Min"]
tlx = ["" "" "K"]
for i in 1:3
    subplot(3, 4, (i - 1) * 4 + 1)
    plot_sa_unit(rg_K, res_K_press, ind[i], false, tlx[i], tly[i])
end 

##----------- aCP
tlx = ["" "" L"a_{CP}"]
for i in 1:3
    subplot(3, 4, (i-1) * 4 + 2)
    plot_sa_unit(rg_aCP, res_aCP_press, ind[i], false, tlx[i], "")
end 

##----------- eCP
tlx = ["" "" L"e_{CP}"]
for i in 1:3
    subplot(3, 4, (i-1) * 4 + 3)
    plot_sa_unit(rg_eCP, res_eCP_press, ind[i], false, tlx[i], "")
end 

##----------- mP
tlx = ["" "" L"m_P"]
for i in 1:3
    subplot(3, 4, (i-1) * 4 + 4)
    plot_sa_unit_rev(rg_mP, res_mP_press, ind[i], i == 1, tlx[i], "")
    val = [0.15, 0.20, 0.25, 0.30]
    xticks(val, reverse(val))
end 

tight_layout()
savefig("fig/fig_S2_parallel_reduced.svg")
# savefig("fig/fig_S2_parallel_all.svg")

