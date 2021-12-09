include("basic_omnivory_module.jl")
include("press.jl")
using PyPlot

# PRESS SIMULATION
par = ModelPar()
# NB: for FigS4, use
# par = ModelPar(e_CP = 0.5, e_RC = 0.5, e_RP = 0.5)
res = press(par, 0.1, 1.2)
# check equilibria 
rgs = [
    280:0.01:299,
    300:0.01:320,
    320:0.01:550,
    880:0.01:900
    ]
check_doo(res, rgs)

# plot helper
function plot_illustration(sol, eq, t_pr, ttl_id, leg = true, y_max = 5)
    RCP_cols = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    labs = ["R", "C", "P"]
    ttls = ["Food Chain", "Passive Omnivory", "Active Omnivory"]

    for i in 1:3
        hlines(eq[i], t_pr, sol.t[end], color = RCP_cols[i], alpha = 0.5)
    end 
    
    for i in 1:3
        plot(sol.t, sol[i, :], color = RCP_cols[i], label = labs[i])
    end 
    
    if leg
        legend()
    end 
    
    xlim(sol.t[1], sol.t[end],)
    ylim(0, y_max)
    title(ttls[ttl_id])
    ylabel("Density")
end 
# Layout
fig = figure(figsize = (8, 9))
st_t = 220
en_t = 900
t_pr = 300
rg_t = st_t:0.1:en_t
# Dynamics
## FC
subplot(3, 2, 1)
plot_illustration(
    res[1][14](rg_t),
    res[1][12],
    t_pr,
    1
)
## PO
subplot(3, 2, 3)
plot_illustration(
    res[2][14](rg_t),
    res[2][12],
    t_pr,
    2
)
## RO
subplot(3, 2, 5)
plot_illustration(
    res[3][14](rg_t),
    res[3][12],
    t_pr,
    3
)
xlabel("Time")
# Barplots (second Column)
## Stability 
subplot(3, 2, 2)
x_loc = [0, 1, 2]
x_labels = ["FC", "PO", "RO"]
plt.bar(x_loc, [res[i][2] for i in 1:3])
plt.xticks(x_loc, x_labels)
plt.tick_params(axis = "x", which = "both", length = 0)
ylabel("Local Return Time")
## Overshoot 
subplot(3, 2, 4)
### Set position of bar on X axis
### set width of bar
bar_width = 0.25
r1 = 1:3
r2 = [x + bar_width for x in r1]
r3 = [x + bar_width for x in r2]
### Make the plot
plt.bar(r1, [res[i][7] for i in 1:3], width = bar_width, edgecolor = "white", label = "R")
plt.bar(r2, [res[i][6] for i in 1:3], width = bar_width, edgecolor = "white", label = "C")
plt.bar(r3, [res[i][5] for i in 1:3], width = bar_width, edgecolor = "white", label = "P")
plt.xticks(r2, x_labels)
plt.tick_params(axis = "x", which = "both", length = 0)
ylabel("Degree of Overshoot")
plt.legend()
## Range (max-min)
subplot(3, 2, 6)
plt.bar(r1, [res[i][10] for i in 1:3], width = bar_width, edgecolor = "white", label = "R")
plt.bar(r2,  [res[i][9] for i in 1:3], width = bar_width, edgecolor = "white", label = "C")
plt.bar(r3,  [res[i][8] for i in 1:3], width = bar_width, edgecolor = "white", label = "P")
plt.xticks(r2, x_labels)
plt.tick_params(axis = "x", which = "both", length = 0)
ylabel("Max - Min")

plt.legend()

tight_layout()
savefig("fig/fig_3.svg")
# savefig("fig/fig_S3.svg")
