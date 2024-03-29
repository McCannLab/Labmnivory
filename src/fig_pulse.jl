include("basic_omnivory_module.jl")
include("pulse.jl")
using PyPlot
pygui(true)

# PULSE SIMULATIONS
# returns dynamics and metrics in different panels, we actually dispatch the 
# different pannels between figure 3 and 4 with Inkscape (we did the same
# for figures S3 and S4).

# we pass "true" to this ARGS for supplementary figures
if isempty(ARGS) || ARGS[1] != "true"
    # default for main figures 3 and 4
    par = ModelPar()
    fl_fig = "fig/fig_pulse.svg"
else 
    # parameters for SI figures S3 and S4
    par = ModelPar(e_CP = 0.5, e_RC = 0.5, e_RP = 0.5)
    fl_fig = "fig/fig_pulse_si.svg"
end

# run pulse simulation
res = pulse(par, 0.1, 2.0, 2.0)

# temporal ranges of the different phases
rgs = [
    180:0.01:199,
    200:0.01:205,
    205:0.01:325,
    325:0.01:350
    ]
# check degree of omnivory for the different phases
check_doo(res, rgs)

# plot helper
function plot_illustration(sol, eq, ttl_id, leg = true, y_max = 5)
    RCP_cols = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    labs = ["R", "C", "P"]
    ttls = ["Food Chain", "Passive Omnivory", "Active Omnivory"]

    for i in 1:3
        hlines(eq[i], sol.t[1], sol.t[end], color = RCP_cols[i], alpha = 0.5)
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
st_t = 175
en_t = 350
rg_t = st_t:0.1:en_t
# Dynamics
## FC
subplot(3, 2, 1)
plot_illustration(
    res[1][14](rg_t),
    res[1][12],
    1
)
## PO
subplot(3, 2, 3)
plot_illustration(
    res[2][14](rg_t),
    res[2][12],
    2
)
## RO
subplot(3, 2, 5)
plot_illustration(
    res[3][14](rg_t),
    res[3][12],
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
savefig(fl_fig)

