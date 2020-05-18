include("basic_omnivory_module.jl")
include("top_heavy.jl")

function top_heavy_K_changes(K_vals, p)
    top_heavy_vals = Array{Float64}(undef, length(K_vals), 4)
    for (K_i, K_val) in enumerate(K_vals)
        u0 = [1.0, 0.5, 0.1]
        t_span = (0.0, 1000.0)
        par = p
        par.K = K_val
        prob = ODEProblem(model!, u0, t_span, par)
        sol = solve(prob)
        top_heavy_vals[K_i,:] = top_heavy(sol)
    end
    return top_heavy_vals
end

function top_heavy_plot(K_vals, p)
    top_heavy_vals = top_heavy_K_changes(K_vals, p)
    plot(K_vals, top_heavy_vals[:, 1], label = "Pred/Con+Res")
    plot(K_vals, top_heavy_vals[:, 2], label = "Pred/Res")
    plot(K_vals, top_heavy_vals[:, 3], label = "Pred/Con")
    plot(K_vals, top_heavy_vals[:, 4], label = "Res/Pred")
    xlabel("K")
    ylabel("Top Heavy Measure")
    ylim(0,12)
    legend()
    return
end

let
    chain_unforced_par = OmnPar(ω = 0.0, a_RP = 0.0, A = 0.0)
    chain_forced_par = OmnPar(ω = 0.0, a_RP = 0.0, A = 1.0)
    omn_unforced_par = OmnPar(ω = 0.5, a_RP = 0.8, A = 0.0)
    omn_forced_par = OmnPar(ω = 0.5, a_RP = 0.8, A = 1.0)
    fig_top_heavy = figure()
    subplot(2,2,1)
    title("Unforced", fontsize = 15)
    top_heavy_plot(K_vals = 1.5:0.1:3, chain_unforced_par)
    subplot(2,2,2)
    title("Forced", fontsize = 15)
    top_heavy_plot(K_vals = 1.5:0.1:3, chain_forced_par)
    subplot(2,2,3)
    top_heavy_plot(K_vals = 1.5:0.1:3, omn_unforced_par)
    subplot(2,2,4)
    top_heavy_plot(K_vals = 1.5:0.1:3, omn_forced_par)
    tight_layout(w_pad = 1.3)
    annotate("Food Chain", (470, 220), xycoords = "figure points", fontsize = 15, rotation = 90)
    annotate("Omnivory", (470, 75), xycoords = "figure points", fontsize = 15, rotation = 90)
    return fig_top_heavy
end


#NOTE when K < 1.5 returnng with error - Warning: dt <= dtmin. Aborting. There is either an error in your model specification or the true solution is unstable.
