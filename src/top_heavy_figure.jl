include("basic_omnivory_module.jl")
include("top_heavy.jl")

function top_heavy_K_changes(Ks, p)
    top_heavy_vals = fill(0.0, length(Ks), 4)
    for (i, K) in enumerate(Ks)
        u0 = [1.0, 0.5, 0.1]
        t_span = (0.0, 1000.0)
        par = deepcopy(p)
        par.K = K
        prob = ODEProblem(model!, u0, t_span, par)
        sol = solve(prob, reltol = 1e-8, abstol = 1e-8)
        top_heavy_vals[i, :] = top_heavy(sol)
    end

    return top_heavy_vals
end

function top_heavy_plot(Ks, p)
    top_heavy_vals = top_heavy_K_changes(Ks, p)

    plot(Ks, top_heavy_vals[:, 1], label = "P:(C + R)")
    plot(Ks, top_heavy_vals[:, 2], label = "P:R")
    plot(Ks, top_heavy_vals[:, 3], label = "P:C")
    plot(Ks, top_heavy_vals[:, 4], label = "R:P")
    xlabel("K")
    ylabel("Top Heavy Measure")
    legend()

    return
end

let
    # Setup up Model Configurations
    chain_unforced_par = ModelPar(ω = 0.0, a_RP = 0.0, A = 0.0)
    chain_forced_par = ModelPar(ω = 0.0, a_RP = 0.0, A = 1.0)
    omn_unforced_par = ModelPar(ω = 0.5, a_RP = 0.8, A = 0.0)
    omn_forced_par = ModelPar(ω = 0.5, a_RP = 0.8, A = 1.0)
    Ks = range(1.5, 3, length = 50)

    fig_top_heavy = figure(figsize = (12, 7))
    subplot(2, 2, 1)
    top_heavy_plot(Ks, chain_unforced_par)
    title("Food Chain Unforced")

    subplot(2, 2, 2)
    top_heavy_plot(Ks, chain_forced_par)
    title("Food Chain Forced")

    subplot(2, 2, 3)
    top_heavy_plot(Ks, omn_unforced_par)
    title("Omnivory Unforced")

    subplot(2, 2, 4)
    top_heavy_plot(Ks, omn_forced_par)
    title("Omnivory Forced")

    tight_layout(pad = 1.6)

    return fig_top_heavy
end
