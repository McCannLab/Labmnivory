include("basic_omnivory_module.jl")
using PyPlot

let
    u0 = [1.0, 0.5, 0.1]
    t_span = (0.0, 100.0)

    chain_par = ModelPar(ω = 0.0, a_CP = 0.6, K = 2.0)
    prob = ODEProblem(model!, u0, t_span, chain_par)
    sol_chain = solve(prob, abstol = 1e-8, reltol = 1e-8)

    omn_par = ModelPar(ω = 0.3, a_CP = 0.6, K = 2.0, pref = adapt_pref)
    prob = ODEProblem(model!, u0, t_span, omn_par)
    sol_omn = solve(prob, abstol = 1e-8, reltol = 1e-8)

    fig = figure()
    subplot(2, 1, 1)
    plot(sol_chain.t, sol_chain.u)
    title("Chain")

    subplot(2, 1, 2)
    plot(sol_omn.t, sol_omn.u)
    title("Omnivory")

    tight_layout()

    return fig
end
