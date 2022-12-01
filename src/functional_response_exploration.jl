include("basic_omnivory_module.jl")
using PyPlot

2/1
0.05/1

let
    ratio = 0.05:0.05:2.0
    con = 1
    par = ModelPar()
    omega = zeros(length(ratio))
    for i in eachindex(ratio)
        u = [ratio[i], 1]
        omega[i] = adapt_pref(u, par, 0)
    end
    test = figure()
    plot(ratio, omega, color = "blue")
    hlines(par.Ω, 0.00, 2.00, color = "red")
    xlabel("Resource:Consumer")
    ylabel("Ω (on resource)")
    # return test
    savefig("fig/omega_response.svg")
end
