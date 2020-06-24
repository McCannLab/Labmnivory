using NLsolve, PyPlot
using DataFrames, CSV
include("basic_omnivory_module.jl")
include("top_heavy.jl")
include("asynchrony.jl")
include("masting_event.jl")

grad_K = 2.2:.1:4
nK = length(grad_K)
grad_eRP = .02:.02:.5
neRP = length(grad_eRP)
u0 = [1.0, 1.5, 1.5]

# The global basic level of "Omnivory" we are looking at:
Ω = 0.1

# K LOOP
res_K = DataFrame(
    foodweb = repeat(["chain", "omn_fixed", "omn"], inner = nK),
    K = repeat(grad_K, 3),
    cv_P = 0.0,
    return_time = 0.0,
    asyn_RP = 0.0,
    auc_P = 0.0
    )

for (i, k) in enumerate(grad_K)
    # Parameters
    par_chain = ModelPar(a_CP = 0.25, Ω = 0.0, K_base = k, K = k)
    par_omn_fixed = ModelPar(a_CP = 0.25, Ω = Ω, pref = fixed_pref, K_base = k, K = k)
    pars = [par_chain, par_omn_fixed]
    sols = Any[]
    for j in 1:3
        if j == 3
            # need to use sols[2] for active Omnivory
            eq = nlsolve((du, u) -> model!(du, u, pars[2], 0.0), sols[2][end]).zero
            ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
            par = ModelPar(a_CP = 0.25, Ω = Ω, pref = adapt_pref, K_base = k, K = k, ω = ω)
        else par = pars[j]
        end
        prob_eq = ODEProblem(model!, u0, t_span, deepcopy(par))
        sol_eq = solve(prob_eq, reltol = 1e-8, abstol = 1e-8)
        prob = ODEProblem(model!, u0, t_span, deepcopy(par), tstops = mast_event_times)
        sol = solve(prob, reltol = 1e-8, abstol = 1e-8, callback = cb)
        push!(sols, sol)
        # Write results
        id = nK*(j-1) + i
        println(id)
        s_all = sol(t_all)
        res_K[id, :cv_P] = global_cv(s_all)[3]
        res_K[id, :asyn_RP] = global_asyn(s_all)[1]
        res_K[id, :auc_P] = global_auc(s_all, sol_eq(t_all), dt)[3]
        res_K[id, :return_time] = return_time(s_all, 3, dt, 1e-10)
    end
end
# Save results
# CSV.write("../res/res_K.csv", res_K)
# res_K = DataFrame(CSV.file("../res/res_K.csv"))





# eRP LOOP
res_eRP = DataFrame(
    foodweb = repeat(["chain", "omn_fixed", "omn"], inner = neRP),
    e_RP = repeat(grad_eRP, 3),
    cv_P = 0.0,
    return_time = 0.0,
    asyn_RP = 0.0,
    auc_P = 0.0
    )
# NB: resuls for chain are useless
for (i, e) in enumerate(grad_eRP)
    # Parameters
    par_chain = ModelPar(a_CP = 0.25, Ω = 0.0, e_RP = e)
    par_omn_fixed = ModelPar(a_CP = 0.25, Ω = Ω, pref = fixed_pref, e_RP = e)
    pars = [par_chain, par_omn_fixed]
    sols = Any[]
    for j in 1:3
        if j == 3
            eq = nlsolve((du, u) -> model!(du, u, pars[2], 0.0), sols[2][end]).zero
            ω = Ω * eq[2] / (Ω * eq[2] + (1 - Ω) * eq[1])
            par = ModelPar(a_CP = 0.25, Ω = Ω, pref = adapt_pref, e_RP = e, ω = ω)
        else par = pars[j]
        end
        prob_eq = ODEProblem(model!, u0, t_span, deepcopy(par))
        sol_eq = solve(prob_eq, reltol = 1e-8, abstol = 1e-8)
        prob = ODEProblem(model!, u0, t_span, deepcopy(par), tstops = mast_event_times)
        sol = solve(prob, reltol = 1e-8, abstol = 1e-8, callback = cb)
        push!(sols, sol)
        # Write results
        id = neRP*(j-1) + i
        s_all = sol(t_all)
        println(id)
        res_eRP[id, :cv_P] = global_cv(sol(t_aft))[3]
        res_eRP[id, :return_time] = return_time(s_all, 3, dt, 1e-10)
        res_eRP[id, :asyn_RP] = global_asyn(s_all)[1]
        res_eRP[id, :auc_P] = global_auc(s_all, sol_eq(t_all), dt)[3]
    end
end
# Save results
# CSV.write("../res/res_eRP.csv", res_eRP)





# Figure

tls = ["SD (P)", "Return time", "Asynchrony (C-R)", "AUC (P)"]
lss = ["--", "-", "-"]
cls = ["black", "black", "grey"]
lbs = ["Food Chain", "Omnivory [Passive]", "Omnivory [Active]"]

dfs = [
    res_K[res_K.foodweb .== "chain", :],
    res_K[res_K.foodweb .== "omn_fixed", :],
    res_K[res_K.foodweb .== "omn", :]
    ]
for (k, j) in enumerate([3 4 5 6])
    for i in 1:3
        subplot(4, 2, 2*k - 1)
        plot(dfs[i][:, :K], dfs[i][:, j], color = cls[i], linestyle = lss[i], label = lbs[i])
    end
    if k == 1
        legend()
    end
    if j ∈ [6]
        xlabel("Carrying capacity (K)")
    end
    title(tls[k])
end

dfs = [
    res_eRP[res_eRP.foodweb .== "chain", :],
    res_eRP[res_eRP.foodweb .== "omn_fixed", :],
    res_eRP[res_eRP.foodweb .== "omn", :]
    ]
for (k, j) in enumerate([3 4 5 6])
    for i in 1:3
        subplot(4, 2, 2*k)
        plot(dfs[i][:, :e_RP], dfs[i][:, j], color = cls[i], linestyle = lss[i])
    end
    if j ∈ [6]
        xlabel("Conversoin efficiency R-P")
    end
    title(tls[k])
end