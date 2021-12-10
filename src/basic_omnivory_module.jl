using Parameters: @with_kw, @unpack
using LinearAlgebra: eigvals
using ForwardDiff
using QuadGK: quadgk

# Omnivory preference functions
function adapt_pref(u, p, t)
    return p.ω * u[1] / (p.ω * u[1] + (1 - p.ω) * u[2])
end

function fixed_pref(u, p, t)
    return p.Ω
end

# Parameters need to be defined after above functions so they can have defaults
@with_kw mutable struct ModelPar
    # Logistic Parameters
    r = 2.0
    # `K_base` measures the underyling K outside of any forcing applied
    K_base = 3
    K = 3
    # Consumer Parameters
    a_RC = 1
    h_RC = 0.4
    e_RC = 0.8
    m_C = 0.4
    # Predator Parameters
    a_CP = 0.5
    h_CP = 0.3
    e_CP = 0.6
    m_P = 0.2
    # Omnivory Parameters
    a_RP = 0.2
    h_RP = 0.6
    e_RP = 0.4
    # if Ω = 0 then we have a food chain
    Ω = 0.1
    # Forcing Function
    # 'fixed_pref' for passive omnivory
    # 'adapt_pref' for active omnivory
    pref::Function = fixed_pref
    # Used in the adaptive forcing to bias towards C or R
    ω = 0.5
end



#NOTE: we could use these in the model, but I am scared of all the function call
#      overhead, likely could be fixed with inlining, but will just leave it for now
function f_RP(u, p, t)
    R, C, P = u
    Ω = p.pref(u, p, t)
    return Ω * p.a_RP * R * P / (1 + Ω * p.a_RP * p.h_RP * R + (1 - Ω) * p.a_CP * p.h_CP * C)
end

function f_CP(u, p, t)
    R, C, P = u
    Ω = p.pref(u, p, t)
    return (1 - Ω) * p.a_CP * C * P / (1 + Ω * p.a_RP * p.h_RP * R + (1 - Ω) * p.a_CP * p.h_CP * C)
end

degree_omnivory(u, p) = f_RP(u, p, 0.0) / (f_RP(u, p, 0.0) + f_CP(u, p, 0.0))



function model!(du, u, p, t)
    @unpack r, K = p
    @unpack a_RC, h_RC, e_RC, m_C = p
    @unpack a_CP, h_CP, e_CP, m_P = p
    @unpack a_RP, h_RP, e_RP = p
    R, C, P = u

    # setup the density dependent preference
    Ω = p.pref(u, p, t)

    # compute those once for all
    int_RC = a_RC * R * C / (1 + a_RC * h_RC * R)
    denom_RCP = 1 + Ω * a_RP * h_RP * R + (1 - Ω) * a_CP * h_CP * C
    num_CP = (1 - Ω) * a_CP * C * P
    num_RP = Ω * a_RP * R * P
    
    # ODE
    du[1] = r * R * (1 - R / K) - int_RC - num_RP / denom_RCP
    du[2] = e_RC * int_RC - num_CP / denom_RCP - m_C * C
    du[3] = (e_RP * num_RP + e_CP * num_CP) / denom_RCP - m_P * P

    return du
end


# Utilities for doing eigenvalue analysis
function rhs(u, p)
    du = similar(u)
    model!(du, u, p, zero(u))
    return du
end

find_eq(u, p) = nlsolve((du, u) -> model!(du, u, p, zero(u)), u).zero
cmat(u, p) = ForwardDiff.jacobian(x -> rhs(x, p), u)

"""M is the community matrix, we can be calculated with `cmat(u, p)`"""
λ1_stability(M) = maximum(real.(eigvals(M)))

"""M is the community matrix, we can be calculated with `cmat(u, p)`"""
function λ1_stability_imag(M) 
    ev = eigvals(M)
    imag.(ev)[findall(real.(ev) .== maximum(real.(ev)))[1]]
end 
    

"""M is the community matrix, we can be calculated with `cmat(u, p)`
Note: `\nu` is the what to input `ν` which looks a bit too much like `v` for my taste
"""
ν_stability(M) = λ1_stability((M + M') / 2)


# overshoot and oscillation range

abs_sol(sol, t, eq) = abs.(sol(t) .- eq)

function overshoot(sol, eq, spc, t_beg, t_end)
    return quadgk(t -> abs_sol(sol, t, eq)[spc], t_beg, t_end)[1]
end


# Calculate max-min metric
function min_max(sol, spc, t_beg, t_end, len = 100000)
    return maximum(sol(range(t_beg, t_end, length = len))[spc, :]) - 
    minimum(sol(range(t_beg, t_end, length = len))[spc, :])
end

# find first time equilibrium is hit
function find_times_hit_equil_press(res)
    eq = res[1, end], res[2, end], res[3, end]
    times = zeros(3)
    for spc in 1:3
        for i in 20:length(res)
            # cannot be too strict here otherwise the value of the 
            # first ht time varies a lort which will have serious 
            # impact on min and max (overshoot too) leading to major oscillations lenth of the ts must be high enough too.
            if isapprox(res[spc, i], eq[spc], atol = 0.01)
                times[spc] = res.t[i]
                break
            end
        end
    end
    return times
end

# Basic plot for Sensitivity analysis (well akin to sa)
function plot_sa_unit(rg, res, id, leg = false, xlb = "", ylb = "", ann = 1)
    # az = 'a':'z'
    cols = ["#000000" "#555555" "#cccccc"]
    labs = ["FC" "PO" "RO"]
    lty = ["solid" "dashed" "solid"]
    for j in 1:3
        plot(rg, [res[i][j][id] for i in eachindex(rg)], 
            label = labs[j], color = cols[j], linestyle = lty[j])
    end 
    if leg 
        legend()
    end 
    xlabel(xlb)
    ylabel(ylb)
    # did not work
    # annotate(string(az[ann]), [.1;.1], xycoords="figure fraction", 
    #     annotation_clip = false)
end 

# Basic plot for Sensitivity analysis (well akin to sa)
function plot_sa_unit_rev(rg, res, id, leg = false, xlb = "", ylb = "")
    # az = 'a':'z'
    cols = ["#000000" "#555555" "#cccccc"]
    labs = ["FC" "PO" "RO"]
    lty = ["solid" "dashed" "solid"]
    for j in 1:3
        plot(rg, [res[i][j][id] for i in reverse(eachindex(rg))], label = labs[j], color = cols[j], linestyle = lty[j])
    end 
    if leg 
        legend()
    end 
    xlabel(xlb)
    ylabel(ylb)
end 


# Helper function to check out degree of omnivory
function check_doo(res, rgs)
    type = ["Fixed" "Responsive"]
    phase = ["Equilibrium", "Pulse", "Transient", "New Equilibrium"]
    # Check values of degree of omnivory
    for i in 2:3
        printstyled(type[i - 1], "---------\n", color = :green)
        for j in eachindex(rgs)
            printstyled(phase[j], " --> ", color = :blue)
            # see press_unit output
            sol = res[i][14](rgs[j])
            println(
                maximum(
                    [degree_omnivory(
                    sol[:, k], 
                    res[i][13]
                    ) for k in eachindex(sol)]
                )
            )
            # maximum degree of omnivory when the C is min (food web is top heavy)
            printstyled("where topheaviest --> ", color = :blue)
            idmx = argmin([sol[2, k] for k in eachindex(sol)])
            # println(idmx)
            print(
                degree_omnivory(
                    sol[idmx], 
                    res[i][13]
                )
            )
            ti = rgs[j][idmx]
            println("(time = $ti)")
        end
    end
end