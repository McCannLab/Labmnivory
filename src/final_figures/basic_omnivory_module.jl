using Parameters: @with_kw, @unpack

using LinearAlgebra: eigvals
using ForwardDiff

# # Forcing Functions
force(u, p, t) = p.A * sin(2 * π * t / p.B + p.𝛗 * π)

# # Omnivory preference functions
function adapt_pref(u, p, t)
    return p.ω * u[1] / (p.ω * u[1] + (1 - p.ω) * u[2])
end

function fixed_pref(u, p, t)
    return p.Ω
end

# # Parameters need to be defined after above functions so they can have defaults
@with_kw mutable struct ModelPar{F <: Function}
    # Logistic Parameters
    r = 2.0
    ## `K_base` measures the underyling K outside of any forcing applied
    K_base = 3.0
    K = 3.0
    # Consumer Parameters
    a_RC = 1.2
    h_RC = 0.8
    e_RC = 0.7
    m_C = 0.4
    # Predator Parameters
    a_CP = 0.8
    h_CP = 0.6
    e_CP = 0.6
    m_P = 0.2
    # Omnivory Parameters
    a_RP = 0.8
    h_RP = 0.9
    e_RP = 0.4
    Ω = 0.1
    # Forcing Function
    pref::F = fixed_pref
    ## Used in the adaptive forcing to bias towards C or R
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

    du[1] = r * R * (1 - R / K) - a_RC * R * C / (1 + a_RC * h_RC * R) - Ω * a_RP * R * P / (1 + Ω * a_RP * h_RP * R + (1 - Ω) * a_CP * h_CP * C)
    du[2] = e_RC * a_RC * R * C / (1 + a_RC * h_RC * R) - (1 - Ω) * a_CP * C * P / (1 + Ω * a_RP * h_RP * R + (1 - Ω) * a_CP * h_CP * C) - m_C * C
    du[3] = (e_RP * Ω * a_RP * R * P + e_CP * (1 - Ω) * a_CP * C * P) / (1 + Ω * a_RP * h_RP * R + (1 - Ω) * a_CP * h_CP * C) - m_P * P

    return du
end


# Utilities for doing eigenvalue analysis using autodiff
function rhs(u, p)
    du = similar(u)
    model!(du, u, p, zero(u))
    return du
end

find_eq(u, p) = nlsolve((du, u) -> model!(du, u, p, zero(u)), u).zero
cmat(u, p) = ForwardDiff.jacobian(x -> rhs(x, p), u)

"""M is the community matrix, we can be calculated with `cmat(u, p)`"""
λ1_stability(M) = maximum(real.(eigvals(M)))

"""M is the community matrix, we can be calculated with `cmat(u, p)`
Note: `\nu` is the what to input `ν` which looks a bit to much like `v` for my taste
"""
ν_stability(M) = λ1_stability((M + M') / 2)
