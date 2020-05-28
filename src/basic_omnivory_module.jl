using Parameters
using DifferentialEquations

# # Forcing Functions
force(u, p, t) = p.A * sin(2 * π * t / p.B + p.𝛗 * π)

# # Omnivory preference functions
adapt_pref(u, p, t) = p.ω * u[1] / (p.ω * u[1] + (1 - p.ω) * u[2])
fixed_pref(u, p, t) = p.ω

# # Parameters need to be defined after above functions so they can have defaults
@with_kw mutable struct ModelPar{F <: Function}
    # Logistic Parameters
    r = 2.0
    ## `K_base` measings the underyling K outside of any forcing applied
    K_base = 3.0
    K = 3.0
    # Consumer Parameters
    a_RC = 1.1
    h_RC = 0.8
    e_RC = 0.7
    m_C = 0.4
    # Predator Parameters
    a_CP = 1.1
    h_CP = 0.6
    e_CP = 0.6
    m_P = 0.2
    # Omnivory Parameters
    a_RP = 0.8
    h_RP = 0.9
    e_RP = 0.4
    ω = 0.5
    # Forcing Function
    pref::F = fixed_pref
end


function model!(du, u, p, t)
    @unpack r, K = p
    @unpack a_RC, h_RC, e_RC, m_C = p
    @unpack a_CP, h_CP, e_CP, m_P = p
    @unpack a_RP, h_RP, e_RP, ω = p
    R, C, P = u

    # setup the density dependent preference
    Ω = p.pref(u, p, t)

    du[1] = r * R * (1 - R / K) - a_RC * R * C / (1 + a_RC * h_RC * R) - Ω * a_RP * R * P / (1 + Ω * a_RP * h_RP * R + (1 - Ω) * a_CP * h_CP * C)
    du[2] = e_RC * a_RC * R * C / (1 + a_RC * h_RC * R) - (1 - Ω) * a_CP * C * P / (1 + Ω * a_RP * h_RP * R + (1 - Ω) * a_CP * h_CP * C) - m_C * C
    du[3] = (e_RP * Ω * a_RP * R * P + e_CP * (1 - Ω) * a_CP * C * P) / (1 + Ω * a_RP * h_RP * R + (1 - Ω) * a_CP * h_CP * C) - m_P * P

    return du
end
