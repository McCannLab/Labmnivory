using Parameters
using DifferentialEquations
using PyPlot

@with_kw mutable struct OmnPar
    # Logistic Parameters
    r = 2.0
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
end

function model!(du, u, p, t)
    @unpack r, K = p
    @unpack a_RC, h_RC, e_RC, m_C = p
    @unpack a_CP, h_CP, e_CP, m_P = p
    @unpack a_RP, h_RP, e_RP, ω = p

    R, C, P = u
    du[1] = r * R * (1 - R / K) - a_RC * R * C / (1 + a_RC * h_RC * R) - ω * a_RP * R * P / (1 + a_RP * h_RP * R + a_CP * h_CP * C)
    du[2] = e_RC * a_RC * R * C / (1 + a_RC * h_RC * R) - (1 - ω) * a_CP * C * P / (1 + a_RP * h_RP * R + a_CP * h_CP * C) - m_C * C
    du[3] = e_CP * (1 - ω) * a_CP * C * P / (1 + a_CP * h_CP * C) + e_RP * ω * a_RP * R * P / (1 + a_RP * h_RP * R) - m_P * P

    return du
end

let
    u0 = [1.0, 0.5, 0.1]
    t_span = (0.0, 1000.0)
    par = OmnPar(ω = 0.0)

    prob = ODEProblem(model!, u0, t_span, par)
    sol = solve(prob)

    plot(sol.t, sol.u)
end
