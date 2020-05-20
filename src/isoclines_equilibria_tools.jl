include("basic_omnivory_module.jl")
include("eigenvalue_utilities.jl")
using SymPy

#### Calculating equilibria
@vars u1 u2 u3
@vars r K a_RC h_RC e_RC m_C a_CP h_CP e_CP m_P a_RP h_RP e_RP ω
# We need to make a symbolic parameter list, as `LVPar` is numeric
chain_par = Dict(
    :r => r,
    :K => K,
    :a_RC => a_RC,
    :h_RC => h_RC,
    :e_RC => e_RC,
    :m_C => m_C,
    :a_CP => a_CP,
    :h_CP => h_CP,
    :e_CP => e_CP,
    :m_P => m_P,
    :a_RP => 0,
    :h_RP => h_RP,
    :e_RP => e_RP,
    :ω => 0);

omni_par = Dict(
    :r => r,
    :K => K,
    :a_RC => a_RC,
    :h_RC => h_RC,
    :e_RC => e_RC,
    :m_C => m_C,
    :a_CP => a_CP,
    :h_CP => h_CP,
    :e_CP => e_CP,
    :m_P => m_P,
    :a_RP => a_RP,
    :h_RP => h_RP,
    :e_RP => e_RP,
    :ω => ω);

function model_noforce(du, u, p, t)
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

function model_noforce(u, p)
    du = similar(u)
    model_noforce(du, u, p, zero(u))
    return du
end

## Food Chain
fc1, fc2, fc3 = model_noforce([u1, u2, u3], chain_par)
#-
SymPy.solve(fc1, u1)
#-
SymPy.solve(fc1, u2)
#-
SymPy.solve(fc2, u1)
#-
SymPy.solve(fc2, u2)
#-
SymPy.solve(fc2, u3)
#-
SymPy.solve(fc3, u2)
#-
SymPy.solve(fc3, u3)

#NOTE - no solution for fc1 in terms of u3 and fc3 in terms of u1, due to no omnivory
## Omnivory model
o1, o2, o3 = model_noforce([u1, u2, u3], omni_par)
#-
SymPy.solve(o1, u1)
#-
SymPy.solve(o1, u2)
#-
SymPy.solve(o1, u3)
#-
SymPy.solve(o2, u1)
#-
SymPy.solve(o2, u2)
#-
SymPy.solve(o2, u3)
#-
SymPy.solve(o3, u1)
#-
SymPy.solve(o3, u2)
#-
SymPy.solve(o3, u3)

#### Calculating and plotting isoclines

## Food Chain

## 3d
function res_iso_fc(u1, p)
    @unpack r, K = p
    @unpack a_RC, h_RC = p
    return r * (K * a_RC * h_RC * u1 + K - a_RC * h_RC * u1^2 - u1) / (K * a_RC)
end
#-
function con_iso_fc(u1, u2, p)
    @unpack r, K = p
    @unpack a_RC, h_RC, e_RC, m_C = p
    @unpack a_CP, h_CP = p
    return (-a_CP*a_RC*h_RC*u1*u3 - a_CP*u3 + a_RC*e_RC*u1 - a_RC*h_RC*m_C*u1 - m_C)/(a_CP*h_CP*(-a_RC*e_RC*u1 + a_RC*h_RC*m_C*u1 + m_C))
end

function pred_iso_fc(p)
    @unpack a_CP, h_CP, e_CP, m_P = p
    return m_P/(a_CP*(e_CP - h_CP*m_P))
end

pred_iso_fc(ModelPar())

u1s = range(0, 5, length = 100)
u2s = range(0, 3, length = 100)
xgrid = repeat(u1s',100,1)
ygrid = repeat(u2s,1,100)

z = zeros(100,100)

for (i,u1) in enumerate(u1s)
    for (j,u2) in enumerate(u2s)
        z[i:i,j:j] .=  con_iso_fc(u1, u2, ModelPar())
    end
end

let
    u1s = range(0, 5, length = 100)
    u2s = range(0, 3, length = 100)
    test = figure()
    plot_surface(xgrid, ygrid, z)
    plot(u1s, [res_iso_fc(u1, ModelPar()) for u1 in u1s])
    plot(repeat([pred_iso_fc(ModelPar())], 100), u2s)
    ylim(0,3)
    # zlim(0,1)
    xlabel("Resource")
    ylabel("Consumer")
    return test
end


let
    p = ModelPar()
    u1s = range(0, 5, length = 100)
    u2s = range(0, 3, length = 100)

    fig = figure()
    plot_surface(u1s, u2s, con_iso_fc(u1s, u2s, ModelPar()), st=:surface) #, pred_iso_fc(p), st=:surface)
    #plot!(u2s, [iso2(u2, p) for u2 in u2s], label = "u2 = 0")
    #xlims!(0, 3)
    #ylims!(0, 5)
    #xlabel!("u2")
    #ylabel!("u1")
    return fig
end

## 2d where P is a parameter (taking snapshots of isoclines at different P values)

## Omnivory model
