using LinearAlgebra
using ForwardDiff

function rhs(u, p)
    du = similar(u)
    model!(du, u, p, zero(u))
    return du
end

cmat(u, p) = ForwardDiff.jacobian(x -> rhs(x, p), u)
λ1_stability(p) = maximum(real.(eigvals(cmat(p))))
