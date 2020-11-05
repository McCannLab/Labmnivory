using LinearAlgebra: eigvals
using ForwardDiff

function rhs(u, p)
    du = similar(u)
    model!(du, u, p, zero(u))
    return du
end

cmat(u, p) = ForwardDiff.jacobian(x -> rhs(x, p), u)
λ1_stability(u, p) = maximum(real.(eigvals(cmat(u, p))))
