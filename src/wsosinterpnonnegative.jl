#=
interpolant-basis weighted sum-of-squares polynomials (real or real-valued complex)
parametrized by vector of matrices `Ps` derived from interpolant basis and domain constraints

primal cone
w : ∃ Sₗ ⪰ 0 ∀ l, w = ∑ₗ diag(Pₗ Sₗ Pₗ')

dual cone
w : Pₗ' Diagonal(w) Pₗ ⪰ 0 ∀ l
=#

mutable struct WSOSInterpNonnegative{D <: PrimDual, C <: RealCompF} <: Cone
    oa_s::Vector{AE}
    Ps::Vector{Matrix{C}}
    d::Int
    WSOSInterpNonnegative{D, C}() where {D <: PrimDual, C <: RealCompF} = new{D, C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.WSOSInterpNonnegativeCone{RealF, C},
    ::Bool,
) where {C <: RealCompF}
    D = primal_or_dual(cone.use_dual)
    cache = WSOSInterpNonnegative{D, C}()
    cache.oa_s = oa_s
    cache.Ps = cone.Ps
    d = cache.d = MOI.dimension(cone)
    @assert length(oa_s) == d
    @assert size(cache.Ps[1], 1) == d
    return cache
end

# primal cone functions

function MOIPajarito.Cones.add_init_cuts(
    cache::WSOSInterpNonnegative{Prim},
    oa_model::JuMP.Model,
)
    # wᵢ ≥ 0
    JuMP.@constraint(oa_model, cache.oa_s .>= 0)
    return cache.d
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::WSOSInterpNonnegative{Prim},
    oa_model::JuMP.Model,
)
    cut = dot_expr(z, cache.oa_s, oa_model)
    return [cut]
end

function MOIPajarito.Cones.get_sep_cuts(
    ::Vector{RealF},
    ::WSOSInterpNonnegative{Prim},
    ::JuMP.Model,
)
    @warn("no separation oracle implemented for WSOSInterpNonnegative", maxlog = 1)
    return AE[]
end

# dual cone functions


#=
interpolation-based weighted-sum-of-squares (multivariate) polynomial cone parametrized by interpolation points ipwt
cuts are <V*V^T, Λ> for eigenvectors V corresponding to negative eigenvalues of Λ
Λ is given by # TODO
=#

# struct WSOSPolyInterpCone <: MOI.AbstractVectorSet
#     dimension::Int
#     ipwt::Vector{Matrix{Float64}}
#     is_dual::Bool
# end
# WSOSPolyInterpCone(dimension::Int, ipwt::Vector{Matrix{Float64}}) = WSOSPolyInterpCone(dimension, ipwt, false)
# dimension(cone::WSOSPolyInterpCone) = cone.dimension

# function get_init_cuts(cone::WSOSPolyInterpCone)
#     cuts = Vector{Float64}[]
#     U = cone.dimension

#     for w in eachindex(cone.ipwt)
#         P = cone.ipwt[w]
#         L = size(P, 2)

#         PiPj = Matrix{Vector{Float64}}(undef, L, L)
#         for i in 1:L, j in 1:i
#             PiPj[i, j] = P[:, i] .* P[:, j]
#         end

#         # diagonal of Λ is nonnegative
#         for i in 1:L
#             push!(cuts, PiPj[i, i])
#         end

#         # 2x2 principal minors of Λ satisfy linearizations of PSD/RSOC condition
#         for i in 1:L, j in 1:(i - 1)
#             cut = PiPj[i, i] + PiPj[j, j] + 2.0 * PiPj[i, j]
#             push!(cuts, cut)
#             cut = PiPj[i, i] + PiPj[j, j] - 2.0 * PiPj[i, j]
#             push!(cuts, cut)
#         end
#     end

#     return cuts
# end

# # TODO maybe faster if store the PiPj matrix calculated for initial cuts and re-use here
# function get_sep_cuts(x::Vector{Float64}, cone::WSOSPolyInterpCone)
#     cuts = Vector{Float64}[]
#     U = cone.dimension

#     for w in eachindex(cone.ipwt)
#         P = cone.ipwt[w]
#         L = size(P, 2)
#         X = Symmetric(P' * Diagonal(x) * P)
#         # X = Matrix{Float64}(undef, L, L)
#         # for i in 1:L, j in 1:i
#         #     X[i, j] = sum(P[u, i] * P[u, j] * x[u] for u in eachindex(x))
#         # end

#         F = eigen!(X, 1:min(L, 5)) # TODO try only getting negative eigenvalues (what lower bound?) vs smallest k eigenvalues
#         # F = eigen!(X, -Inf, -1e-7)
#         # F = eigen!(X)

#         @show F.values

#         for k in eachindex(F.values)
#             if F.values[k] >= -1e-5
#                 continue
#             end

#             V = F.vectors[:, k]
#             PV = P * V
#             cut = abs2.(PV) # V' * P' * Diagonal(vars) * P * V

#             @show dot(x, cut)

#             if dot(x, cut) < -1e-5 # TODO tolerance option
#                 push!(cuts, cut)
#             end
#         end
#     end

#     return cuts
# end