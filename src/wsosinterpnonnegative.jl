#=
interpolant-basis weighted sum-of-squares polynomials (real or real-valued complex)
parametrized by vector of matrices `Ps` derived from interpolant basis and domain constraints

primal cone
w : ∃ Sₗ ⪰ 0 ∀ l, w = ∑ₗ diag(Pₗ Sₗ Pₗ')

dual cone
w : Pₗ' Diagonal(w) Pₗ ⪰ 0 ∀ l
=#

mutable struct WSOSInterpNonnegative{C <: RealCompF} <: Cone
    Ps::Vector{Matrix{C}}
    d::Int
    oa_s::Vector{AE}
    s::Vector{RealF}
    WSOSInterpNonnegative{C}() where {C <: RealCompF} = new{C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.WSOSInterpNonnegativeCone{RealF, C},
    ::Bool,
) where {C <: RealCompF}
    cache = WSOSInterpNonnegative{C}()
    @assert !cone.use_dual # TODO
    cache.oa_s = oa_s
    cache.Ps = cone.Ps
    d = cache.d = MOI.dimension(cone)
    @assert length(oa_s) == d
    @assert size(cache.Ps[1], 1) == d
    return cache
end

function MOIPajarito.Cones.add_init_cuts(cache::WSOSInterpNonnegative, oa_model::JuMP.Model)
    # wᵢ ≥ 0
    JuMP.@constraint(oa_model, cache.oa_s .>= 0)
    return cache.d
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::WSOSInterpNonnegative,
    oa_model::JuMP.Model,
)
    # TODO ????
    cut = dot_expr(z, cache.oa_s, oa_model)
    return [cut]
end

function MOIPajarito.Cones.get_sep_cuts(::WSOSInterpNonnegative, ::JuMP.Model)
    # check s ∉ K
    # TODO don't know a fast separation oracle
    @warn("no separation oracle implemented for WSOSInterpNonnegative", maxlog = 1)
    return AE[]
end
