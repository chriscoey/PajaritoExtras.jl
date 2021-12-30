#=
epigraph of Euclidean norm (AKA second-order cone)
(u, w) : u ≥ ‖w‖
self-dual

extended formulation
∃ λ, u ≥ 2 Σᵢ λᵢ, (u, λᵢ, wᵢ) ∈ EpiPerSquare
i.e. λᵢ ≥ 0, 2 u λᵢ ≥ wᵢ²
=#

mutable struct EpiNormEucl{E <: NatExt} <: Cone
    oa_s::Vector{AE}
    s::Vector{RealF}
    d::Int
    λ::Vector{VR}
    EpiNormEucl{E}() where {E <: NatExt} = new{E}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiNormEuclCone{RealF},
    extend::Bool,
)
    dim = MOI.dimension(cone)
    @assert dim == length(oa_s)
    d = dim - 1
    E = extender(extend, d)
    cache = EpiNormEucl{E}()
    cache.oa_s = oa_s
    cache.d = d
    return cache
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::EpiNormEucl,
    oa_model::JuMP.Model,
)
    return _get_cuts(z[2:end], cache, oa_model)
end

function MOIPajarito.Cones.get_sep_cuts(cache::EpiNormEucl, oa_model::JuMP.Model)
    s = cache.s
    us = s[1]
    @views ws = s[2:end]
    ws_norm = LinearAlgebra.norm(ws)
    # check s ∉ K
    if us - ws_norm > -1e-7 # TODO option
        return AE[]
    end

    # gradient cut is (1, -ws / ‖ws‖)
    r = -inv(ws_norm) * ws
    return _get_cuts(r, cache, oa_model)
end

# unextended formulation

function MOIPajarito.Cones.add_init_cuts(cache::EpiNormEucl{Nat}, oa_model::JuMP.Model)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end] # TODO cache?
    d = cache.d
    # u ≥ 0, u ≥ |wᵢ|
    JuMP.@constraints(oa_model, begin
        u >= 0
        [i in 1:d], u >= w[i]
        [i in 1:d], u >= -w[i]
    end)
    return 1 + 2d
end

function _get_cuts(r::Vector{RealF}, cache::EpiNormEucl{Nat}, oa_model::JuMP.Model)
    # strengthened cut is (‖r‖, r)
    clean_array!(r) && return AE[]
    p = LinearAlgebra.norm(r)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    cut = JuMP.@expression(oa_model, p * u + JuMP.dot(r, w))
    return [cut]
end

# extended formulation

MOIPajarito.Cones.num_ext_variables(cache::EpiNormEucl{Ext}) = cache.d

function MOIPajarito.Cones.extend_start(cache::EpiNormEucl{Ext}, s_start::Vector{RealF})
    u_start = s_start[1]
    w_start = s_start[2:end]
    @assert u_start - LinearAlgebra.norm(w_start) >= -1e-7 # TODO
    if u_start < 1e-8
        return zeros(cache.d)
    end
    return [w_i / 2u_start * w_i for w_i in w_start]
end

function MOIPajarito.Cones.setup_auxiliary(cache::EpiNormEucl{Ext}, oa_model::JuMP.Model)
    @assert cache.d >= 2
    λ = cache.λ = JuMP.@variable(oa_model, [1:(cache.d)], lower_bound = 0)
    u = cache.oa_s[1]
    JuMP.@constraint(oa_model, u >= 2 * sum(λ))
    return λ
end

function MOIPajarito.Cones.add_init_cuts(cache::EpiNormEucl{Ext}, oa_model::JuMP.Model)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    d = cache.d
    λ = cache.λ
    # u ≥ 0, u ≥ |wᵢ|
    # disaggregated cut on (u, λᵢ, wᵢ) is (1, 2, ±2)
    JuMP.@constraints(oa_model, begin
        u >= 0
        [i in 1:d], u + 2 * λ[i] + 2 * w[i] >= 0
        [i in 1:d], u + 2 * λ[i] - 2 * w[i] >= 0
    end)
    return 1 + 2d
end

function _get_cuts(r::Vector{RealF}, cache::EpiNormEucl{Ext}, oa_model::JuMP.Model)
    clean_array!(r) && return AE[]
    p = LinearAlgebra.norm(r)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    λ = cache.λ
    cuts = AE[]
    for i in 1:(cache.d)
        r_i = r[i]
        iszero(r_i) && continue
        # strengthened disaggregated cut on (u, λᵢ, wᵢ) is (rᵢ² / 2‖r‖, ‖r‖, rᵢ)
        cut = JuMP.@expression(oa_model, r_i^2 / 2p * u + p * λ[i] + r_i * w[i])
        push!(cuts, cut)
    end
    return cuts
end
