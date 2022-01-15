#=
epigraph of perspective of squared Euclidean norm (AKA rotated second-order cone)
(u, v, w) : v ≥ 0, u ≥ v / 2 * ‖w / v‖² = ‖w‖² / 2v
self-dual

extended formulation
∃ λ, u ≥ Σᵢ λᵢ, (λᵢ, v, wᵢ) ∈ EpiPerSquare
i.e. λᵢ ≥ 0, 2 λᵢ v ≥ wᵢ²
=#

mutable struct EpiPerSquare{E <: NatExt} <: Cache
    oa_s::Vector{AE}
    d::Int
    λ::Vector{VR}
    EpiPerSquare{E}() where {E <: NatExt} = new{E}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiPerSquareCone{RealF},
    extend::Bool,
)
    dim = MOI.dimension(cone)
    @assert dim == length(oa_s)
    d = dim - 2
    E = nat_or_ext(extend, d)
    cache = EpiPerSquare{E}()
    cache.oa_s = oa_s
    cache.d = d
    return cache
end

per_square(v::RealF, w::AbstractVector{RealF}) = sum(w_i / 2v * w_i for w_i in w)

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::EpiPerSquare,
    opt::Optimizer,
)
    return _get_cuts(z[2], z[3:end], cache, opt)
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::EpiPerSquare,
    opt::Optimizer,
)
    us = s[1]
    vs = s[2]
    @views ws = s[3:end]
    rhs = per_square(vs, ws)
    if us - rhs > -1e-7 # TODO option
        return AE[]
    end

    # gradient cut is (1, rhs / vs, -ws / vs)
    q = rhs / vs
    r = ws / -vs
    return _get_cuts(q, r, cache, opt)
end

# unextended formulation

function MOIPajarito.Cones.add_init_cuts(cache::EpiPerSquare{Nat}, opt::Optimizer)
    u = cache.oa_s[1]
    v = cache.oa_s[2]
    @views w = cache.oa_s[3:end] # TODO cache?
    d = cache.d
    # variable bounds u ≥ 0, v ≥ 0 and cut (1, 2, ±2eᵢ)
    JuMP.@constraints(opt.oa_model, begin
        u >= 0
        v >= 0
        [i in 1:d], u + 2v + 2w[i] >= 0
        [i in 1:d], u + 2v - 2w[i] >= 0
    end)
    return 2 + 2d
end

function _get_cuts(q::RealF, r::Vector{RealF}, cache::EpiPerSquare{Nat}, opt::Optimizer)
    clean_array!(r) && return AE[]
    # strengthened cut is (‖r‖² / 2q, q, r)
    p = per_square(q, r)
    u = cache.oa_s[1]
    v = cache.oa_s[2]
    @views w = cache.oa_s[3:end]
    cut = JuMP.@expression(opt.oa_model, p * u + q * v + JuMP.dot(r, w))
    return [cut]
end

# extended formulation

MOIPajarito.Cones.num_ext_variables(cache::EpiPerSquare{Ext}) = cache.d

function MOIPajarito.Cones.extend_start(cache::EpiPerSquare{Ext}, s_start::Vector{RealF})
    u_start = s_start[1]
    v_start = s_start[2]
    w_start = s_start[3:end]
    @assert u_start - per_square(v_start, w_start) >= -1e-7 # TODO
    if max(u_start, v_start) < 1e-8
        return zeros(cache.d)
    end
    return [w_i / 2u_start * w_i for w_i in w_start]
end

function MOIPajarito.Cones.setup_auxiliary(cache::EpiPerSquare{Ext}, opt::Optimizer)
    @assert cache.d >= 2
    λ = cache.λ = JuMP.@variable(opt.oa_model, [1:(cache.d)], lower_bound = 0)
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, u >= sum(λ))
    return λ
end

function MOIPajarito.Cones.add_init_cuts(cache::EpiPerSquare{Ext}, opt::Optimizer)
    u = cache.oa_s[1]
    v = cache.oa_s[2]
    @views w = cache.oa_s[3:end]
    d = cache.d
    λ = cache.λ
    # variable bounds u ≥ 0, v ≥ 0 and cut (1, 2, ±2eᵢ)
    # disaggregated cut on (λᵢ, v, wᵢ) is (1, 2, ±2) since u ≥ λᵢ
    JuMP.@constraints(opt.oa_model, begin
        u >= 0
        v >= 0
        [i in 1:d], λ[i] + 2v + 2w[i] >= 0
        [i in 1:d], λ[i] + 2v - 2w[i] >= 0
    end)
    return 1 + 2d
end

function _get_cuts(q::RealF, r::Vector{RealF}, cache::EpiPerSquare{Ext}, opt::Optimizer)
    clean_array!(r) && return AE[]
    p = per_square(q, r)
    v = cache.oa_s[2]
    @views w = cache.oa_s[3:end]
    λ = cache.λ
    cuts = AE[]
    for i in 1:(cache.d)
        r_i = r[i]
        iszero(r_i) && continue
        # p = ‖r‖² / 2q, strengthened disaggregated cut on (λᵢ, v, wᵢ) is (p, rᵢ² / 2p, rᵢ)
        # TODO check math
        q_i = r_i / 2p * r_i
        cut = JuMP.@expression(opt.oa_model, p * λ[i] + q_i * v + r_i * w[i])
        push!(cuts, cut)
    end
    return cuts
end
