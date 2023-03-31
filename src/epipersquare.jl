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

function Pajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiPerSquareCone{RealF},
    opt::Optimizer,
)
    dim = MOI.dimension(cone)
    @assert dim == length(oa_s)
    d = dim - 2
    E = nat_or_ext(opt, d)
    cache = EpiPerSquare{E}()
    cache.oa_s = oa_s
    cache.d = d
    return cache
end

function per_square(v::RealF, w::AbstractVector{RealF})
    den = 2 * max(v, 1e-9)
    return sum(w_i / den * w_i for w_i in w)
end

function Pajarito.Cones.get_subp_cuts(z::Vector{RealF}, cache::EpiPerSquare, opt::Optimizer)
    return _get_cuts(z[2], z[3:end], cache, opt)
end

function Pajarito.Cones.get_sep_cuts(s::Vector{RealF}, cache::EpiPerSquare, opt::Optimizer)
    us = s[1]
    vs = s[2]
    @views ws = s[3:end]
    if us - per_square(vs, ws) > -1e-7 # TODO option
        return AE[]
    end

    # gradient cut is (1, ‖w‖² / (2 * vs²), -ws / vs)
    r = ws / -max(vs, 1e-9)
    q = 0.5 * sum(abs2, r)
    return _get_cuts(q, r, cache, opt)
end

# unextended formulation

function Pajarito.Cones.add_init_cuts(cache::EpiPerSquare{Nat}, opt::Optimizer)
    # add variable bounds
    u = cache.oa_s[1]
    v = cache.oa_s[2]
    JuMP.@constraint(opt.oa_model, u >= 0)
    JuMP.@constraint(opt.oa_model, v >= 0)
    opt.use_init_fixed_oa || return

    # add cuts (1, 2, ±2eᵢ) on (u, v, wᵢ)
    @views w = cache.oa_s[3:end]
    d = cache.d
    JuMP.@constraints(opt.oa_model, begin
        [i in 1:d], u + 2v + 2w[i] >= 0
        [i in 1:d], u + 2v - 2w[i] >= 0
    end)
    return
end

function _get_cuts(q::RealF, r::Vector{RealF}, cache::EpiPerSquare{Nat}, opt::Optimizer)
    clean_array!(r) && return AE[]
    # strengthened cut is (‖r‖² / 2q, q, r)
    p = per_square(q, r)
    u = cache.oa_s[1]
    v = cache.oa_s[2]
    @views w = cache.oa_s[3:end]
    cut = JuMP.@expression(opt.oa_model, p * u + q * v + dot(r, w))
    return [cut]
end

# extended formulation

Pajarito.Cones.num_ext_variables(cache::EpiPerSquare{Ext}) = cache.d

function Pajarito.Cones.extend_start(
    cache::EpiPerSquare{Ext},
    s_start::Vector{RealF},
    opt::Optimizer,
)
    u_start = s_start[1]
    v_start = s_start[2]
    w_start = s_start[3:end]
    if max(u_start, v_start) < 1e-8
        return zeros(cache.d)
    end
    return [w_i / 2u_start * w_i for w_i in w_start]
end

function Pajarito.Cones.setup_auxiliary(cache::EpiPerSquare{Ext}, opt::Optimizer)
    @assert cache.d >= 2
    λ = cache.λ = JuMP.@variable(opt.oa_model, [1:(cache.d)], lower_bound = 0)
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, u >= sum(λ))
    return λ
end

function Pajarito.Cones.add_init_cuts(cache::EpiPerSquare{Ext}, opt::Optimizer)
    # add variable bounds
    u = cache.oa_s[1]
    v = cache.oa_s[2]
    JuMP.@constraint(opt.oa_model, u >= 0)
    JuMP.@constraint(opt.oa_model, v >= 0)
    opt.use_init_fixed_oa || return

    # add disaggregated cuts (1, 2, ±2) on (λᵢ, v, wᵢ)
    @views w = cache.oa_s[3:end]
    d = cache.d
    λ = cache.λ
    JuMP.@constraints(opt.oa_model, begin
        [i in 1:d], λ[i] + 2v + 2w[i] >= 0
        [i in 1:d], λ[i] + 2v - 2w[i] >= 0
    end)
    return
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
