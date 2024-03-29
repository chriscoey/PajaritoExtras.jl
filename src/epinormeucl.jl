#=
epigraph of Euclidean norm (AKA second-order cone)
(u, w) : u ≥ ‖w‖
self-dual

extended formulation
∃ λ, u ≥ 2 Σᵢ λᵢ, (u, λᵢ, wᵢ) ∈ EpiPerSquare
i.e. λᵢ ≥ 0, 2 u λᵢ ≥ wᵢ²
=#

mutable struct EpiNormEucl{E <: NatExt} <: Cache
    oa_s::Vector{AE}
    d::Int
    λ::Vector{VR}
    EpiNormEucl{E}() where {E <: NatExt} = new{E}()
end

function Pajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiNormEuclCone{RealF},
    opt::Optimizer,
)
    dim = MOI.dimension(cone)
    @assert dim == length(oa_s)
    d = dim - 1
    E = nat_or_ext(opt, d)
    cache = EpiNormEucl{E}()
    cache.oa_s = oa_s
    cache.d = d
    return cache
end

function Pajarito.Cones.get_subp_cuts(z::Vector{RealF}, cache::EpiNormEucl, opt::Optimizer)
    return _get_cuts(z[2:end], cache, opt)
end

function Pajarito.Cones.get_sep_cuts(s::Vector{RealF}, cache::EpiNormEucl, opt::Optimizer)
    us = s[1]
    @views ws = s[2:end]
    ws_norm = norm(ws)
    if us - ws_norm > -opt.tol_feas
        return AE[]
    end

    # gradient cut is (1, -ws / ‖ws‖)
    r = ws / -ws_norm
    return _get_cuts(r, cache, opt)
end

# unextended formulation

function Pajarito.Cones.add_init_cuts(cache::EpiNormEucl{Nat}, opt::Optimizer)
    # add variable bound
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, u >= 0)
    opt.use_init_fixed_oa || return

    # add cuts u ≥ |wᵢ|
    @views w = cache.oa_s[2:end]
    JuMP.@constraints(opt.oa_model, begin
        u .>= w
        u .>= -w
    end)
    return
end

function _get_cuts(r::Vector{RealF}, cache::EpiNormEucl{Nat}, opt::Optimizer)
    # strengthened cut is (‖r‖, r)
    clean_array!(r) && return AE[]
    p = norm(r)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    cut = JuMP.@expression(opt.oa_model, p * u + dot(r, w))
    return [cut]
end

# extended formulation

Pajarito.Cones.num_ext_variables(cache::EpiNormEucl{Ext}) = cache.d

function Pajarito.Cones.extend_start(
    cache::EpiNormEucl{Ext},
    s_start::Vector{RealF},
    opt::Optimizer,
)
    u_start = s_start[1]
    w_start = s_start[2:end]
    if u_start < 1e-9
        return zeros(cache.d)
    end
    return [w_i / 2u_start * w_i for w_i in w_start]
end

function Pajarito.Cones.setup_auxiliary(cache::EpiNormEucl{Ext}, opt::Optimizer)
    @assert cache.d >= 2
    λ = cache.λ = JuMP.@variable(opt.oa_model, [1:(cache.d)], lower_bound = 0)
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, u >= 2 * sum(λ))
    return λ
end

function Pajarito.Cones.add_init_cuts(cache::EpiNormEucl{Ext}, opt::Optimizer)
    # add variable bound
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, u >= 0)
    opt.use_init_fixed_oa || return

    # add disaggregated cuts (1, 2, ±2) on (u, λᵢ, wᵢ), implying u ≥ |wᵢ|
    @views w = cache.oa_s[2:end]
    d = cache.d
    λ = cache.λ
    JuMP.@constraints(opt.oa_model, begin
        [i in 1:d], u + 2λ[i] + 2w[i] >= 0
        [i in 1:d], u + 2λ[i] - 2w[i] >= 0
    end)
    return
end

function _get_cuts(r::Vector{RealF}, cache::EpiNormEucl{Ext}, opt::Optimizer)
    clean_array!(r) && return AE[]
    p = norm(r)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    λ = cache.λ
    cuts = AE[]
    for i in 1:(cache.d)
        r_i = r[i]
        iszero(r_i) && continue
        # strengthened disaggregated cut on (u, λᵢ, wᵢ) is (rᵢ² / 2‖r‖, ‖r‖, rᵢ)
        cut = JuMP.@expression(opt.oa_model, r_i^2 / 2p * u + p * λ[i] + r_i * w[i])
        push!(cuts, cut)
    end
    return cuts
end
