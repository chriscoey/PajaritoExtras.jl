#=
hypograph of geometric mean
(u, w) : wᵢ ≥ 0, u ≤ geom(w)
where geom(w) = (∏ᵢ wᵢ)^(1/d)

dual cone is
(u, w) : u ≤ 0, wᵢ ≥ 0, u ≥ -d * geom(w)

extended formulation
NOTE slightly different to EF in Hypatia solver paper
∃ θ ≥ u, ∃ λ, Σᵢ λᵢ ≥ 0, (λᵢ, θ, wᵢ) ∈ HypoPerLog
i.e. θ ≥ 0, wᵢ ≥ 0, λᵢ ≤ θ log(wᵢ / θ)
dual of HypoPerLog is (p, q, r) : p ≤ 0, r ≥ 0, q ≥ p * (log(-r / p) + 1)
=#

mutable struct HypoGeoMean{E <: NatExt} <: Cache
    oa_s::Vector{AE}
    d::Int
    θ::VR
    λ::Vector{VR}
    HypoGeoMean{E}() where {E <: NatExt} = new{E}()
end

function Pajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.HypoGeoMeanCone{RealF},
    opt::Optimizer,
)
    @assert !cone.use_dual # TODO
    dim = MOI.dimension(cone)
    @assert dim == length(oa_s)
    d = dim - 1
    E = nat_or_ext(opt, d)
    cache = HypoGeoMean{E}()
    cache.oa_s = oa_s
    cache.d = d
    return cache
end

function Pajarito.Cones.get_subp_cuts(z::Vector{RealF}, cache::HypoGeoMean, opt::Optimizer)
    return _get_cuts(z[2:end], cache, opt)
end

function Pajarito.Cones.get_sep_cuts(s::Vector{RealF}, cache::HypoGeoMean, opt::Optimizer)
    us = s[1]
    @views ws = s[2:end]
    if us < opt.tol_feas || geomean(ws) - us > -opt.tol_feas
        return AE[]
    end

    # gradient cut is (-1, geom(w) / d ./ w)
    w_pos = max.(ws, 1e-7)
    c1 = geomean(w_pos) / cache.d
    r = c1 ./ w_pos
    return _get_cuts(r, cache, opt)
end

# unextended formulation

function Pajarito.Cones.add_init_cuts(cache::HypoGeoMean{Nat}, opt::Optimizer)
    # add variable bounds wᵢ ≥ 0
    @views w = cache.oa_s[2:end]
    JuMP.@constraint(opt.oa_model, w .>= 0)
    opt.use_init_fixed_oa || return

    # add cut (-d, e) on (u, w), a linearization at w = e
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, -cache.d * u + sum(w) >= 0)
    return
end

function _get_cuts(r::Vector{RealF}, cache::HypoGeoMean{Nat}, opt::Optimizer)
    # strengthened cut is (-d * geom(r), r)
    clean_array!(r) && return AE[]
    p = -cache.d * geomean(r)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    cut = JuMP.@expression(opt.oa_model, p * u + dot(r, w))
    return [cut]
end

# extended formulation

Pajarito.Cones.num_ext_variables(cache::HypoGeoMean{Ext}) = 1 + cache.d

function Pajarito.Cones.extend_start(
    cache::HypoGeoMean{Ext},
    s_start::Vector{RealF},
    opt::Optimizer,
)
    @views w_start = s_start[2:end]
    if any(<(1e-9), w_start)
        return zeros(1 + cache.d)
    end
    θ_start = geomean(w_start)
    λ_start = [θ_start * log(w_i / θ_start) for w_i in w_start]
    return vcat(θ_start, λ_start)
end

function Pajarito.Cones.setup_auxiliary(cache::HypoGeoMean{Ext}, opt::Optimizer)
    @assert cache.d >= 2
    θ = cache.θ = JuMP.@variable(opt.oa_model, lower_bound = 0)
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, θ >= u)
    λ = cache.λ = JuMP.@variable(opt.oa_model, [1:(cache.d)])
    JuMP.@constraint(opt.oa_model, sum(λ) >= 0)
    return vcat(θ, λ)
end

function Pajarito.Cones.add_init_cuts(cache::HypoGeoMean{Ext}, opt::Optimizer)
    # add variable bounds wᵢ ≥ 0
    @views w = cache.oa_s[2:end]
    JuMP.@constraint(opt.oa_model, w .>= 0)
    opt.use_init_fixed_oa || return

    # add disaggregated cuts (-1, -1, 1) on (λᵢ, θ, wᵢ) TODO check
    θ = cache.θ
    λ = cache.λ
    JuMP.@constraint(opt.oa_model, [i in 1:(cache.d)], -λ[i] - θ + w[i] >= 0)
    return
end

function _get_cuts(r::Vector{RealF}, cache::HypoGeoMean{Ext}, opt::Optimizer)
    clean_array!(r) && return AE[]
    p = -geomean(r)
    @views w = cache.oa_s[2:end]
    θ = cache.θ
    λ = cache.λ
    cuts = AE[]
    for i in 1:(cache.d)
        r_i = r[i]
        iszero(r_i) && continue # TODO ?
        # strengthened disaggregated cut on (λᵢ, θ, wᵢ) is
        # p = -geom(r), (p, p * (log(-rᵢ / p) + 1), rᵢ)
        # TODO check math
        q = p * (log(-r_i / p) + 1)
        cut = JuMP.@expression(opt.oa_model, p * λ[i] + q * θ + r_i * w[i])
        push!(cuts, cut)
    end
    return cuts
end
