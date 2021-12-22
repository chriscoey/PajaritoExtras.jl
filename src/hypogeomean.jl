#=
hypograph of geometric mean
(u, w) : wᵢ ≥ 0, u ≤ (∏ᵢ wᵢ)^(1/d)

dual cone
(u, w) : u ≤ 0, wᵢ ≥ 0, u ≥ -d (∏ᵢ wᵢ)^(1/d)

extended formulation
NOTE slightly different to EF in solver paper
∃ θ ≥ u, ∃ λ, Σᵢ λᵢ ≥ 0, (λᵢ, θ, wᵢ) ∈ HypoPerLog
i.e. θ ≥ 0, wᵢ ≥ 0, λᵢ ≤ θ log(wᵢ / θ)
dual of HypoPerLog is (p, q, r) : p ≤ 0, r ≥ 0, q ≥ p * (log(-r / p) + 1)
=#

mutable struct HypoGeoMeanCache{E <: Extender} <: ConeCache
    cone::Hypatia.HypoGeoMeanCone{Float64}
    oa_s::Vector{AE}
    s::Vector{Float64}
    d::Int
    θ::VR
    λ::Vector{VR}
    HypoGeoMeanCache{E}() where {E <: Extender} = new{E}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.HypoGeoMeanCone{Float64},
    extend::Bool,
)
    dim = MOI.dimension(cone)
    @assert dim == length(oa_s)
    d = dim - 1
    E = extender(extend, d)
    cache = HypoGeoMeanCache{E}()
    cache.cone = cone
    cache.oa_s = oa_s
    cache.d = d
    return cache
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::HypoGeoMeanCache,
    oa_model::JuMP.Model,
)
    return _get_cuts(z[2:end], cache, oa_model)
end

function MOIPajarito.Cones.get_sep_cuts(cache::HypoGeoMeanCache, oa_model::JuMP.Model)
    s = cache.s
    us = s[1]
    @views ws = s[2:end]
    ws_geom = geomean(ws)
    # check s ∉ K
    if ws_geom - us > -1e-7 # TODO option
        return AE[]
    end

    # gradient cut is (-1, -ws_geom / d ./ w)
    c1 = ws_geom / -d
    r = [c1 / max(w_i, 1e-7) for w_i in ws]
    return _get_cuts(r, cache, oa_model)
end

# unextended formulation

function MOIPajarito.Cones.add_init_cuts(
    cache::HypoGeoMeanCache{Unextended},
    oa_model::JuMP.Model,
)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end] # TODO cache?
    d = cache.d
    # variable bounds wᵢ ≥ 0 and cut (-d, e)
    JuMP.@constraints(oa_model, begin
        [i in 1:d], w[i] >= 0
        -d * u + sum(w) >= 0
    end)
    return d + 1
end

function _get_cuts(
    r::Vector{Float64},
    cache::HypoGeoMeanCache{Unextended},
    oa_model::JuMP.Model,
)
    # strengthened cut is (-d * geom(r), r)
    clean_array!(r) && return AE[]
    p = -cache.d * geomean(r)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    cut = JuMP.@expression(oa_model, p * u + JuMP.dot(r, w))
    return [cut]
end

# extended formulation

MOIPajarito.Cones.num_ext_variables(cache::HypoGeoMeanCache{Extended}) = 1 + cache.d

function MOIPajarito.Cones.extend_start(
    cache::HypoGeoMeanCache{Extended},
    s_start::Vector{Float64},
)
    u_start = s_start[1]
    w_start = s_start[2:end]
    w_geom = geomean(w_start)
    @assert w_geom - u_start >= -1e-7 # TODO
    if u_start < 1e-8
        return zeros(1 + cache.d)
    end
    λ_start = [u_start * log(w_i / u_start) for w_i in w_start]
    @assert sum(λ_start) >= -eps()
    return vcat(u_start, λ_start)
end

function MOIPajarito.Cones.setup_auxiliary(
    cache::HypoGeoMeanCache{Extended},
    oa_model::JuMP.Model,
)
    @assert cache.d >= 2
    θ = cache.θ = JuMP.@variable(oa_model, lower_bound = 0)
    u = cache.oa_s[1]
    JuMP.@constraint(oa_model, θ >= u)
    λ = cache.λ = JuMP.@variable(oa_model, [1:(cache.d)])
    JuMP.@constraint(oa_model, sum(λ) >= 0)
    return vcat(θ, λ)
end

function MOIPajarito.Cones.add_init_cuts(
    cache::HypoGeoMeanCache{Extended},
    oa_model::JuMP.Model,
)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    d = cache.d
    θ = cache.θ
    λ = cache.λ
    # variable bounds wᵢ ≥ 0 and cut (-d, e)
    # disaggregated cut on (λᵢ, θ, wᵢ) is (-1, -1, 1)
    # TODO check implication in math
    JuMP.@constraints(oa_model, begin
        [i in 1:d], w[i] >= 0
        [i in 1:d], -λ[i] - θ + w[i] >= 0
    end)
    return 2d
end

function _get_cuts(
    r::Vector{Float64},
    cache::HypoGeoMeanCache{Extended},
    oa_model::JuMP.Model,
)
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
        cut = JuMP.@expression(oa_model, p * λ[i] + q * θ + r_i * w[i])
        push!(cuts, cut)
    end
    return cuts
end
