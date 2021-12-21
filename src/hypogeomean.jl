#=
hypograph of geometric mean
(u, w) : wᵢ ≥ 0, u ≤ (∏ᵢ wᵢ)^(1/d)

extended formulation
∃ θ ≥ 0, ∃ λ, Σᵢ λᵢ ≥ 0, (λᵢ, u + θ, wᵢ) ∈ HypoPerLog
i.e. u + θ ≥ 0, wᵢ ≥ 0, λᵢ ≤ (u + θ) log(wᵢ / (u + θ))

dual cone
(u, w) : u ≤ 0, wᵢ ≥ 0, u ≥ -d (∏ᵢ wᵢ)^(1/d)
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

# # extended formulation

# MOIPajarito.Cones.num_ext_variables(cache::HypoGeoMeanCache{Extended}) = cache.d

# function MOIPajarito.Cones.extend_start(
#     cache::HypoGeoMeanCache{Extended},
#     s_start::Vector{Float64},
# )
#     u_start = s_start[1]
#     w_start = s_start[2:end]
#     @assert u_start - LinearAlgebra.norm(w_start) >= -1e-7 # TODO
#     if u_start < 1e-8
#         return zeros(cache.d)
#     end
#     return [w_i / 2u_start * w_i for w_i in w_start]
# end

# function MOIPajarito.Cones.setup_auxiliary(
#     cache::HypoGeoMeanCache{Extended},
#     oa_model::JuMP.Model,
# )
#     @assert cache.d >= 2
#     λ = cache.λ = JuMP.@variable(oa_model, [1:(cache.d)], lower_bound = 0)
#     u = cache.oa_s[1]
#     JuMP.@constraint(oa_model, u >= 2 * sum(λ))
#     return λ
# end

# function MOIPajarito.Cones.add_init_cuts(
#     cache::HypoGeoMeanCache{Extended},
#     oa_model::JuMP.Model,
# )
#     u = cache.oa_s[1]
#     @views w = cache.oa_s[2:end]
#     d = cache.d
#     λ = cache.λ
#     # u ≥ 0, u ≥ |wᵢ|
#     # disaggregated cut on (u, λᵢ, wᵢ) is (1, 2, ±2)
#     JuMP.@constraints(oa_model, begin
#         u >= 0
#         [i in 1:d], u + 2 * λ[i] + 2 * w[i] >= 0
#         [i in 1:d], u + 2 * λ[i] - 2 * w[i] >= 0
#     end)
#     return 1 + 2d
# end

# function _get_cuts(
#     r::Vector{Float64},
#     cache::HypoGeoMeanCache{Extended},
#     oa_model::JuMP.Model,
# )
#     clean_array!(r) && return AE[]
#     p = LinearAlgebra.norm(r)
#     u = cache.oa_s[1]
#     @views w = cache.oa_s[2:end]
#     λ = cache.λ
#     cuts = AE[]
#     for i in 1:(cache.d)
#         r_i = r[i]
#         iszero(r_i) && continue
#         # strengthened disaggregated cut on (u, λᵢ, wᵢ) is (rᵢ² / 2‖r‖, ‖r‖, rᵢ)
#         cut = JuMP.@expression(oa_model, r_i^2 / 2p * u + p * λ[i] + r_i * w[i])
#         push!(cuts, cut)
#     end
#     return cuts
# end
