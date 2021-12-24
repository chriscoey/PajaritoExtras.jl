# real vector domain ℝᵈ₊

mutable struct VectorEpiPerSepSpectralCache{E <: Extender} <: ConeCache
    cone::Hypatia.EpiPerSepSpectralCone{Float64}
    oa_s::Vector{AE}
    s::Vector{Float64}
    h::SepSpectralFun
    d::Int
    # θ::VR
    # λ::Vector{VR}
    VectorEpiPerSepSpectralCache{E}() where {E <: Extender} = new{E}()
end

function create_sepspectral_cache(::Type{VectorCSqr{Float64}}, d::Int, extend::Bool)
    E = extender(extend, d)
    return VectorEpiPerSepSpectralCache{E}()
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::VectorEpiPerSepSpectralCache,
    oa_model::JuMP.Model,
)
    @show z
return _get_cuts(z[1], z[3:end], cache, oa_model)
    cut = dot_expr(z, cache.oa_s, oa_model)
    return [cut]
end

function MOIPajarito.Cones.get_sep_cuts(
    cache::VectorEpiPerSepSpectralCache,
    oa_model::JuMP.Model,
)
return AE[]

    s = cache.s
    us = s[1]
    vs = s[2]
    ws = s[3:end]
    @assert vs > -1e-7
    @assert all(>(-1e-7), ws)
    # check s ∉ K
    if (vs < 1e-7 && us > -1e-7) || us - per_sepspec(h_val, cache.h, vs, ws) > -1e-7
        return AE[]
    end

    # gradient cut is (1, h⋆(r), r) at r = -h'(ws / vs)
    v_pos = max(vs, 1e-7)
    w_pos = max.(ws, 1e-7)
    p = 1.0
    r = h_grad(w_pos / v_pos, cache.h)
    r .*= -1
    return _get_cuts(p, r, cache, oa_model)
end

# unextended formulation

function MOIPajarito.Cones.add_init_cuts(
    cache::VectorEpiPerSepSpectralCache{Unextended},
    oa_model::JuMP.Model,
)
    u = cache.oa_s[1]
    v = cache.oa_s[2]
    @views w = cache.oa_s[3:end] # TODO cache?
    d = cache.d
    h = cache.h
    # variable bounds v ≥ 0, wᵢ ≥ 0
    JuMP.@constraints(oa_model, begin
        v >= 0
        [i in 1:d], w[i] >= 0
    end)
    # cuts using values of p = 1 and r = r₀ e
    r_vals = init_r_vals(h)
    # for r0 in r_vals
    #     r = fill(r0, d)
    #     q = h_conj(r, h)
    #     JuMP.@constraint(oa_model, u + q * v + JuMP.dot(r, w) >= 0)
    # end
    return 1 + d + length(r_vals)
end

function _get_cuts(
    p::Float64,
    r::Vector{Float64},
    cache::VectorEpiPerSepSpectralCache{Unextended},
    oa_model::JuMP.Model,
)
    # strengthened cut is (p, p * h⋆(r / p), r)
    @assert p > 1e-12
    q = per_sepspec(h_conj, cache.h, p, r)
    @show p
    @show q
    @show r
    println()
    if h_conj_dom_pos(cache.h)
        @assert all(>(1e-12), r)
    end
    z = vcat(p, q, r)
    cut = dot_expr(z, cache.oa_s, oa_model)
    return [cut]
end

# # extended formulation

# MOIPajarito.Cones.num_ext_variables(cache::VectorEpiPerSepSpectralCache{Extended}) = 1 + cache.d

# function MOIPajarito.Cones.extend_start(
#     cache::VectorEpiPerSepSpectralCache{Extended},
#     s_start::Vector{Float64},
# )
#     u_start = s_start[1]
#     w_start = s_start[3:end]
#     w_geom = geomean(w_start)
#     @assert w_geom - u_start >= -1e-7 # TODO
#     if u_start < 1e-8
#         return zeros(1 + cache.d)
#     end
#     λ_start = [u_start * log(w_i / u_start) for w_i in w_start]
#     @assert sum(λ_start) >= -eps()
#     return vcat(u_start, λ_start)
# end

# function MOIPajarito.Cones.setup_auxiliary(
#     cache::VectorEpiPerSepSpectralCache{Extended},
#     oa_model::JuMP.Model,
# )
#     @assert cache.d >= 2
#     θ = cache.θ = JuMP.@variable(oa_model, lower_bound = 0)
#     u = cache.oa_s[1]
#     JuMP.@constraint(oa_model, θ >= u)
#     λ = cache.λ = JuMP.@variable(oa_model, [1:(cache.d)])
#     JuMP.@constraint(oa_model, sum(λ) >= 0)
#     return vcat(θ, λ)
# end

# function MOIPajarito.Cones.add_init_cuts(
#     cache::VectorEpiPerSepSpectralCache{Extended},
#     oa_model::JuMP.Model,
# )
#     @views w = cache.oa_s[3:end]
#     d = cache.d
#     θ = cache.θ
#     λ = cache.λ
#     # variable bounds wᵢ ≥ 0 and cut (-d, e)
#     # disaggregated cut on (λᵢ, θ, wᵢ) is (-1, -1, 1)
#     # TODO check implication in math
#     JuMP.@constraints(oa_model, begin
#         [i in 1:d], w[i] >= 0
#         [i in 1:d], -λ[i] - θ + w[i] >= 0
#     end)
#     return 2d
# end

# function _get_cuts(
#     r::Vector{Float64},
#     cache::VectorEpiPerSepSpectralCache{Extended},
#     oa_model::JuMP.Model,
# )
#     clean_array!(r) && return AE[]
#     p = -geomean(r)
#     @views w = cache.oa_s[3:end]
#     θ = cache.θ
#     λ = cache.λ
#     cuts = AE[]
#     for i in 1:(cache.d)
#         r_i = r[i]
#         iszero(r_i) && continue # TODO ?
#         # strengthened disaggregated cut on (λᵢ, θ, wᵢ) is
#         # p = -geom(r), (p, p * (log(-rᵢ / p) + 1), rᵢ)
#         # TODO check math
#         q = p * (log(-r_i / p) + 1)
#         cut = JuMP.@expression(oa_model, p * λ[i] + q * θ + r_i * w[i])
#         push!(cuts, cut)
#     end
#     return cuts
# end
