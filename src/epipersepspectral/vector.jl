#=
real vector domain ℝᵈ₊

extended formulation
∃ λ, u ≥ Σᵢ λᵢ, (λᵢ, v, wᵢ) ∈ EpiPerSepSpectralCone
=#

mutable struct VectorEpiPerSepSpectral{D <: PrimDual, E <: NatExt} <: Cone
    oa_s::Vector{AE}
    s::Vector{RealF}
    h::SepSpectralFun
    d::Int
    λ::Vector{VR}
    function VectorEpiPerSepSpectral{D, E}() where {D <: PrimDual, E <: NatExt}
        return new{D, E}()
    end
end

function create_sepspectral_cache(
    ::Type{VectorCSqr{RealF}},
    use_dual::Bool,
    d::Int,
    extend::Bool,
)
    D = primal_or_dual(use_dual)
    E = extender(extend, d)
    return VectorEpiPerSepSpectral{D, E}()
end

# primal cone functions

function MOIPajarito.Cones.add_init_cuts(
    cache::VectorEpiPerSepSpectral{Primal},
    oa_model::JuMP.Model,
)
    v = cache.oa_s[2]
    @views w = cache.oa_s[3:end]
    d = cache.d
    # variable bounds v ≥ 0, wᵢ ≥ 0
    JuMP.@constraints(oa_model, begin
        v >= 0
        [i in 1:d], w[i] >= 0
    end)
    num_cuts = 1 + d

    # cuts using values of p = 1 and r = r₀ e
    r_vals = init_r_vals(cache.h)
    for r0 in r_vals
        r = fill(r0, d)
        cuts = _get_cuts(1.0, r, cache, oa_model)
        JuMP.@constraint(oa_model, cuts .>= 0)
        num_cuts += length(cuts)
    end
    return num_cuts
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::VectorEpiPerSepSpectral{Primal},
    oa_model::JuMP.Model,
)
    return _get_cuts(z[1], z[3:end], cache, oa_model)
end

function MOIPajarito.Cones.get_sep_cuts(
    cache::VectorEpiPerSepSpectral{Primal},
    oa_model::JuMP.Model,
)
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
    r = -h_grad(w_pos / v_pos, cache.h)
    return _get_cuts(1.0, r, cache, oa_model)
end

# primal unextended formulation

function _get_cuts(
    p::RealF,
    r::Vector{RealF},
    cache::VectorEpiPerSepSpectral{Primal, Nat},
    oa_model::JuMP.Model,
)
    # strengthened cut is (p, p * h⋆(r / p), r)
    @assert p > 1e-12
    if h_conj_dom_pos(cache.h)
        @assert all(>(1e-12), r)
    end
    q = per_sepspec(h_conj, cache.h, p, r)
    z = vcat(p, q, r)
    cut = dot_expr(z, cache.oa_s, oa_model)
    return [cut]
end

# primal extended formulation

function MOIPajarito.Cones.num_ext_variables(cache::VectorEpiPerSepSpectral{Primal, Ext})
    return cache.d
end

function MOIPajarito.Cones.extend_start(
    cache::VectorEpiPerSepSpectral{Primal, Ext},
    s_start::Vector{RealF},
)
    v_start = s_start[2]
    w_start = s_start[3:end]
    # @assert s_start[1] - per_sepspec(h_val, cache.h, v_start, w_start) >= -1e-7
    return [per_sepspec(h_val, cache.h, v_start, [w_i]) for w_i in w_start]
end

function MOIPajarito.Cones.setup_auxiliary(
    cache::VectorEpiPerSepSpectral{Primal, Ext},
    oa_model::JuMP.Model,
)
    @assert cache.d >= 2
    u = cache.oa_s[1]
    λ = cache.λ = JuMP.@variable(oa_model, [1:(cache.d)])
    JuMP.@constraint(oa_model, u >= sum(λ))
    return λ
end

function _get_cuts(
    p::RealF,
    r::Vector{RealF},
    cache::VectorEpiPerSepSpectral{Primal, Ext},
    oa_model::JuMP.Model,
)
    @assert p > 1e-12
    if h_conj_dom_pos(cache.h)
        @assert all(>(1e-12), r)
    end
    v = cache.oa_s[2]
    @views w = cache.oa_s[3:end]
    λ = cache.λ
    cuts = AE[]
    for i in 1:(cache.d)
        r_i = r[i]
        # r_i = min(r[i], 1e-7) # iszero(r_i) && continue
        # strengthened disaggregated cut on (λᵢ, v, wᵢ) is (p, p * h⋆(rᵢ / p), rᵢ)
        # TODO check math
        q = per_sepspec(h_conj, cache.h, p, [r_i])
        cut = JuMP.@expression(oa_model, p * λ[i] + q * v + r_i * w[i])
        push!(cuts, cut)
    end
    return cuts
end

# dual cone functions
