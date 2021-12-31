#=
real vector domain ℝᵈ₊

extended formulation
∃ λ, u ≥ Σᵢ λᵢ, (λᵢ, v, wᵢ) ∈ EpiPerSepSpectralCone
=#

mutable struct VectorEpiPerSepSpectral{D <: PrimDual, E <: NatExt} <: Cone
    oa_s::Vector{AE}
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
    E = nat_or_ext(extend, d)
    return VectorEpiPerSepSpectral{D, E}()
end

function MOIPajarito.Cones.add_init_cuts(
    cache::VectorEpiPerSepSpectral{D},
    oa_model::JuMP.Model,
) where {D}
    (_, v) = swap_epiper(D, cache.oa_s[1:2]...)
    @views w = cache.oa_s[3:end]

    # variable bounds
    JuMP.@constraint(oa_model, v >= 0)
    if dom_pos(D, cache.h)
        JuMP.@constraint(oa_model, w .>= 0)
    end
    num_cuts = 1 + cache.d

    # cuts using values of p = 1 and r = r₀ e
    r_vals = init_r_vals(D, cache.h)
    for r0 in r_vals
        r = fill(r0, cache.d)
        cuts = _get_cuts(1.0, r, cache, oa_model)
        JuMP.@constraint(oa_model, cuts .>= 0)
        num_cuts += length(cuts)
    end
    return num_cuts
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::VectorEpiPerSepSpectral{D},
    oa_model::JuMP.Model,
) where {D}
    (p, _) = swap_epiper(D, z[1:2]...)
    return _get_cuts(p, z[3:end], cache, oa_model)
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::VectorEpiPerSepSpectral{D},
    oa_model::JuMP.Model,
) where {D}
    h = cache.h
    ws = s[3:end]
    if dom_pos(D, h)
        ws = max.(ws, 1e-7)
    end

    (us, vs) = swap_epiper(D, s[1:2]...)
    v_pos = max(vs, 1e-7)
    if us - per_sepspec(val_or_conj(D), h, v_pos, ws) > -1e-7
        return AE[]
    end

    # gradient cut is (1, h⋆(r), r) at r = -h'(ws / vs)
    r = -h_grad(D, h, ws / v_pos)
    return _get_cuts(1.0, r, cache, oa_model)
end

# primal unextended formulation

function _get_cuts(
    p::RealF,
    r::Vector{RealF},
    cache::VectorEpiPerSepSpectral{D, Nat},
    oa_model::JuMP.Model,
) where {D}
    # strengthened cut is (p, p * h⋆(r / p), r)
    q = per_sepspec(conj_or_val(D), cache.h, p, r)
    (p, q) = swap_epiper(D, p, q)

    z = vcat(p, q, r)
    cut = dot_expr(z, cache.oa_s, oa_model)
    return [cut]
end

# primal extended formulation

function MOIPajarito.Cones.num_ext_variables(
    cache::VectorEpiPerSepSpectral{<:PrimDual, Ext},
)
    return cache.d
end

function MOIPajarito.Cones.extend_start(
    cache::VectorEpiPerSepSpectral{D, Ext},
    s_start::Vector{RealF},
) where {D}
    (_, v_start) = swap_epiper(D, s_start[1:2]...)
    @views w_start = s_start[3:end]
    return [per_sepspec(val_or_conj(D), cache.h, v_start, [w_i]) for w_i in w_start]
end

function MOIPajarito.Cones.setup_auxiliary(
    cache::VectorEpiPerSepSpectral{D, Ext},
    oa_model::JuMP.Model,
) where {D}
    @assert cache.d >= 2
    (u, _) = swap_epiper(D, cache.oa_s[1:2]...)
    λ = cache.λ = JuMP.@variable(oa_model, [1:(cache.d)])
    JuMP.@constraint(oa_model, u >= sum(λ))
    return λ
end

function _get_cuts(
    p::RealF,
    r::Vector{RealF},
    cache::VectorEpiPerSepSpectral{D, Ext},
    oa_model::JuMP.Model,
) where {D}
    (_, v) = swap_epiper(D, cache.oa_s[1:2]...)
    @views w = cache.oa_s[3:end]
    λ = cache.λ

    cuts = AE[]
    for i in 1:(cache.d)
        r_i = r[i]
        # r_i = min(r[i], 1e-7) # iszero(r_i) && continue
        # strengthened disaggregated cut on (λᵢ, v, wᵢ) is (p, p * h⋆(rᵢ / p), rᵢ)
        # TODO check math
        q = per_sepspec(conj_or_val(D), cache.h, p, [r_i])
        cut = JuMP.@expression(oa_model, p * λ[i] + q * v + r_i * w[i])
        push!(cuts, cut)
    end
    return cuts
end
