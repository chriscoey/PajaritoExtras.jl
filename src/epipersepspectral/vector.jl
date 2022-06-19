#=
real vector domain ℝᵈ₊

extended formulation
∃ λ, u ≥ Σᵢ λᵢ, (λᵢ, v, wᵢ) ∈ EpiPerSepSpectralCone
=#

mutable struct VectorEpiPerSepSpectral{D <: PrimDual, E <: NatExt} <: Cache
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
    opt::Optimizer,
)
    D = primal_or_dual(use_dual)
    E = nat_or_ext(opt, d)
    return VectorEpiPerSepSpectral{D, E}()
end

function Pajarito.Cones.add_init_cuts(
    cache::VectorEpiPerSepSpectral{D},
    opt::Optimizer,
) where {D}
    (_, v) = swap_epiper(D, cache.oa_s[1:2]...)
    @views w = cache.oa_s[3:end]

    # add variable bounds
    JuMP.@constraint(opt.oa_model, v >= 0)
    if dom_pos(D, cache.h)
        JuMP.@constraint(opt.oa_model, w .>= 0)
    end
    opt.use_init_fixed_oa || return

    # add cuts at p = 1 and R = r₀ e
    r_vals = init_r_vals(D, cache.h)
    for r0 in r_vals
        r = fill(r0, cache.d)
        cuts = _get_cuts(1.0, r, cache, opt)
        JuMP.@constraint(opt.oa_model, cuts .>= 0)
    end
    return
end

function Pajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::VectorEpiPerSepSpectral{D},
    opt::Optimizer,
) where {D}
    (p, _) = swap_epiper(D, z[1:2]...)
    return _get_cuts(p, z[3:end], cache, opt)
end

function Pajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::VectorEpiPerSepSpectral{D},
    opt::Optimizer,
) where {D}
    h = cache.h
    ws = s[3:end]
    if dom_pos(D, h)
        ws = max.(ws, 1e-10)
    end

    (us, vs) = swap_epiper(D, s[1:2]...)
    v_pos = max(vs, 1e-9)
    if us - per_sepspec(val_or_conj(D), h, v_pos, ws) > -opt.tol_feas
        return AE[]
    end

    # gradient cut is (1, h⋆(r), r) at r = -h'(ws / vs)
    r = -h_grad(D, h, ws / v_pos)
    return _get_cuts(1.0, r, cache, opt)
end

# unextended formulation

function _get_cuts(
    p::RealF,
    r::Vector{RealF},
    cache::VectorEpiPerSepSpectral{D, Nat},
    opt::Optimizer,
) where {D}
    # strengthened cut is (p, p * h⋆(r / p), r)
    q = per_sepspec(conj_or_val(D), cache.h, p, r)
    (p, q) = swap_epiper(D, p, q)

    z = vcat(p, q, r)
    cut = dot_expr(z, cache.oa_s, opt)
    return [cut]
end

# extended formulation

function Pajarito.Cones.num_ext_variables(cache::VectorEpiPerSepSpectral{<:PrimDual, Ext})
    return cache.d
end

function Pajarito.Cones.extend_start(
    cache::VectorEpiPerSepSpectral{D, Ext},
    s_start::Vector{RealF},
    opt::Optimizer,
) where {D}
    (_, v_start) = swap_epiper(D, s_start[1:2]...)
    if v_start < 1e-8
        return zeros(cache.d)
    end
    w_start = [max(s_start[i], 1e-9) for i in 3:length(s_start)]
    return [per_sepspec(val_or_conj(D), cache.h, v_start, [w_i]) for w_i in w_start]
end

function Pajarito.Cones.setup_auxiliary(
    cache::VectorEpiPerSepSpectral{D, Ext},
    opt::Optimizer,
) where {D}
    @assert cache.d >= 2
    (u, _) = swap_epiper(D, cache.oa_s[1:2]...)
    λ = cache.λ = JuMP.@variable(opt.oa_model, [1:(cache.d)])
    JuMP.@constraint(opt.oa_model, u >= sum(λ))
    return λ
end

function _get_cuts(
    p::RealF,
    r::Vector{RealF},
    cache::VectorEpiPerSepSpectral{D, Ext},
    opt::Optimizer,
) where {D}
    (_, v) = swap_epiper(D, cache.oa_s[1:2]...)
    @views w = cache.oa_s[3:end]
    λ = cache.λ

    cuts = AE[]
    for i in 1:(cache.d)
        r_i = r[i]
        # r_i = min(r[i], 1e-7) # iszero(r_i) && continue # TODO?
        # strengthened disaggregated cut on (λᵢ, v, wᵢ) is (p, p * h⋆(rᵢ / p), rᵢ)
        # TODO check math
        q = per_sepspec(conj_or_val(D), cache.h, p, [r_i])
        cut = JuMP.@expression(opt.oa_model, p * λ[i] + q * v + r_i * w[i])
        push!(cuts, cut)
    end
    return cuts
end
