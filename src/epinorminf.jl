#=
epigraph of ∞-norm of real or complex vectors
(u, w) : u ≥ ‖Wᵢ‖, ∀i
no extended formulation needed
real case (polyhedral): u ≥ Wᵢ, u ≥ -Wᵢ
complex case: (u, re(Wᵢ), im(Wᵢ)) in EpiNormEucl

dual cone is epigraph of 1-norm
(u, w) : u ≥ ∑ᵢ ‖Wᵢ‖
extended formulation: ∃ λ, u ≥ Σᵢ λᵢ
real case (polyhedral): λᵢ ≥ Wᵢ, λᵢ ≥ -Wᵢ
complex case: (λᵢ, re(Wᵢ), im(Wᵢ)) in EpiNormEucl
=#

mutable struct EpiNormInf{D <: PrimDual, C <: RealCompF, E <: NatExt} <: Cache
    oa_s::Vector{AE}
    d::Int
    λ::Vector{VR}
    W_temp::Vector{C}
    function EpiNormInf{D, C, E}() where {D <: PrimDual, C <: RealCompF, E <: NatExt}
        return new{D, C, E}()
    end
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiNormInfCone{RealF, C},
    opt::Optimizer,
) where {C <: RealCompF}
    dim = MOI.dimension(cone)
    @assert dim == length(oa_s)
    d = (C == CompF ? div(dim - 1, 2) : dim - 1)
    E = (cone.use_dual ? nat_or_ext(opt, d) : Nat)
    D = primal_or_dual(cone.use_dual)
    cache = EpiNormInf{D, C, E}()
    cache.oa_s = oa_s
    cache.d = d
    cache.W_temp = zeros(C, d)
    return cache
end

# primal cone functions

function MOIPajarito.Cones.add_init_cuts(cache::EpiNormInf{Prim, RealF}, opt::Optimizer)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    # polyhedral (no OA)
    JuMP.@constraints(opt.oa_model, begin
        u .>= w
        u .>= -w
    end)
    return 2 * cache.d
end

function MOIPajarito.Cones.add_init_cuts(cache::EpiNormInf{Prim, CompF}, opt::Optimizer)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    d = cache.d
    # cuts on (u, re(Wᵢ), im(Wᵢ)) are (rt2, ±1, ±1) ∈ EpiNormEucl
    for i in 1:d
        re_w = w[2i - 1]
        im_w = w[2i]
        JuMP.@constraints(opt.oa_model, begin
            rt2 * u - re_w - im_w >= 0
            rt2 * u - re_w + im_w >= 0
            rt2 * u + re_w - im_w >= 0
            rt2 * u + re_w + im_w >= 0
        end)
    end
    return 4d
end

function MOIPajarito.Cones.get_subp_cuts(
    ::Vector{RealF},
    ::EpiNormInf{Prim, RealF},
    ::Optimizer,
)
    return AE[]
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::EpiNormInf{Prim, CompF},
    opt::Optimizer,
)
    # strengthened cut on (u, re(Wᵢ), im(Wᵢ)) is (‖Rᵢ‖, re(Rᵢ), im(Rᵢ))
    R = cache.W_temp
    @views vec_copyto!(R, z[2:end])
    cuts = AE[]
    for (i, R_i) in enumerate(R)
        abs(R_i) < 1e-7 && continue
        cut = _get_cut(R_i, i, cache, opt)
        push!(cuts, cut)
    end
    return cuts
end

function MOIPajarito.Cones.get_sep_cuts(
    ::Vector{RealF},
    ::EpiNormInf{Prim, RealF},
    ::Optimizer,
)
    return AE[]
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::EpiNormInf{Prim, CompF},
    opt::Optimizer,
)
    # decomposed gradient cut on (u, Wᵢ) is (1, -Wsᵢ / ‖Wsᵢ‖)
    us = s[1]
    @views Ws = vec_copyto!(cache.W_temp, s[2:end])
    cuts = AE[]
    for (i, Ws_i) in enumerate(Ws)
        abs_i = abs(Ws_i)
        if abs_i < 1e-7 || us - abs_i > -1e-7
            continue
        end
        R_i = Ws_i / -abs_i
        cut = _get_cut(R_i, i, cache, opt)
        push!(cuts, cut)
    end
    return cuts
end

function _get_cut(R_i::CompF, i::Int, cache::EpiNormInf{Prim, CompF}, opt::Optimizer)
    u = cache.oa_s[1]
    re_w = cache.oa_s[2i]
    im_w = cache.oa_s[2i + 1]
    return JuMP.@expression(
        opt.oa_model,
        abs(R_i) * u + real(R_i) * re_w + imag(R_i) * im_w
    )
end

# dual cone functions

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::EpiNormInf{Dual},
    opt::Optimizer,
)
    R = cache.W_temp
    @views vec_copyto!(R, z[2:end])
    return _get_cuts(R, cache, opt)
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::EpiNormInf{Dual},
    opt::Optimizer,
)
    us = s[1]
    Ws = cache.W_temp
    @views vec_copyto!(Ws, s[2:end])
    if us - norm(Ws, 1) > -1e-7
        return AE[]
    end

    # gradient cut is (1, (-Wsᵢ / ‖Wsᵢ‖)ᵢ)
    R = zero(Ws)
    for (i, Ws_i) in enumerate(Ws)
        if abs(Ws_i) > 1e-12
            R[i] = -Ws_i / abs(Ws_i)
        end
    end
    return _get_cuts(R, cache, opt)
end

# unextended formulation

function MOIPajarito.Cones.add_init_cuts(
    cache::EpiNormInf{Dual, RealF, Nat},
    opt::Optimizer,
)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    # polyhedral but needs exponentially many constraints
    JuMP.@constraints(opt.oa_model, begin
        u >= 0
        u >= sum(w)
        u >= -sum(w)
        [w_i in w], u >= w_i
        [w_i in w], u >= -w_i
    end)
    return 3 + 2 * cache.d
end

function MOIPajarito.Cones.add_init_cuts(
    cache::EpiNormInf{Dual, CompF, Nat},
    opt::Optimizer,
)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    d = cache.d
    JuMP.@constraints(opt.oa_model, begin
        u >= 0
        u >= irt2 * sum(w)
        u >= -irt2 * sum(w)
        [i in 1:d], u >= irt2 * (w[2i - 1] + w[2i])
        [i in 1:d], u >= -irt2 * (w[2i - 1] + w[2i])
        [i in 1:d], u >= irt2 * (w[2i - 1] - w[2i])
        [i in 1:d], u >= -irt2 * (w[2i - 1] - w[2i])
    end)
    return 3 + 4 * cache.d
end

function _get_cuts(R::Vector{C}, cache::EpiNormInf{Dual, C, Nat}, opt::Optimizer) where {C}
    # strengthened cut on (u, W) is (‖R‖∞, R)
    # TODO extreme ray decomposition?
    p = norm(R, Inf)
    r = reinterpret(RealF, R)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    cut = JuMP.@expression(opt.oa_model, p * u + JuMP.dot(r, w))
    return [cut]
end

# extended formulation

function MOIPajarito.Cones.num_ext_variables(cache::EpiNormInf{Dual, <:RealCompF, Ext})
    return cache.d
end

function MOIPajarito.Cones.extend_start(
    cache::EpiNormInf{Dual, C, Ext},
    s_start::Vector{RealF},
) where {C}
    s_start[1] < 1e-7 && return zeros(cache.d)
    w_start = s_start[2:end]
    λ_start = abs.(reinterpret(C, w_start))
    return λ_start
end

function MOIPajarito.Cones.setup_auxiliary(
    cache::EpiNormInf{Dual, <:RealCompF, Ext},
    opt::Optimizer,
)
    @assert cache.d >= 2
    u = cache.oa_s[1]
    λ = cache.λ = JuMP.@variable(opt.oa_model, [1:(cache.d)], lower_bound = 0)
    JuMP.@constraint(opt.oa_model, u >= sum(λ))
    return λ
end

function MOIPajarito.Cones.add_init_cuts(
    cache::EpiNormInf{Dual, RealF, Ext},
    opt::Optimizer,
)
    @views w = cache.oa_s[2:end]
    λ = cache.λ
    d = cache.d
    # polyhedral (no OA)
    JuMP.@constraints(opt.oa_model, begin
        [i in 1:d], λ[i] >= w[i]
        [i in 1:d], λ[i] >= -w[i]
    end)
    return 2d
end

function MOIPajarito.Cones.add_init_cuts(
    cache::EpiNormInf{Dual, CompF, Ext},
    opt::Optimizer,
)
    @views w = cache.oa_s[2:end]
    λ = cache.λ
    d = cache.d
    # cuts on (λᵢ, re(Wᵢ), im(Wᵢ)) are (rt2, ±1, ±1) ∈ EpiNormEucl
    for i in 1:d
        λ_i = λ[i]
        re_w = w[2i - 1]
        im_w = w[2i]
        JuMP.@constraints(opt.oa_model, begin
            rt2 * λ_i - re_w - im_w >= 0
            rt2 * λ_i - re_w + im_w >= 0
            rt2 * λ_i + re_w - im_w >= 0
            rt2 * λ_i + re_w + im_w >= 0
        end)
    end
    return 4d
end

function MOIPajarito.Cones.get_subp_cuts(
    ::Vector{RealF},
    ::EpiNormInf{Dual, RealF, Ext},
    ::Optimizer,
)
    return AE[]
end

function MOIPajarito.Cones.get_sep_cuts(
    ::Vector{RealF},
    ::EpiNormInf{Dual, RealF, Ext},
    ::Optimizer,
)
    return AE[]
end

function _get_cuts(R::Vector{CompF}, cache::EpiNormInf{Dual, CompF, Ext}, opt::Optimizer)
    cuts = AE[]
    for (i, R_i) in enumerate(R)
        # strengthened disaggregated cut on (λᵢ, re(Wᵢ), im(Wᵢ)) is (‖Rᵢ‖, re(Rᵢ), im(Rᵢ))
        abs(R_i) < 1e-7 && continue
        cut = _get_cut(R_i, i, cache, opt)
        push!(cuts, cut)
    end
    return cuts
end

function _get_cut(R_i::CompF, i::Int, cache::EpiNormInf{Dual, CompF, Ext}, opt::Optimizer)
    λ_i = cache.λ[i]
    re_w = cache.oa_s[2i]
    im_w = cache.oa_s[2i + 1]
    return JuMP.@expression(
        opt.oa_model,
        abs(R_i) * λ_i + real(R_i) * re_w + imag(R_i) * im_w
    )
end
