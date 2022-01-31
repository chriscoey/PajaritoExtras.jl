#=
epigraph of matrix spectral norm
(u, w) : u ≥ σ₁(W)
TODO equivalent to u ≥ 0, u I - W W' / u ⪰ 0, so from eigendecomposition we can get SOC cuts

dual cone is epigraph of matrix nuclear norm
(u, w) : u ≥ ∑ᵢ σᵢ(mat(w))

W = mat(w), in ℝ(d1, d2) or ℂ(d1, d2), where d1 ≤ d2
SVD is W = U Diag(σ) Vt
subdifferential characterized by e.g.
G.A. Watson, "Characterization of the Subdifferential of Some Matrix Norms"
=#

mutable struct EpiNormSpectral{D <: PrimDual, C <: RealCompF} <: Cache
    oa_s::Vector{AE}
    d1::Int
    d2::Int
    w_temp::Vector{RealF}
    W_temp::Matrix{C}
    EpiNormSpectral{D, C}() where {D <: PrimDual, C <: RealCompF} = new{D, C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiNormSpectralCone{RealF, C},
    ::Optimizer,
) where {C <: RealCompF}
    D = primal_or_dual(cone.use_dual)
    cache = EpiNormSpectral{D, C}()
    cache.oa_s = oa_s
    d1 = cache.d1 = cone.d1
    d2 = cache.d2 = cone.d2
    dim = MOI.dimension(cone)
    @assert dim == 1 + vec_length(C, d1 * d2)
    cache.w_temp = zeros(RealF, dim - 1)
    cache.W_temp = zeros(C, d1, d2)
    return cache
end

function get_svd(sz::Vector{RealF}, cache::EpiNormSpectral)
    @views W = vec_copyto!(cache.W_temp, sz[2:end])
    return svd(W, full = false)
end

function MOIPajarito.Cones.add_init_cuts(cache::EpiNormSpectral, opt::Optimizer)
    # TODO use simple bounds to derive init cuts:
    # frob / rtd1 <= spec
    # opinf / rtd2 <= spec
    # op1 / rtd1 <= spec
    # and
    # frob <= nuc
    # u ≥ 0
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, u >= 0)
    return 1
end

# primal cone functions

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::EpiNormSpectral{Prim},
    opt::Optimizer,
)
    F = get_svd(z, cache)
    cuts = AE[]
    for (i, σ_i) in enumerate(F.S)
        σ_i < 1e-7 && break
        cut = _get_cut(σ_i, i, F.U, F.Vt, cache, opt)
        push!(cuts, cut)
    end
    return cuts
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::EpiNormSpectral{Prim},
    opt::Optimizer,
)
    # decomposed gradient cut is (1, -Uᵢ * Vtᵢ)
    us = s[1]
    F = get_svd(s, cache)
    cuts = AE[]
    for (i, σ_i) in enumerate(F.S)
        if σ_i < 1e-7 || us - σ_i > -1e-7
            break
        end
        cut = _get_cut(1.0, i, F.U, F.Vt, cache, opt)
        push!(cuts, cut)
    end
    return cuts
end

function _get_cut(
    σ_i::RealF,
    i::Int,
    U::Matrix{C},
    Vt::Matrix{C},
    cache::EpiNormSpectral{Prim, C},
    opt::Optimizer,
) where {C}
    u = cache.oa_s[1]
    w = cache.oa_s[2:end] # TODO cache
    R_i = cache.W_temp
    @views mul!(R_i, U[:, i], transpose(Vt[i, :]), σ_i, false)
    R_vec_i = vec_copyto!(cache.w_temp, R_i) # TODO maybe reinterpret
    return JuMP.@expression(opt.oa_model, σ_i * u + JuMP.dot(R_vec_i, w))
end

# dual cone functions

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::EpiNormSpectral{Dual},
    opt::Optimizer,
)
    # strengthened cut is (‖R‖∞, R)
    # TODO extreme ray decomposition?
    F = get_svd(z, cache)
    p = F.S[1]
    p < 1e-7 && return AE[]
    z2 = vcat(p, z[2:end])
    cut = dot_expr(z2, cache.oa_s, opt)
    return [cut]
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::EpiNormSpectral{Dual},
    opt::Optimizer,
)
    # gradient cut is (1, -∑ᵢ Uᵢ * Vtᵢ)
    # TODO check math
    F = get_svd(s, cache)
    (s[1] - sum(F.S) > -1e-7) && return AE[]
    R = -F.U * F.Vt
    R_vec = vec_copyto!(cache.w_temp, R) # TODO maybe reinterpret
    z2 = vcat(1, R_vec)
    cut = dot_expr(z2, cache.oa_s, opt)
    return [cut]
end
