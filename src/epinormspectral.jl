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

mutable struct EpiNormSpectral{D <: PrimDual, C <: RealComp} <: Cone
    oa_s::Vector{AE}
    s::Vector{Float64}
    d1::Int
    d2::Int
    w_temp::Vector{Float64}
    W_temp::Matrix{C}
    EpiNormSpectral{D, C}() where {D <: PrimDual, C <: RealComp} = new{D, C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiNormSpectralCone{Float64, C},
    ::Bool,
) where {C <: RealComp}
    D = primal_or_dual(cone.use_dual)
    cache = EpiNormSpectral{D, C}()
    cache.oa_s = oa_s
    d1 = cache.d1 = cone.d1
    d2 = cache.d2 = cone.d2
    dim = MOI.dimension(cone)
    @assert dim == 1 + vec_length(C, d1 * d2)
    cache.w_temp = zeros(Float64, dim - 1)
    cache.W_temp = zeros(C, d1, d2)
    return cache
end

function get_svd(sz::Vector{Float64}, cache::EpiNormSpectral)
    @views W = vec_copyto!(cache.W_temp, sz[2:end])
    return svd(W, full = false)
end

function MOIPajarito.Cones.add_init_cuts(cache::EpiNormSpectral, oa_model::JuMP.Model)
    # TODO use simple bounds to derive init cuts:
    # frob / rtd1 <= spec
    # opinf / rtd2 <= spec
    # op1 / rtd1 <= spec
    # and
    # frob <= nuc
    # u ≥ 0
    u = cache.oa_s[1]
    JuMP.@constraint(oa_model, u >= 0)
    return 1
end

# primal cone functions

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::EpiNormSpectral{Primal},
    oa_model::JuMP.Model,
)
    F = get_svd(z, cache)
    cuts = AE[]
    for (i, σ_i) in enumerate(F.S)
        σ_i < 1e-7 && break
        cut = _get_cut(σ_i, i, F.U, F.Vt, cache, oa_model)
        push!(cuts, cut)
    end
    return cuts
end

function MOIPajarito.Cones.get_sep_cuts(
    cache::EpiNormSpectral{Primal},
    oa_model::JuMP.Model,
)
    # decomposed gradient cut is (1, -Uᵢ * Vtᵢ)
    us = cache.s[1]
    F = get_svd(cache.s, cache)
    cuts = AE[]
    for (i, σ_i) in enumerate(F.S)
        if σ_i < 1e-7 || us - σ_i > -1e-7
            break
        end
        cut = _get_cut(one(T), i, F.U, F.Vt, cache, oa_model)
        push!(cuts, cut)
    end
    return cuts
end

function _get_cut(
    σ_i::Float64,
    i::Int,
    U::Matrix{C},
    Vt::Matrix{C},
    cache::EpiNormSpectral{Primal, C},
    oa_model::JuMP.Model,
) where {C}
    u = cache.oa_s[1]
    w = cache.oa_s[2:end] # TODO cache
    R_i = cache.W_temp
    @views mul!(R_i, U[:, i], transpose(Vt[i, :]), σ_i, false)
    R_vec_i = vec_copyto!(cache.w_temp, R_i) # TODO maybe reinterpret
    return JuMP.@expression(oa_model, σ_i * u + JuMP.dot(R_vec_i, w))
end

# dual cone functions

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::EpiNormSpectral{Dual},
    oa_model::JuMP.Model,
)
    # strengthened cut is (‖R‖∞, R)
    # TODO extreme ray decomposition?
    F = get_svd(z, cache)
    p = F.S[1]
    p < 1e-7 && return AE[]
    z2 = vcat(p, z[2:end])
    cut = dot_expr(z2, cache.oa_s, oa_model)
    return [cut]
end

function MOIPajarito.Cones.get_sep_cuts(cache::EpiNormSpectral{Dual}, oa_model::JuMP.Model)
    # gradient cut is (1, -∑ᵢ Uᵢ * Vtᵢ)
    # TODO check math
    F = get_svd(cache.s, cache)
    (cache.s[1] - sum(F.S) > -1e-7) && return AE[]
    R = -F.U * F.Vt
    R_vec = vec_copyto!(cache.w_temp, R) # TODO maybe reinterpret
    z2 = vcat(1, R_vec)
    cut = dot_expr(z2, cache.oa_s, oa_model)
    return [cut]
end
