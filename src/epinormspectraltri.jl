#=
epigraph of symmetric/Hermitian matrix spectral norm
(u, w) : u ≥ σ₁(W)

dual cone is epigraph of nuclear norm
(u, w) : u ≥ ∑ᵢ σᵢ(mat(w))

W = mat(w), in 𝕊(d) or ℍ(d)
eigendecomposition is W = V Diag(σ) V'
=#

mutable struct EpiNormSpectralTri{D <: PrimDual, C <: RealCompF} <: Cache
    oa_s::Vector{AE}
    d::Int
    w_temp::Vector{RealF}
    W_temp::Matrix{C}
    EpiNormSpectralTri{D, C}() where {D <: PrimDual, C <: RealCompF} = new{D, C}()
end

function Pajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiNormSpectralTriCone{RealF, C},
    ::Optimizer,
) where {C <: RealCompF}
    D = primal_or_dual(cone.use_dual)
    cache = EpiNormSpectralTri{D, C}()
    cache.oa_s = oa_s
    dim = MOI.dimension(cone)
    d = cache.d = svec_side(C, dim - 1)
    cache.w_temp = zeros(RealF, dim - 1)
    cache.W_temp = zeros(C, d, d)
    return cache
end

function get_eig(sz::Vector{RealF}, cache::EpiNormSpectralTri)
    @views W = svec_to_smat!(cache.W_temp, sz[2:end], rt2)
    return eigen(Hermitian(W, :U))
end

function Pajarito.Cones.add_init_cuts(cache::EpiNormSpectralTri, opt::Optimizer)
    # add variable bound
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, u >= 0)
    # TODO see comment in epinormspectral
    return
end

# primal cone functions

function Pajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::EpiNormSpectralTri{Prim},
    opt::Optimizer,
)
    F = get_eig(z, cache)
    cuts = AE[]
    for (i, σ_i) in enumerate(F.values)
        abs(σ_i) < 1e-8 && continue
        cut = _get_cut(σ_i, i, F.vectors, cache, opt)
        push!(cuts, cut)
    end
    return cuts
end

function Pajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::EpiNormSpectralTri{Prim},
    opt::Optimizer,
)
    # decomposed gradient cut is (1, -sign(σᵢ) Vᵢ Vᵢ')
    # TODO check math
    us = s[1]
    F = get_eig(s, cache)
    cuts = AE[]
    for (i, σ_i) in enumerate(F.values)
        if abs(σ_i) < opt.tol_feas || us - abs(σ_i) > -opt.tol_feas
            continue
        end
        cut = _get_cut(-sign(σ_i), i, F.vectors, cache, opt)
        push!(cuts, cut)
    end
    return cuts
end

function _get_cut(
    σ_i::RealF,
    i::Int,
    V::Matrix{C},
    cache::EpiNormSpectralTri{Prim, C},
    opt::Optimizer,
) where {C}
    u = cache.oa_s[1]
    w = cache.oa_s[2:end] # TODO cache
    R_i = cache.W_temp
    @views V_i = V[:, i]
    mul!(R_i, V_i, V_i', σ_i, false)
    R_vec_i = smat_to_svec!(cache.w_temp, R_i, rt2)
    return JuMP.@expression(opt.oa_model, abs(σ_i) * u + dot(R_vec_i, w))
end

# dual cone functions

function Pajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::EpiNormSpectralTri{Dual},
    opt::Optimizer,
)
    # strengthened cut is (‖R‖∞, R)
    # TODO extreme ray decomposition?
    F = get_eig(z, cache)
    p = maximum(abs, F.values)
    p < 1e-8 && return AE[]
    z2 = vcat(p, z[2:end])
    cut = dot_expr(z2, cache.oa_s, opt)
    return [cut]
end

function Pajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::EpiNormSpectralTri{Dual},
    opt::Optimizer,
)
    # gradient cut is (1, -∑ᵢ sign(σ_i) Vᵢ Vᵢ')
    # TODO check math
    F = get_eig(s, cache)
    if s[1] - sum(abs, F.values) > -opt.tol_feas
        return AE[]
    end
    R = F.vectors * Diagonal(-sign.(F.values)) * F.vectors'
    R_vec = smat_to_svec!(cache.w_temp, R, rt2)
    z2 = vcat(1, R_vec)
    cut = dot_expr(z2, cache.oa_s, opt)
    return [cut]
end
