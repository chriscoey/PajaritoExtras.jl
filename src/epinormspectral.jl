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

mutable struct EpiNormSpectralCache{C <: RealOrComplex} <: ConeCache
    cone::Hypatia.EpiNormSpectralCone{Float64, C}
    is_complex::Bool
    oa_s::Vector{AE}
    s::Vector{Float64}
    d1::Int
    d2::Int
    w_temp::Vector{Float64}
    W_temp::Matrix{C}
    EpiNormSpectralCache{C}() where {C <: RealOrComplex} = new{C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiNormSpectralCone{Float64, C},
    ::Bool,
) where {C <: RealOrComplex}
    @assert !cone.use_dual # TODO separate cache type DualEpiNormSpectralCache?
    cache = EpiNormSpectralCache{C}()
    cache.cone = cone
    cache.is_complex = (C == ComplexF64)
    cache.oa_s = oa_s
    d1 = cache.d1 = cone.d1
    d2 = cache.d2 = cone.d2
    dim = MOI.dimension(cone)
    @assert dim == 1 + vec_length(C, d1 * d2)
    cache.w_temp = zeros(Float64, dim - 1)
    cache.W_temp = zeros(C, d1, d2)
    return cache
end

function MOIPajarito.Cones.add_init_cuts(
    cache::EpiNormSpectralCache{C},
    oa_model::JuMP.Model,
) where {C}
    # TODO use simple bounds to derive init cuts: (could even add a SOC cut if supported, for the frob bound)
    # frob / rtd1 <= spec
    # opinf / rtd2 <= spec
    # op1 / rtd1 <= spec
    u = cache.oa_s[1]
    # w = cache.oa_s[2:end] # TODO cache
    # u ≥ 0
    JuMP.@constraint(oa_model, u >= 0)
    return 1
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::EpiNormSpectralCache,
    oa_model::JuMP.Model,
)
    # strengthened cuts from SVD are σᵢ * Uᵢ * Vᵢ', ∀ i : σᵢ = σ₁
    p = z[1]
    R = cache.W_temp
    @views vec_copyto!(R, z[2:end])
    F = svd!(R, full = false)
    σ = F.S
    @assert issorted(σ, rev = true)
    σ[1] <= 1e-7 && return AE[]
    @assert p - σ[1] >= -1e-7
    return _get_cuts(σ, F.U, F.Vt, cache, oa_model)
end

function MOIPajarito.Cones.get_sep_cuts(cache::EpiNormSpectralCache, oa_model::JuMP.Model)
    # check s ∉ K
    # NOTE only need to compute singular values larger than su, but not possible with LAPACK
    su = cache.s[1]
    sW = cache.W_temp
    @views vec_copyto!(sW, cache.s[2:end])
    F = svd!(sW, full = false)
    σ = F.S
    @assert issorted(σ, rev = true)
    if σ[1] <= 1e-7 || su - σ[1] >= -1e-7
        return AE[]
    end
    σ_viol = ones(count(>(su + 1e-7), σ))
    return _get_cuts(σ_viol, F.U, F.Vt, cache, oa_model)
end

function _get_cuts(
    σ::Vector{Float64},
    U::Matrix{C},
    Vt::Matrix{C},
    cache::EpiNormSpectralCache{C},
    oa_model::JuMP.Model,
) where {C}
    # cuts from svd are (σᵢ, σᵢ * Uᵢ * Vtᵢ)
    u = cache.oa_s[1]
    w = cache.oa_s[2:end] # TODO cache
    cuts = AE[]
    R_vec_i = cache.w_temp
    R_i = cache.W_temp
    for (i, σ_i) in enumerate(σ)
        σ_i < 1e-9 && continue
        @views mul!(R_i, U[:, i], transpose(Vt[i, :]), σ_i, false)
        clean_array!(R_i) && continue
        vec_copyto!(R_vec_i, R_i)
        cut = JuMP.@expression(oa_model, σ_i * u + JuMP.dot(R_vec_i, w))
        push!(cuts, cut)
    end
    return cuts
end
