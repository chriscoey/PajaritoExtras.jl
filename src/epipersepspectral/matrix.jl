# real symmetric or complex Hermitian (svec scaled triangle) domain 𝕊ᵈ₊

mutable struct MatrixEpiPerSepSpectralCache{C <: RealOrComplex} <: ConeCache
    cone::Hypatia.EpiPerSepSpectralCone{Float64}
    is_complex::Bool
    oa_s::Vector{AE}
    s::Vector{Float64}
    d::Int
    w_temp::Vector{Float64}
    W_temp::Matrix{C}
    MatrixEpiPerSepSpectralCache{C}() where {C <: RealOrComplex} = new{C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiPerSepSpectralCone{Float64},
    ::Bool,
) where {C <: RealOrComplex}
    cache = MatrixEpiPerSepSpectralCache{C}()
    cache.cone = cone
    cache.is_complex = (C == ComplexF64)
    cache.oa_s = oa_s
    dim = MOI.dimension(cone)
    d = cache.d = svec_side(C, dim - 1)
    cache.w_temp = zeros(Float64, dim - 1)
    cache.W_temp = zeros(C, d, d)
    return cache
end

function MOIPajarito.Cones.add_init_cuts(
    cache::MatrixEpiPerSepSpectralCache{C},
    oa_model::JuMP.Model,
) where {C}
    d = cache.d
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    # W_ii ≥ 0
    # linearize at W = I, cut on (u, W) is (-d, I)
    W_diag = [w[svec_idx(C, i, i)] for i in 1:d]
    JuMP.@constraints(oa_model, begin
        W_diag .>= 0
        -d * u + sum(W_diag) >= 0
    end)
    return d + 1
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::MatrixEpiPerSepSpectralCache,
    oa_model::JuMP.Model,
)
    p = z[1]
    R = cache.W_temp
    @views svec_to_smat!(R, z[2:end], rt2)
    F = eigen(Hermitian(R, :U))
    V = F.vectors
    ω = F.values
    @assert all(>(-1e-7), ω)
    @assert geomean(ω) - p > -1e-7

    cuts = AE[]
    if p > -1e-7
        # add eigenvector cuts
        num_pos = count(>(1e-7), ω)
        @views V_neg = V[:, 1:num_pos] * Diagonal(sqrt.(ω[1:num_pos]))
        cuts = _get_psd_cuts(V_neg, cache, oa_model)
    end

    # add rootdet cut
    ω_pos = max.(ω, 1e-7)
    cut = _get_rootdet_cut(ω_pos, V, cache, oa_model)
    push!(cuts, cut)
    return cuts
end

function MOIPajarito.Cones.get_sep_cuts(
    cache::MatrixEpiPerSepSpectralCache,
    oa_model::JuMP.Model,
)
    # check s ∉ K
    us = cache.s[1]
    Ws = cache.W_temp
    @views svec_to_smat!(Ws, cache.s[2:end], rt2)
    F = eigen(Hermitian(Ws, :U))
    V = F.vectors
    ω = F.values
    num_neg = count(<(-1e-7), ω)
    iszero(num_neg) && us < 1e-7 && return AE[]
    @assert issorted(ω)

    cuts = AE[]
    if !iszero(num_neg)
        # add eigenvector cuts
        V_neg = V[:, 1:num_neg]
        cuts = _get_psd_cuts(V_neg, cache, oa_model)
    end

    geomean(ω) - us > -1e-7 && return AE[]
    # gradient cut is (-1, V * Diagonal(geom(ω) / d ./ ω) * V')
    ω_pos = max.(ω, 1e-7)
    rω = (geomean(ω_pos) / cache.d) ./ ω_pos
    cut = _get_rootdet_cut(rω, V, cache, oa_model)
    push!(cuts, cut)
    return cuts
end

function _get_rootdet_cut(
    rω::Vector{Float64},
    V::Matrix{C},
    cache::MatrixEpiPerSepSpectralCache{C},
    oa_model::JuMP.Model,
) where {C}
    # strengthened cut is (-d * geom(rω), V * Diagonal(rω) * V')
    p = -cache.d * geomean(rω)
    R = V * Diagonal(rω) * V'
    r = smat_to_svec!(cache.w_temp, R, rt2)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    return JuMP.@expression(oa_model, p * u + JuMP.dot(r, w))
end
