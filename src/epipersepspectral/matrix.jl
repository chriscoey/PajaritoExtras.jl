#=
real symmetric or complex Hermitian (svec scaled triangle) domain 𝕊ᵈ₊
=#

mutable struct MatrixEpiPerSepSpectralCache{D <: PrimalOrDual, C <: RealOrComplex} <:
               ConeCache
    oa_s::Vector{AE}
    s::Vector{Float64}
    h::SepSpectralFun
    d::Int
    w_temp::Vector{Float64}
    W_temp::Matrix{C}
    function MatrixEpiPerSepSpectralCache{
        D,
        C,
    }() where {D <: PrimalOrDual, C <: RealOrComplex}
        return new{D, C}()
    end
end

function create_sepspectral_cache(
    ::Type{MatrixCSqr{Float64, C}},
    use_dual::Bool,
    d::Int,
    ::Bool,
) where {C <: RealOrComplex}
    D = primal_or_dual(use_dual)
    cache = MatrixEpiPerSepSpectralCache{D, C}()
    cache.w_temp = zeros(Float64, svec_length(C, d))
    cache.W_temp = zeros(C, d, d)
    return cache
end

# primal cone functions

function MOIPajarito.Cones.add_init_cuts(
    cache::MatrixEpiPerSepSpectralCache{Primal, C},
    oa_model::JuMP.Model,
) where {C}
    v = cache.oa_s[2]
    @views w = cache.oa_s[3:end]
    d = cache.d
    # variable bounds v ≥ 0, Wᵢᵢ ≥ 0
    W_diag = [w[svec_idx(C, i, i)] for i in 1:d]
    JuMP.@constraints(oa_model, begin
        v >= 0
        W_diag .>= 0
    end)

    # cuts using values of p = 1 and R = r₀ I
    u = cache.oa_s[1]
    r_vals = init_r_vals(cache.h)
    for r0 in r_vals
        R_diag = fill(r0, d)
        q = per_sepspec(h_conj, cache.h, 1.0, R_diag)
        JuMP.@constraint(oa_model, u + q * v + JuMP.dot(R_diag, W_diag) >= 0)
    end
    return 1 + d + length(r_vals)
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::MatrixEpiPerSepSpectralCache{Primal},
    oa_model::JuMP.Model,
)
    R = cache.W_temp
    @views svec_to_smat!(R, z[3:end], rt2)
    F = eigen(Hermitian(R, :U))
    V = F.vectors
    ω = F.values

    cuts = AE[]
    if max(abs(z[2])) < 1e-7
        # TODO decide when to add
        # add eigenvector cuts
        num_pos = count(>(1e-7), ω)
        @views V_neg = V[:, 1:num_pos] * Diagonal(sqrt.(ω[1:num_pos]))
        @views w = cache.oa_s[3:end]
        cuts = _get_psd_cuts(V_neg, w, cache, oa_model)
    end

    # add epigraph cut
    cut = _get_cut(z[1], ω, V, cache, oa_model)
    push!(cuts, cut)
    return cuts
end

function MOIPajarito.Cones.get_sep_cuts(
    cache::MatrixEpiPerSepSpectralCache{Primal},
    oa_model::JuMP.Model,
)
    # check s ∉ K
    Ws = cache.W_temp
    @views svec_to_smat!(Ws, cache.s[3:end], rt2)
    F = eigen(Hermitian(Ws, :U))
    V = F.vectors
    ω = F.values
    num_neg = count(<(-1e-7), ω)
    @assert issorted(ω)

    cuts = AE[]
    if !iszero(num_neg)
        # add eigenvector cuts
        V_neg = V[:, 1:num_neg]
        @views w = cache.oa_s[3:end]
        cuts = _get_psd_cuts(V_neg, w, cache, oa_model)
    end

    us = cache.s[1]
    v_pos = max(cache.s[2], 1e-7)
    ω_pos = max.(ω, 1e-7)
    us - per_sepspec(h_val, cache.h, v_pos, ω_pos) > -1e-7 && return AE[]

    # gradient cut is (1, h⋆(rω), V * Diagonal(rω) * V') at rω = -h'(ω / vs)
    rω = -h_grad(ω_pos / v_pos, cache.h)
    cut = _get_cut(1.0, rω, V, cache, oa_model)
    push!(cuts, cut)
    return cuts
end

function _get_cut(
    p::Float64,
    rω::Vector{Float64},
    V::Matrix{C},
    cache::MatrixEpiPerSepSpectralCache{Primal, C},
    oa_model::JuMP.Model,
) where {C}
    # strengthened cut is (p, p * h⋆(rω / p), V * Diagonal(rω) * V')
    @assert p > 1e-12
    if h_conj_dom_pos(cache.h)
        @assert all(>(1e-12), rω)
    end
    q = per_sepspec(h_conj, cache.h, p, rω)
    R = V * Diagonal(rω) * V'
    r = smat_to_svec!(cache.w_temp, R, rt2)
    z = vcat(p, q, r)
    cut = dot_expr(z, cache.oa_s, oa_model)
    return cut
end

# dual cone functions
