#=
real symmetric or complex Hermitian (svec scaled triangle) domain 𝕊ᵈ₊
=#

mutable struct MatrixEpiPerSepSpectral{D <: PrimDual, C <: RealCompF} <: Cone
    oa_s::Vector{AE}
    h::SepSpectralFun
    d::Int
    w_temp::Vector{RealF}
    W_temp::Matrix{C}
    function MatrixEpiPerSepSpectral{D, C}() where {D <: PrimDual, C <: RealCompF}
        return new{D, C}()
    end
end

function create_sepspectral_cache(
    ::Type{MatrixCSqr{RealF, C}},
    use_dual::Bool,
    d::Int,
    ::Bool,
) where {C <: RealCompF}
    D = primal_or_dual(use_dual)
    cache = MatrixEpiPerSepSpectral{D, C}()
    cache.w_temp = zeros(RealF, svec_length(C, d))
    cache.W_temp = zeros(C, d, d)
    return cache
end

function MOIPajarito.Cones.add_init_cuts(
    cache::MatrixEpiPerSepSpectral{D, C},
    oa_model::JuMP.Model,
) where {D, C}
    h = cache.h
    d = cache.d
    (u, v) = cache.oa_s[epi_per_idxs(D)]
    @views w = cache.oa_s[3:end]
    W_diag = [w[svec_idx(C, i, i)] for i in 1:d]

    # variable bounds
    JuMP.@constraint(oa_model, v >= 0)
    if dom_pos(D, h)
        JuMP.@constraint(oa_model, W_diag .>= 0)
    end

    # cuts using values of p = 1 and R = r₀ I
    r_vals = init_r_vals(D, h)
    for r0 in r_vals
        R_diag = fill(r0, d)
        q = per_sepspec(sepspec_dual(D), h, 1.0, R_diag)
        JuMP.@constraint(oa_model, u + q * v + JuMP.dot(R_diag, W_diag) >= 0)
    end
    return 1 + d + length(r_vals)
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::MatrixEpiPerSepSpectral{D},
    oa_model::JuMP.Model,
) where {D}
    (p, q) = z[epi_per_idxs(D)]
    R = cache.W_temp
    @views svec_to_smat!(R, z[3:end], rt2)
    F = eigen(Hermitian(R, :U))
    V = F.vectors
    ω = F.values

    cuts = AE[]
    if abs(q) < 1e-7
        # TODO only if domain for this part is pos

        # TODO decide when to add
        # add eigenvector cuts
        num_pos = count(>(1e-7), ω)
        @views V_neg = V[:, 1:num_pos] * Diagonal(sqrt.(ω[1:num_pos]))
        @views w = cache.oa_s[3:end]
        cuts = _get_psd_cuts(V_neg, w, cache, oa_model)
    end

    # add epigraph cut
    cut = _get_cut(p, ω, V, cache, oa_model)
    push!(cuts, cut)
    return cuts
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::MatrixEpiPerSepSpectral{Primal},
    oa_model::JuMP.Model,
)
    Ws = cache.W_temp
    @views svec_to_smat!(Ws, s[3:end], rt2)
    F = eigen(Hermitian(Ws, :U))
    V = F.vectors
    ω = F.values

    # TODO only if domain for this part is pos
    num_neg = count(<(-1e-7), ω)
    @assert issorted(ω)

    cuts = AE[]
    if !iszero(num_neg)
        # add eigenvector cuts
        V_neg = V[:, 1:num_neg]
        @views w = cache.oa_s[3:end]
        cuts = _get_psd_cuts(V_neg, w, cache, oa_model)
    end

    (us, vs) = s[epi_per_idxs(D)]
    v_pos = max(vs, 1e-7)

    # TODO only if domain for this part is pos
    ω_pos = max.(ω, 1e-7)
    us - per_sepspec(h_val, cache.h, v_pos, ω_pos) > -1e-7 && return AE[]

    # gradient cut is (1, h⋆(rω), V * Diagonal(rω) * V') at rω = -h'(ω / vs)
    rω = -h_grad(ω_pos / v_pos, cache.h)
    cut = _get_cut(1.0, rω, V, cache, oa_model)
    push!(cuts, cut)
    return cuts
end

function _get_cut(
    p::RealF,
    rω::Vector{RealF},
    V::Matrix{C},
    cache::MatrixEpiPerSepSpectral{Primal, C},
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
