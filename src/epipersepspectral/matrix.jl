#=
real symmetric or complex Hermitian (svec scaled triangle) domain 𝕊ᵈ₊
=#

mutable struct MatrixEpiPerSepSpectral{D <: PrimDual, C <: RealCompF} <: Cache
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
    ::Optimizer,
) where {C <: RealCompF}
    D = primal_or_dual(use_dual)
    cache = MatrixEpiPerSepSpectral{D, C}()
    cache.w_temp = zeros(RealF, svec_length(C, d))
    cache.W_temp = zeros(C, d, d)
    return cache
end

function MOIPajarito.Cones.add_init_cuts(
    cache::MatrixEpiPerSepSpectral{D, C},
    opt::Optimizer,
) where {D, C}
    h = cache.h
    d = cache.d
    (u, v) = swap_epiper(D, cache.oa_s[1:2]...)
    @views w = cache.oa_s[3:end]
    W_diag = [w[svec_idx(C, i, i)] for i in 1:d]

    # variable bounds
    JuMP.@constraint(opt.oa_model, v >= 0)
    if dom_pos(D, h)
        JuMP.@constraint(opt.oa_model, W_diag .>= 0)
    end

    # cuts using values of p = 1 and R = r₀ I
    r_vals = init_r_vals(D, h)
    for r0 in r_vals
        R_diag = fill(r0, d)
        q = per_sepspec(conj_or_val(D), h, 1.0, R_diag)
        JuMP.@constraint(opt.oa_model, u + q * v + JuMP.dot(R_diag, W_diag) >= 0)
    end
    return 1 + d + length(r_vals)
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::MatrixEpiPerSepSpectral{D},
    opt::Optimizer,
) where {D}
    (p, q) = swap_epiper(D, z[1:2]...)
    R = cache.W_temp
    @views svec_to_smat!(R, z[3:end], rt2)
    F = eigen(Hermitian(R, :U))
    V = F.vectors
    ω = F.values
    cuts = AE[]

    if dom_pos(D, cache.h) && abs(q) < 1e-7
        # TODO decide when to add
        # add eigenvector cuts
        num_pos = count(>(1e-8), ω)
        @views V_neg = V[:, 1:num_pos] * Diagonal(sqrt.(ω[1:num_pos]))
        @views w = cache.oa_s[3:end]
        cuts = _get_psd_cuts(V_neg, w, cache, opt)
    end

    # add epigraph cut
    cut = _get_cut(p, ω, V, cache, opt)
    push!(cuts, cut)
    return cuts
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::MatrixEpiPerSepSpectral{D},
    opt::Optimizer,
) where {D}
    h = cache.h
    Ws = cache.W_temp
    @views svec_to_smat!(Ws, s[3:end], rt2)
    F = eigen(Hermitian(Ws, :U))
    V = F.vectors
    ω = F.values
    cuts = AE[]

    if dom_pos(D, h)
        num_neg = count(<(-opt.tol_feas), ω)
        if !iszero(num_neg)
            # add eigenvector cuts
            V_neg = V[:, 1:num_neg]
            @views w = cache.oa_s[3:end]
            cuts = _get_psd_cuts(V_neg, w, cache, opt)
        end

        ω = max.(ω, 1e-7)
    end

    (us, vs) = swap_epiper(D, s[1:2]...)
    v_pos = max(vs, 1e-7)
    if us - per_sepspec(val_or_conj(D), h, v_pos, ω) > -opt.tol_feas
        return cuts
    end

    # gradient cut is (1, h⋆(rω), V * Diagonal(rω) * V') at rω = -h'(ω / vs)
    rω = -h_grad(D, h, ω / v_pos)
    cut = _get_cut(1.0, rω, V, cache, opt)
    push!(cuts, cut)
    return cuts
end

function _get_cut(
    p::RealF,
    rω::Vector{RealF},
    V::Matrix{C},
    cache::MatrixEpiPerSepSpectral{D, C},
    opt::Optimizer,
) where {D, C}
    # strengthened cut is (p, p * h⋆(rω / p), V * Diagonal(rω) * V')
    q = per_sepspec(conj_or_val(D), cache.h, p, rω)
    (p, q) = swap_epiper(D, p, q)

    R = V * Diagonal(rω) * V'
    r = smat_to_svec!(cache.w_temp, R, rt2)

    z = vcat(p, q, r)
    cut = dot_expr(z, cache.oa_s, opt)
    return cut
end
