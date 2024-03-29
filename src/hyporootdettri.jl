#=
hypograph of root-determinant of real symmetric or complex Hermitian
positive semidefinite matrices (svec scaled triangle)
(u, w) : W ⪰ 0, u ≤ rtdet(W)
where rtdet(W) = (∏ᵢ λᵢ(W))^(1/d)

dual cone is
(u, w) : u ≤ 0, W ⪰ 0, u ≥ -d * rtdet(W)
=#

mutable struct HypoRootdetTri{C <: RealCompF} <: Cache
    oa_s::Vector{AE}
    d::Int
    w_temp::Vector{RealF}
    W_temp::Matrix{C}
    HypoRootdetTri{C}() where {C <: RealCompF} = new{C}()
end

function Pajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.HypoRootdetTriCone{RealF, C},
    ::Optimizer,
) where {C <: RealCompF}
    cache = HypoRootdetTri{C}()
    cache.oa_s = oa_s
    dim = MOI.dimension(cone)
    d = cache.d = svec_side(C, dim - 1)
    cache.w_temp = zeros(RealF, dim - 1)
    cache.W_temp = zeros(C, d, d)
    return cache
end

function Pajarito.Cones.add_init_cuts(cache::HypoRootdetTri{C}, opt::Optimizer) where {C}
    # add variable bounds W_ii ≥ 0
    d = cache.d
    @views w = cache.oa_s[2:end]
    W_diag = [w[svec_idx(C, i, i)] for i in 1:d]
    JuMP.@constraint(opt.oa_model, W_diag .>= 0)
    opt.use_init_fixed_oa || return

    # add cut (-d, I) on (u, W), a linearization at W = I
    u = cache.oa_s[1]
    JuMP.@constraint(opt.oa_model, -d * u + sum(W_diag) >= 0)
    return
end

function Pajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::HypoRootdetTri,
    opt::Optimizer,
)
    p = z[1]
    R = cache.W_temp
    @views svec_to_smat!(R, z[2:end], rt2)
    F = eigen(Hermitian(R, :U))
    V = F.vectors
    ω = F.values

    cuts = AE[]
    if p > -1e-7
        # add eigenvector cuts
        num_pos = count(>(1e-7), ω)
        @views V_neg = V[:, 1:num_pos] * Diagonal(sqrt.(ω[1:num_pos]))
        @views w = cache.oa_s[2:end]
        cuts = _get_psd_cuts(V_neg, w, cache, opt)
    end

    # add rootdet cut
    ω_pos = max.(ω, 1e-7)
    cut = _get_cut(ω_pos, V, cache, opt)
    push!(cuts, cut)
    return cuts
end

function Pajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::HypoRootdetTri,
    opt::Optimizer,
)
    us = s[1]
    Ws = cache.W_temp
    @views svec_to_smat!(Ws, s[2:end], rt2)
    F = eigen(Hermitian(Ws, :U))
    V = F.vectors
    ω = F.values
    num_neg = count(<(-opt.tol_feas), ω)
    if iszero(num_neg) && us < opt.tol_feas
        return AE[]
    end
    @assert issorted(ω)

    cuts = AE[]
    if !iszero(num_neg)
        # add eigenvector cuts
        V_neg = V[:, 1:num_neg]
        @views w = cache.oa_s[2:end]
        cuts = _get_psd_cuts(V_neg, w, cache, opt)
    end

    geomean(ω) - us > -opt.tol_feas && return AE[]
    # gradient cut is (-1, V * Diagonal(geom(ω) / d ./ ω) * V')
    ω_pos = max.(ω, 1e-7)
    rω = (geomean(ω_pos) / cache.d) ./ ω_pos
    cut = _get_cut(rω, V, cache, opt)
    push!(cuts, cut)
    return cuts
end

function _get_cut(
    rω::Vector{RealF},
    V::Matrix{C},
    cache::HypoRootdetTri{C},
    opt::Optimizer,
) where {C}
    # strengthened cut is (-d * geom(rω), V * Diagonal(rω) * V')
    p = -cache.d * geomean(rω)
    R = V * Diagonal(rω) * V'
    r = smat_to_svec!(cache.w_temp, R, rt2)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    return JuMP.@expression(opt.oa_model, p * u + dot(r, w))
end
