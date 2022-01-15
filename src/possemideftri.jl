#=
real symmetric or complex Hermitian positive semidefinite matrices (svec scaled triangle)
W ⪰ 0
self-dual
=#

mutable struct PosSemidefTri{C <: RealCompF} <: Cache
    oa_s::Vector{AE}
    d::Int
    w_temp::Vector{RealF}
    W_temp::Matrix{C}
    PosSemidefTri{C}() where {C <: RealCompF} = new{C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.PosSemidefTriCone{RealF, C},
    ::Bool,
) where {C <: RealCompF}
    cache = PosSemidefTri{C}()
    cache.oa_s = oa_s
    dim = MOI.dimension(cone)
    d = cache.d = svec_side(C, dim)
    cache.w_temp = zeros(RealF, dim)
    cache.W_temp = zeros(C, d, d)
    return cache
end

# initial OA polyhedron is the dual cone of diagonally dominant matrices
function MOIPajarito.Cones.add_init_cuts(cache::PosSemidefTri{C}, opt::Optimizer) where {C}
    d = cache.d
    w = cache.oa_s
    # W_ii ≥ 0
    JuMP.@constraint(opt.oa_model, [i in 1:d], w[svec_idx(C, i, i)] >= 0)
    # a linearization of W_ii * W_jj ≥ |W_ij|^2
    for j in 1:d, i in 1:(j - 1)
        _add_init_cuts(i, j, cache, opt)
    end
    return (C == CompF ? 2d^2 - d : d^2)
end

# real: cuts on (w_ii, w_jj, w_ij) are (1, 1, ±rt2), ∀i != j
function _add_init_cuts(i::Int, j::Int, cache::PosSemidefTri{RealF}, opt::Optimizer)
    w = cache.oa_s
    w_ii = w[svec_idx(RealF, i, i)]
    w_jj = w[svec_idx(RealF, j, j)]
    w_ij = w[svec_idx(RealF, i, j)]
    JuMP.@constraints(opt.oa_model, begin
        w_ii + w_jj + rt2 * w_ij >= 0
        w_ii + w_jj - rt2 * w_ij >= 0
    end)
    return
end

# complex: cuts on smat (w_ii, w_jj, wr_ij, wi_ij) are (1, 1, ±1, ±1), ∀i != j
function _add_init_cuts(i::Int, j::Int, cache::PosSemidefTri{CompF}, opt::Optimizer)
    w = cache.oa_s
    w_ii = w[svec_idx(CompF, i, i)]
    w_jj = w[svec_idx(CompF, j, j)]
    ij_idx = svec_idx(CompF, i, j)
    wr_ij = w[ij_idx]
    wi_ij = w[ij_idx + 1]
    JuMP.@constraints(opt.oa_model, begin
        w_ii + w_jj + wr_ij + wi_ij >= 0
        w_ii + w_jj + wr_ij - wi_ij >= 0
        w_ii + w_jj - wr_ij + wi_ij >= 0
        w_ii + w_jj - wr_ij - wi_ij >= 0
    end)
    return
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::PosSemidefTri,
    opt::Optimizer,
)
    # strengthened cuts from eigendecomposition are λᵢ * rᵢ * rᵢ'
    R = cache.W_temp
    svec_to_smat!(R, z, rt2)
    F = eigen(Hermitian(R, :U), 1e-8, Inf) # TODO tune
    isempty(F.values) && return AE[]
    R_eig = F.vectors * Diagonal(sqrt.(F.values))
    return _get_psd_cuts(R_eig, cache.oa_s, cache, opt)
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::PosSemidefTri,
    opt::Optimizer,
)
    Ws = cache.W_temp
    svec_to_smat!(Ws, s, rt2)
    F = eigen(Hermitian(Ws, :U), -Inf, -1e-7)
    isempty(F.values) && return AE[]
    return _get_psd_cuts(F.vectors, cache.oa_s, cache, opt)
end
