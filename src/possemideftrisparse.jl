#=
real symmetric or complex Hermitian sparse positive semidefinite cone (svec format)
sparsity pattern given by lower triangle row and column indices, with all diagonal elements
dual cone is the PSD-completable matrices
=#

mutable struct PosSemidefTriSparse{D <: PrimDual, C <: RealCompF} <: Cache
    oa_s::Vector{AE}
    side::Int
    row_idxs::Vector{Int}
    col_idxs::Vector{Int}
    sep_constr::CR
    PosSemidefTriSparse{D, C}() where {D <: PrimDual, C <: RealCompF} = new{D, C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.PosSemidefTriSparseCone{<:Hypatia.Cones.PSDSparseImpl, RealF, C},
    opt::Optimizer,
) where {C <: RealCompF}
    D = primal_or_dual(cone.use_dual)
    cache = PosSemidefTriSparse{D, C}()
    cache.oa_s = oa_s
    cache.side = cone.side
    @assert length(oa_s) == MOI.dimension(cone)
    cache.row_idxs = cone.row_idxs
    cache.col_idxs = cone.col_idxs
    if D == Dual
        cache.sep_constr = get_sep_constr(cone, opt)
    end
    return cache
end

# TODO maybe define Hypatia MOI cone equal/hash in Hypatia
function hash_cone(cone::Hypatia.PosSemidefTriSparseCone{I, RealF, C}) where {I, C}
    @assert cone.use_dual
    return hash(I) + hash(C) + hash(cone.side) + hash(cone.row_idxs) + hash(cone.col_idxs)
end

function MOIPajarito.Cones.add_init_cuts(
    cache::PosSemidefTriSparse{<:PrimDual, C},
    opt::Optimizer,
) where {C}
    # diagonal nonnegative
    # TODO for each offdiag, add linearization of 2x2 principal matrix condition
    num_cuts = 0
    incr = (C == CompF ? 2 : 1)
    i = 1
    for (r, c) in zip(cache.row_idxs, cache.col_idxs)
        if r == c
            JuMP.@constraint(opt.oa_model, cache.oa_s[i] >= 0)
            num_cuts += 1
            i += 1
        else
            i += incr
        end
    end
    @assert num_cuts == cache.side # all diagonal elements should be present
    return num_cuts
end

# primal cone functions

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::PosSemidefTriSparse{Prim},
    opt::Optimizer,
)
    cut = dot_expr(z, cache.oa_s, opt)
    return [cut]
end

function MOIPajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::PosSemidefTriSparse{Prim, C},
    opt::Optimizer,
) where {C}
    # TODO use Arpack for non-dense implementation
    # @warn("no separation oracle implemented for PosSemidefTriSparse", maxlog = 1)
    Λ = svec_to_smat_sparse(s, cache)
    F = eigen!(Hermitian(Λ, :L), -Inf, -1e-7)
    isempty(F.values) && return AE[]

    cuts = AE[]
    R_i = similar(Λ)
    for i in 1:length(F.values)
        @views r_i = F.vectors[:, i]
        # cuts from eigendecomposition are svec(rᵢ * rᵢ')
        mul!(R_i, r_i, r_i')
        clean_array!(R_i) && continue
        R_vec_i = smat_to_svec_sparse(R_i, cache)
        cut = dot_expr(R_vec_i, oa_w, opt)
        push!(cuts, cut)
    end
    return cuts
end

# dual cone functions

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::PosSemidefTriSparse{Dual},
    opt::Optimizer,
)
    cut = dot_expr(z, cache.oa_s, opt)
    return [cut]
end

# helpers
# TODO maybe refac with similar Hypatia functions (which should not have the cone as an argument)

function svec_to_smat_sparse(
    vec::AbstractVector{RealF},
    cache::PosSemidefTriSparse{<:PrimDual, RealF},
)
    @assert length(vec) == length(cache.row_idxs)
    mat = zeros(RealF, cache.side, cache.side)
    @inbounds for (idx, (i, j)) in enumerate(zip(cache.row_idxs, cache.col_idxs))
        x = vec[idx]
        if i != j
            x /= rt2
        end
        mat[i, j] = x
    end
    return mat
end

function svec_to_smat_sparse(
    vec::AbstractVector{RealF},
    cache::PosSemidefTriSparse{<:PrimDual, CompF},
)
    mat = zeros(CompF, cache.side, cache.side)
    idx = 1
    @inbounds for (i, j) in zip(cache.row_idxs, cache.col_idxs)
        if i == j
            mat[i, j] = vec[idx]
            idx += 1
        else
            mat[i, j] = Complex(vec[idx], vec[idx + 1]) / rt2
            idx += 2
        end
    end
    @assert idx == length(vec) + 1
    return mat
end

function smat_to_svec_sparse(
    mat::AbstractMatrix{RealF},
    cache::PosSemidefTriSparse{<:PrimDual, RealF},
)
    vec = zeros(length(cone.row_idxs))
    @inbounds for (idx, (i, j)) in enumerate(zip(cache.row_idxs, cache.col_idxs))
        x = mat[i, j]
        if i != j
            x *= rt2
        end
        vec[idx] = x
    end
    return vec
end

function smat_to_svec_sparse(
    mat::AbstractMatrix{CompF},
    cache::PosSemidefTriSparse{<:PrimDual, CompF},
)
    vec = zeros(2 * length(cache.row_idxs) - cache.side)
    idx = 1
    @inbounds for (i, j) in zip(cache.row_idxs, cache.col_idxs)
        x = mat[i, j]
        if i == j
            vec[idx] = real(x)
            idx += 1
        else
            x *= rt2
            vec[idx] = real(x)
            vec[idx + 1] = imag(x)
            idx += 2
        end
    end
    @assert idx == length(vec) + 1
    return vec
end
