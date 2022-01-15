#=
real symmetric or complex Hermitian sparse positive semidefinite cone (svec format)
sparsity pattern given by lower triangle row and column indices, with all diagonal elements
dual cone is the PSD-completable matrices

TODO use Arpack
=#

mutable struct PosSemidefTriSparse{D <: PrimDual, C <: RealCompF} <: Cache
    oa_s::Vector{AE}
    side::Int
    row_idxs::Vector{Int}
    col_idxs::Vector{Int}
    PosSemidefTriSparse{D, C}() where {D <: PrimDual, C <: RealCompF} = new{D, C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.PosSemidefTriSparseCone{<:Hypatia.Cones.PSDSparseImpl, RealF, C},
    ::Bool,
) where {C <: RealCompF}
    D = primal_or_dual(cone.use_dual)
    cache = PosSemidefTriSparse{D, C}()
    cache.oa_s = oa_s
    cache.side = cone.side
    @assert length(oa_s) == MOI.dimension(cone)
    cache.row_idxs = cone.row_idxs
    cache.col_idxs = cone.col_idxs
    return cache
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
    # TODO eig cuts?
    cut = dot_expr(z, cache.oa_s, opt)
    return [cut]
end

# function MOIPajarito.Cones.get_sep_cuts(
#     ::Vector{RealF},
#     ::PosSemidefTriSparse{Prim},
#     ::Optimizer,
# )
#     # TODO eig cuts
#     @warn("no separation oracle implemented for PosSemidefTriSparse", maxlog = 1)
#     return AE[]
# end

# dual cone functions

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::PosSemidefTriSparse{Dual},
    opt::Optimizer,
)
    cut = dot_expr(z, cache.oa_s, opt)
    return [cut]
end

# function MOIPajarito.Cones.get_sep_cuts(
#     s::Vector{RealF},
#     cache::PosSemidefTriSparse{Dual},
#     opt::Optimizer,
# )
#     # TODO think about completion algorithms
#     @warn("no separation oracle implemented for dual PosSemidefTriSparse", maxlog = 1)
#     return AE[]
# end
