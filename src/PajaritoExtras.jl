
module PajaritoExtras

using LinearAlgebra

import JuMP
const MOI = JuMP.MOI
const VR = JuMP.VariableRef
const AE = JuMP.AffExpr

import Hypatia
import Hypatia.Cones: vec_length, vec_copyto!, svec_length, svec_side
import Hypatia.Cones: smat_to_svec!, svec_to_smat!

import MOIPajarito
import MOIPajarito: Cache, Optimizer
import MOIPajarito.Cones: NatExt, Nat, Ext, nat_or_ext, clean_array!, dot_expr

const RealF = Float64
const CompF = ComplexF64
const RealCompF = Union{RealF, CompF}

abstract type PrimDual end
struct Prim <: PrimDual end
struct Dual <: PrimDual end
primal_or_dual(use_dual::Bool) = (use_dual ? Dual : Prim)

const rt2 = sqrt(2.0)
const irt2 = inv(rt2)

include("possemideftri.jl")
include("possemideftrisparse.jl")
include("epinormeucl.jl")
include("epipersquare.jl")
include("epinorminf.jl")
include("epinormspectral.jl")
include("hypogeomean.jl")
include("hyporootdettri.jl")
include("epipersepspectral/epipersepspectral.jl")
include("wsosinterpnonnegative.jl")

# supported cones for outer approximation
const OACone = Union{
    Hypatia.PosSemidefTriCone{RealF, <:RealCompF},
    # TODO maybe shouldn't have the PSDSparseImpl param in Hypatia (sep spectral doesn't do it):
    Hypatia.PosSemidefTriSparseCone{<:Hypatia.Cones.PSDSparseImpl, RealF, <:RealCompF},
    Hypatia.EpiNormEuclCone{RealF},
    Hypatia.EpiPerSquareCone{RealF},
    Hypatia.EpiNormInfCone{RealF, <:RealCompF},
    Hypatia.EpiNormSpectralCone{RealF, <:RealCompF},
    Hypatia.HypoGeoMeanCone{RealF},
    Hypatia.HypoRootdetTriCone{RealF, <:RealCompF},
    Hypatia.EpiPerSepSpectralCone{RealF},
    Hypatia.WSOSInterpNonnegativeCone{RealF, <:RealCompF},
}

# cone must be supported by both Pajarito and the conic solver
function MOI.supports_constraint(
    opt::MOIPajarito.Optimizer,
    F::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{RealF}}},
    S::Type{<:OACone},
)
    return MOI.supports_constraint(MOIPajarito.get_conic_opt(opt), F, S)
end

# fallback for separation cuts solves a separation subproblem
function MOIPajarito.Cones.get_sep_cuts(s::Vector{RealF}, cache::Cache, opt::Optimizer)
    # @warn("solving separation subproblem for $(typeof(cache))")
    # TODO only update model, only make unique models
    # TODO Hypatia options etc?
    # sep_model = JuMP.Model(Hypatia.Optimizer)
    # JuMP.set_optimizer_attribute(sep_model, MOI.Silent(), true)
    # ci = JuMP.@constraint(sep_model, s in cache.cone)
    # JuMP.optimize!(sep_model)

    # stat = JuMP.termination_status(sep_model)
    # @show stat
    # if stat in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
    #     return AE[]
    # elseif stat in (MOI.INFEASIBLE, MOI.ALMOST_INFEASIBLE)
    #     z = JuMP.getdual(ci)
    #     return MOIPajarito.Cones.get_subp_cuts(z, cache, opt)
    # end
    # @warn("separation subproblem status was $stat")
    return AE[]
end

# eigenvector cuts for a PSD constraint W ⪰ 0
function _get_psd_cuts(
    R_eig::Matrix{C},
    oa_w::AbstractVector{AE},
    cache::Cache,
    opt::Optimizer,
) where {C}
    cuts = AE[]
    R_vec_i = cache.w_temp
    R_i = cache.W_temp
    for i in 1:size(R_eig, 2)
        @views r_i = R_eig[:, i]
        # cuts from eigendecomposition are rᵢ * rᵢ'
        mul!(R_i, r_i, r_i')
        clean_array!(R_i) && continue
        smat_to_svec!(R_vec_i, R_i, rt2)
        cut = dot_expr(R_vec_i, oa_w, opt)
        push!(cuts, cut)
    end
    return cuts
end

# TODO move to Hypatia array utilities?
svec_idx(::Type{RealF}, row::Int, col::Int) = Hypatia.Cones.svec_idx(row, col)

function svec_idx(::Type{CompF}, row::Int, col::Int)
    if row < col
        (row, col) = (col, row)
    end
    return (row - 1) * row + col
end

function geomean(w::AbstractVector{RealF}; min_val::RealF = 1e-12)
    w = max.(w, min_val)
    any(<=(eps()), w) && return 0.0
    return exp(sum(log, w) / length(w))
end

end
