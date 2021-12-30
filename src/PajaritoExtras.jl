
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
import MOIPajarito.Cones: NatExt, Nat, Ext, extender
import MOIPajarito.Cones: Cone, clean_array!, dot_expr

const RealComp = Union{Float64, ComplexF64}

abstract type PrimDual end
struct Primal <: PrimDual end
struct Dual <: PrimDual end
primal_or_dual(use_dual::Bool) = (use_dual ? Dual : Primal)

const rt2 = sqrt(2.0)
const irt2 = inv(rt2)

include("possemideftri.jl")
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
    Hypatia.PosSemidefTriCone{Float64, <:RealComp},
    Hypatia.EpiNormEuclCone{Float64},
    Hypatia.EpiPerSquareCone{Float64},
    Hypatia.EpiNormInfCone{Float64, <:RealComp},
    Hypatia.EpiNormSpectralCone{Float64, <:RealComp},
    Hypatia.HypoGeoMeanCone{Float64},
    Hypatia.HypoRootdetTriCone{Float64, <:RealComp},
    Hypatia.EpiPerSepSpectralCone{Float64},
    Hypatia.WSOSInterpNonnegativeCone{Float64, <:RealComp},
}

# cone must be supported by both Pajarito and the conic solver
function MOI.supports_constraint(
    opt::MOIPajarito.Optimizer,
    F::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64}}},
    S::Type{<:OACone},
)
    return MOI.supports_constraint(MOIPajarito.get_conic_opt(opt), F, S)
end

# eigenvector cuts for a PSD constraint W ⪰ 0
function _get_psd_cuts(
    R_eig::Matrix{C},
    oa_w::AbstractVector{AE},
    cache::Cone,
    oa_model::JuMP.Model,
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
        cut = dot_expr(R_vec_i, oa_w, oa_model)
        push!(cuts, cut)
    end
    return cuts
end

# TODO move to Hypatia array utilities?
svec_idx(::Type{Float64}, row::Int, col::Int) = Hypatia.Cones.svec_idx(row, col)

function svec_idx(::Type{ComplexF64}, row::Int, col::Int)
    if row < col
        (row, col) = (col, row)
    end
    return (row - 1) * row + col
end

function geomean(w::AbstractVector{Float64}; min_val::Float64 = 1e-12)
    w = max.(w, min_val)
    any(<=(eps()), w) && return 0.0
    return exp(sum(log, w) / length(w))
end

end
