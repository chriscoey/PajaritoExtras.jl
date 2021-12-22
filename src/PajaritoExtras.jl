
module PajaritoExtras

const RealOrComplex = Union{Float64, ComplexF64}
const rt2 = sqrt(2.0)
const irt2 = inv(rt2)

using LinearAlgebra

import JuMP
const MOI = JuMP.MOI
const VR = JuMP.VariableRef
const AE = JuMP.AffExpr

import Hypatia
import Hypatia.Cones: vec_copyto!, svec_length, svec_side, smat_to_svec!, svec_to_smat!

import MOIPajarito
import MOIPajarito.Cones: Extender, Unextended, Extended, extender
import MOIPajarito.Cones: ConeCache, clean_array!, dot_expr

include("possemideftri.jl")
include("epinormeucl.jl")
include("epipersquare.jl")
include("hypogeomean.jl")

# supported cones for outer approximation
const OACone = Union{
    Hypatia.PosSemidefTriCone{Float64, <:RealOrComplex},
    Hypatia.EpiNormEuclCone{Float64},
    Hypatia.EpiPerSquareCone{Float64},
    Hypatia.HypoGeoMeanCone{Float64},
}

# cone must be supported by both Pajarito and the conic solver
function MOI.supports_constraint(
    opt::MOIPajarito.Optimizer,
    F::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64}}},
    S::Type{<:OACone},
)
    return MOI.supports_constraint(MOIPajarito.get_conic_opt(opt), F, S)
end

# TODO move to Hypatia array utilities?
svec_idx(::Type{Float64}, row::Int, col::Int) = Hypatia.Cones.svec_idx(row, col)

function svec_idx(::Type{ComplexF64}, row::Int, col::Int)
    if row < col
        (row, col) = (col, row)
    end
    return (row - 1) * row + col
end

function geomean(w::AbstractVector{Float64})
    any(<=(eps()), w) && return eps()
    return exp(sum(log, w) / length(w))
end

end
