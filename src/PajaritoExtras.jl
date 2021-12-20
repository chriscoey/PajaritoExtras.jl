
module PajaritoExtras

const RealOrComplex = Union{Float64, ComplexF64}
const rt2 = sqrt(2.0)
const irt2 = inv(rt2)

using LinearAlgebra

import JuMP
const MOI = JuMP.MOI
const VR = JuMP.VariableRef
const AE = JuMP.AffExpr

# to create symmetric containers of AbstractJuMPScalars
LinearAlgebra.hermitian_type(::Type{T}) where {T <: JuMP.AbstractJuMPScalar} = T
LinearAlgebra.hermitian(scalar::JuMP.AbstractJuMPScalar, ::Symbol) = scalar
Base.isreal(::JuMP.AbstractJuMPScalar) = true
Base.real(scalar::JuMP.AbstractJuMPScalar) = scalar
# for linear algebra operations involving transposes
# LinearAlgebra.adjoint(scalar::Complex{<:JuMP.AbstractJuMPScalar}) = ?

import Hypatia
import Hypatia.Cones: vec_copyto!, svec_length, svec_side, smat_to_svec!, svec_to_smat!

import MOIPajarito
# import MOIPajarito.Cones: Extender, Unextended, Extended, extender
import MOIPajarito.Cones: ConeCache, clean_array!, dot_expr

include("possemideftri.jl")

# supported cones for outer approximation
const OACone = Union{Hypatia.PosSemidefTriCone}

# cone must be supported by both Pajarito and the conic solver
function MOI.supports_constraint(
    opt::MOIPajarito.Optimizer,
    F::Type{<:Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64}}},
    S::Type{<:OACone},
)
    return MOI.supports_constraint(MOIPajarito.get_conic_opt(opt), F, S)
end

end
