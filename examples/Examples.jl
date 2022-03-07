"""
PajaritoExtras examples and script utilities.
"""
module Examples

using LinearAlgebra
# delete later, affects qr. see https://github.com/JuliaLang/julia/pull/40623
if VERSION < v"1.7.0-DEV.1188"
    const ColumnNorm = Val{true}
end
using SparseArrays

import Random
# this is a workaround for randn's lack of support for BigFloat
Random.randn(R::Type{BigFloat}, dims::Vararg{Int, N} where {N}) = R.(randn(dims...))
function Random.randn(R::Type{Complex{BigFloat}}, dims::Vararg{Int, N} where {N})
    return R.(randn(ComplexF64, dims...))
end

using Test
using Printf
import DataFrames
import CSV
import DataStructures: OrderedDict

import Hypatia
import Hypatia.PolyUtils
import Hypatia.Cones
import Hypatia.Solvers

import MOIPajarito
import PajaritoExtras

import JuMP
const MOI = JuMP.MOI
const MOIU = MOI.Utilities
const JuMPScalar = JuMP.AbstractJuMPScalar

abstract type ExampleInstance end

include("benchmark_utils.jl")
include("JuMP_utils.jl")

const OPT_OR_LIMIT = [MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.ITERATION_LIMIT]

# TODO maybe move to Hypatia
include("WSOS_utils.jl")
include("specnuc_utils.jl")
include(joinpath(pkgdir(Hypatia), "examples", "spectral_functions_JuMP.jl"))
include("SOS2_utils.jl")

include(joinpath(@__DIR__, "..", "test", "model_utils.jl"))

# load all examples
function load_examples(examples::Vector{String})
    for ex in examples
        include(joinpath(@__DIR__, ex, "model.jl"))
    end
end

end
