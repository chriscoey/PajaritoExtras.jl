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

include("JuMP_utils.jl")
include("benchmark_utils.jl")

# TODO maybe move to Hypatia
include("WSOS_utils.jl")
include("specnuc_utils.jl")

include(joinpath(pkgdir(Hypatia), "examples", "spectral_functions_JuMP.jl"))

include(joinpath(@__DIR__, "..", "test", "model_utils.jl"))

# list of names of JuMP examples to run
const JuMP_examples = [
    "completablepsd",
    "experimentdesign",
    "inversecovariance",
    "matrixcompletion",
    "matrixdecomposition",
    "matrixregression",
    "polyfacilitylocation",
    "polyregression",
    "twostagestochastic",
    "vectorregression",
]

# load all examples
for ex in JuMP_examples
    include(joinpath(@__DIR__, ex, "model.jl"))
end

# build ordered dictionary of all test instances
function get_test_instances()
    test_insts = OrderedDict()
    for ex in JuMP_examples
        test_insts[ex] = include(joinpath(@__DIR__, ex, "instances.jl"))
    end
    return test_insts
end

end
