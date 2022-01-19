"""
PajaritoExtras examples and script utilities.
"""
module Examples

using LinearAlgebra
# delete later, affects qr. see https://github.com/JuliaLang/julia/pull/40623
if VERSION < v"1.7.0-DEV.1188"
    const ColumnNorm = Val{true}
end

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
# import Hypatia.Models
import Hypatia.Solvers

import MOIPajarito
import PajaritoExtras

import JuMP
const MOI = JuMP.MOI
const MOIU = MOI.Utilities

abstract type ExampleInstance end

include("JuMP_utils.jl")
include("benchmark_utils.jl")
include("WSOS_JuMP.jl") # TODO maybe move to Hypatia
include(joinpath(pkgdir(Hypatia), "examples", "spectral_functions_JuMP.jl"))

# list of names of JuMP examples to run
const JuMP_examples = [
    # "experimentdesign",
    # "matrixdecomposition",
    "polyfacilitylocation",
    # TODO
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
