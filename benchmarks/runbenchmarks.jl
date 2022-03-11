#=
run from the benchmarks folder

~/julia/julia runbenchmarks.jl &> raw/bench.txt
=#

include(joinpath(@__DIR__, "..", "examples", "Examples.jl"))
using Test

# uncomment path for writing to results CSV
# csv = nothing
csv = joinpath(mkpath(joinpath(@__DIR__, "raw")), "bench.csv")

import MathOptInterface
const MOI = MathOptInterface

import Gurobi
gurobi = MOI.OptimizerWithAttributes(
    Gurobi.Optimizer,
    MOI.Silent() => true,
    "IntFeasTol" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "MIPGap" => 1e-9,
    "DualReductions" => 0, # fix infeasible or unbounded statuses
    "InfUnbdInfo" => 1, # get ray of the primal OA continuous relaxation
)

# default MOIPajarito options
options = (;
    iteration_limit = 10000,
    time_limit = 600.0,
    verbose = true,
    oa_solver = gurobi,
    use_iterative_method = true,
    use_init_fixed_oa = false,
)

# instance sets to run
inst_sets = [
    # generic nat vs ext:
    "nat",
    "nat_noext",
    "ext",
    # # experiment design:
    # "nat_rtdet",
    # "ext_rtdet",
    # "nat_entr",
    # "ext_entr",
    # # PWL formulations:
    # "sos2",
    # "logib",
    # "cc",
]

# list of names of JuMP examples to run
examples = [
    # # PSD:
    # "completablepsd",
    # # WSOS:
    # "polyfacilitylocation",
    # "polyregression",
    # "twostagestochastic",
    # norm:
    "matrixcompletion",
    # "matrixdecomposition",
    # "matrixregression",
    # # spectral function:
    # "experimentdesign",
    # "inversecovariance",
    # "vectorregression",
    # # nonconvex:
    # "ballpacking",
    # "modulardesign",
]

Examples.load_examples(examples)

function run_benchmarks()
    @testset "benchmarks" begin
        # Examples.run_examples(examples, inst_sets, options, csv, false) # TODO final run
        Examples.run_examples(examples, inst_sets, options, csv, true)
    end
    println()
    return
end

run_benchmarks()
