#=
run from the benchmarks folder
~/julia/julia run.jl &> raw/bench.txt
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

# default Pajarito options
options = (;
    iteration_limit = 100000,
    time_limit = 600.0,
    verbose = true,
    oa_solver = gurobi,
    use_iterative_method = true,
    use_init_fixed_oa = false,
)

# instance sets to run
inst_sets = [
    # test instances (comment out usually):
    # "test",
    # generic formulations:
    "nat",
    "noext",
    "ext",
    # knapsack:
    "cont_geo",
    "cont_noext_geo",
    "cont_log",
    "cont_noext_log",
    "cont_inv",
    "cont_noext_inv",
    "nat_geo",
    "noext_geo",
    "ext_geo",
    # experiment design:
    "nat_rtdet",
    "ext_rtdet",
    "nat_entr",
    "ext_entr",
    # PWL formulations:
    "sos2",
    "logib",
    "cc",
]

# list of names of JuMP examples to run
examples = [
    # PSD:
    "completablepsd",
    # WSOS:
    "polyfacilitylocation",
    "polyregression",
    "twostagestochastic",
    # spectral norm:
    "matrixcompletion",
    "matrixdecomposition",
    "matrixregression",
    # spectral function:
    "experimentdesign",
    "inversecovariance",
    "knapsack",
    "vectorregression",
    # MIP formulations:
    "ballpacking",
    "modulardesign",
]

skip_limit = false
# skip_limit = true

Examples.load_examples(examples)

function run_benchmarks()
    @testset verbose = true "benchmarks" begin
        Examples.run_examples(examples, inst_sets, options, csv, skip_limit)
    end
    println()
    return
end

run_benchmarks()
