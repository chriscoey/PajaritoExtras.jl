#=
run examples tests from the examples folder
=#

# uncomment path for writing to results CSV
results_path = joinpath(mkpath(joinpath(@__DIR__, "..", "benchmarks", "raw")), "bench.csv")
# results_path = nothing

import MathOptInterface
const MOI = MathOptInterface
import Gurobi
gurobi = MOI.OptimizerWithAttributes(
    Gurobi.Optimizer,
    MOI.Silent() => true,
    "IntFeasTol" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "MIPGap" => 1e-10,
    "DualReductions" => 0, # fixes infeasible or unbounded status
)

# default MOIPajarito options
default_options = (;
    verbose = true,
    # verbose = false,
    oa_solver = gurobi,
    use_extended_form = true,
    use_iterative_method = true,
    # use_iterative_method = false,
    # iteration_limit = 30,
    # time_limit = 120.0,
    # debug_cuts = use_iterative_method,
)

# instance sets to run
inst_sets = [
    "test",
    # benchmarks:
    # "nat",
    # "ext",
]

include(joinpath(@__DIR__, "Examples.jl"))
perf = Examples.run_examples(inst_sets, default_options, results_path)
println()
