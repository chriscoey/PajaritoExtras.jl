#=
run examples tests from the examples folder
=#

# uncomment path for writing to results CSV
# results_path = nothing
results_path = joinpath(mkpath(joinpath(@__DIR__, "..", "benchmarks", "raw")), "bench.csv")

import MathOptInterface
const MOI = MathOptInterface

# import GLPK
# glpk = MOI.OptimizerWithAttributes(
#     GLPK.Optimizer,
#     MOI.Silent() => true,
#     "tol_int" => 1e-10,
#     "tol_bnd" => 1e-10,
#     "mip_gap" => 1e-10,
# )

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

default_options = (;
    verbose = true,
    # verbose = false,
    oa_solver = gurobi,
    # oa_solver = glpk,
    use_extended_form = true,
    use_iterative_method = true,
    # use_iterative_method = false,
    solve_relaxation = true,
    # solve_subproblems = true,
    # solve_relaxation = false,
    solve_subproblems = false,
    iteration_limit = 500,
    time_limit = 120.0,
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
