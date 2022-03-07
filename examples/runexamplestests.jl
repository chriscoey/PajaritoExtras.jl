#=
run examples tests from the examples folder
=#

include(joinpath(@__DIR__, "Examples.jl"))

# uncomment path for writing to results CSV
results_path = nothing
# results_path = joinpath(mkpath(joinpath(@__DIR__, "..", "benchmarks", "raw")), "tests.csv")

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
    iteration_limit = 500,
    time_limit = 120.0,
    # verbose = false,
    oa_solver = gurobi,
    # oa_solver = glpk,
    use_iterative_method = true,
    # use_extended_form = false,
    # solve_relaxation = false,
    # solve_subproblems = false,
    # use_init_fixed_oa = false,
)

perf = Examples.run_examples(["test"], default_options, results_path)
println()
