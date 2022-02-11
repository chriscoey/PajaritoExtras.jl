# all tests

using Test
import MathOptInterface
const MOI = MathOptInterface

import Gurobi
oa_solver = MOI.OptimizerWithAttributes(
    Gurobi.Optimizer,
    MOI.Silent() => true,
    "IntFeasTol" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "MIPGap" => 1e-10,
    "DualReductions" => 0, # fixes infeasible or unbounded status
)

# import GLPK
# oa_solver = MOI.OptimizerWithAttributes(
#     GLPK.Optimizer,
#     MOI.Silent() => true,
#     "tol_int" => 1e-9,
#     "tol_bnd" => 1e-9,
#     "mip_gap" => 1e-9,
# )

println("starting PajaritoExtras tests")
@testset "PajaritoExtras tests" begin
    include("JuMP_tests.jl")
    TestJuMP.runtests(oa_solver)
end
