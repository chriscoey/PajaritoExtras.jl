# all tests

import Test
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

import Hypatia
hypatia = MOI.OptimizerWithAttributes(
    Hypatia.Optimizer,
    MOI.Silent() => true,
    "near_factor" => 1000,
    "tol_feas" => 1e-10,
    "tol_rel_opt" => 1e-9,
    "tol_abs_opt" => 1e-8,
    "tol_illposed" => 1e-9,
    "tol_slow" => 2e-2,
    "tol_inconsistent" => 1e-7,
)

# println("starting PajaritoExtras tests")
# Test.@testset "PajaritoExtras tests" begin
#     include("JuMP_tests.jl")
#     TestJuMP.runtests(gurobi, hypatia)
# end

println("starting examples tests")
Test.@testset "examples tests" begin
    include("../examples/polyfacilitylocation/JuMP.jl")
end
