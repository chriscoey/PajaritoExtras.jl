# all tests

using Test
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

# import GLPK
# glpk = MOI.OptimizerWithAttributes(
#     GLPK.Optimizer,
#     MOI.Silent() => true,
#     "tol_int" => 1e-10,
#     "tol_bnd" => 1e-10,
#     "mip_gap" => 1e-10,
# )

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

println("starting PajaritoExtras tests")
@testset "PajaritoExtras tests" begin
    include("JuMP_tests.jl")
    TestJuMP.runtests(gurobi, hypatia)
    # TestJuMP.runtests(glpk, hypatia)
end

# TODO
import LinearAlgebra
import JuMP
import Hypatia
import MOIPajarito
import PajaritoExtras
# TODO include spectral_functions_JuMP.jl file from Hypatia examples

opt = JuMP.optimizer_with_attributes(
    MOIPajarito.Optimizer,
    "verbose" => true,
    # "verbose" => false,
    "oa_solver" => gurobi,
    "conic_solver" => hypatia,
    "use_extended_form" => true,
    "use_iterative_method" => true,
    # "debug_cuts" => use_iterative_method,
    # "iteration_limit" => 30,
    # "time_limit" => 120.0,
)

println("starting examples tests")
@testset "examples tests" begin
    # include("../examples/polyfacilitylocation/JuMP.jl")
    # include("../examples/experimentdesign/JuMP.jl")
    # include("../examples/matrixdecomposition/JuMP.jl")
end
