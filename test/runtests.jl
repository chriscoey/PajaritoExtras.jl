# all tests

using Test
import MathOptInterface
const MOI = MathOptInterface

import GLPK
glpk = MOI.OptimizerWithAttributes(
    GLPK.Optimizer,
    MOI.Silent() => true,
    "tol_int" => 1e-10,
    "tol_bnd" => 1e-10,
    "mip_gap" => 1e-10,
)

println("starting PajaritoExtras tests")
@testset "PajaritoExtras tests" begin
    include("JuMP_tests.jl")
    TestJuMP.runtests(glpk)
end
