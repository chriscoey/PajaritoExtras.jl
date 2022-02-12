# JuMP instance tests

module TestJuMP

using Test
using LinearAlgebra
import MathOptInterface
const MOI = MathOptInterface
import JuMP
import MOIPajarito
import PajaritoExtras
import PajaritoExtras: Prim, Dual, svec_idx
import Hypatia
import Hypatia.Cones: vec_length, svec_length, vec_copyto!
import Hypatia.PolyUtils

include("model_utils.jl")

# all instances
inst_all = String[
    "possemideftri1",
    "possemideftri2",
    "possemideftrisparse1",
    "possemideftrisparse2",
    "epinormeucl1",
    "epinormeucl2",
    "epipersquare1",
    "epipersquare2",
    "epinorminf1",
    "epinorminf2",
    "epinormspectraltri1",
    "epinormspectraltri2",
    "epinormspectral1",
    "epinormspectral2",
    "hypogeomean1",
    "hypogeomean2",
    "hyporootdettri1",
    "hyporootdettri2",
    "vector_epipersepspectral1",
    "vector_epipersepspectral2",
    "matrix_epipersepspectral1",
    "matrix_epipersepspectral2",
    "wsosinterpnonnegative1",
    "wsosinterpnonnegative2",
]

# instances to test with extended formulations option off
inst_noextend = String[
    "epinormeucl1",
    "epinormeucl2",
    "epipersquare1",
    "epipersquare2",
    "epinorminf1",
    "epinorminf2",
    "hypogeomean1",
    "hypogeomean2",
    "vector_epipersepspectral1",
    "vector_epipersepspectral2",
]

function runtests(oa_solver)
    @testset "solve conic" begin
        @testset "iterative" begin
            @testset "all" begin
                run_jump_tests(inst_all, true, true, true, oa_solver)
            end
            @testset "no ext" begin
                run_jump_tests(inst_noextend, false, true, true, oa_solver)
            end
        end
        @testset "one tree" begin
            @testset "all" begin
                run_jump_tests(inst_all, true, false, true, oa_solver)
            end
            @testset "no ext" begin
                run_jump_tests(inst_noextend, false, false, true, oa_solver)
            end
        end
    end
    @testset "sep only" begin
        @testset "iterative" begin
            @testset "all" begin
                run_jump_tests(inst_all, true, true, false, oa_solver)
            end
            @testset "no ext" begin
                run_jump_tests(inst_noextend, false, true, false, oa_solver)
            end
        end
        @testset "one tree" begin
            @testset "all" begin
                run_jump_tests(inst_all, true, false, false, oa_solver)
            end
            @testset "no ext" begin
                run_jump_tests(inst_noextend, false, false, false, oa_solver)
            end
        end
    end
    return
end

function run_jump_tests(
    insts::Vector{String},
    use_extended_form::Bool,
    use_iterative_method::Bool,
    solve_relax_subp::Bool,
    oa_solver,
)
    hypatia = MOI.OptimizerWithAttributes(
        Hypatia.Optimizer,
        MOI.Silent() => true,
        "near_factor" => 1000,
        "tol_feas" => 1e-10,
        "tol_infeas" => 1e-12,
        "tol_rel_opt" => 1e-9,
        "tol_abs_opt" => 1e-8,
        "tol_illposed" => 1e-9,
        "tol_slow" => 2e-2,
        "tol_inconsistent" => 1e-7,
    )

    opt = JuMP.optimizer_with_attributes(
        MOIPajarito.Optimizer,
        "verbose" => true,
        "oa_solver" => oa_solver,
        "conic_solver" => hypatia,
        "sep_solver" => hypatia,
        "use_extended_form" => use_extended_form,
        "use_iterative_method" => use_iterative_method,
        "solve_relaxation" => solve_relax_subp,
        "solve_subproblems" => solve_relax_subp,
        "iteration_limit" => 200,
        # "time_limit" => 120.0,
    )

    @testset "$inst" for inst in insts
        println()
        @info inst
        eval(Symbol(inst))(opt)
    end
    println()
    return
end

include("JuMP_instances.jl")

end
