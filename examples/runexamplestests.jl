#=
run examples tests from the examples folder
=#

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

import MOIPajarito
default_options = (;
    verbose = true,
    oa_solver = gurobi,
    conic_solver = hypatia,
    use_extended_form = true,
    use_iterative_method = true,
    # "debug_cuts" => use_iterative_method,
    # "iteration_limit" => 30,
    # "time_limit" => 120.0,
)

# instance sets to run and corresponding time limits (seconds)
inst_sets = [
    ("minimal", 60),
    # ("various", 120),
]

using Test
import DataFrames
include(joinpath(@__DIR__, "Examples.jl"))

perf = Examples.setup_benchmark_dataframe()

@testset "examples tests" begin
    @testset "$ex" for (ex, (ex_type, ex_insts)) in Examples.get_test_instances()
        @testset "$inst_set, $time_limit" for (inst_set, time_limit) in inst_sets
            haskey(ex_insts, inst_set) || continue
            inst_subset = ex_insts[inst_set]
            isempty(inst_subset) && continue

            info_perf = (; inst_set, :example => ex)
            new_default_options = (; default_options..., time_limit = time_limit)

            str = "$ex $inst_set"
            println("\nstarting $str tests")
            @testset "$str" begin
                Examples.run_instance_set(
                    inst_subset,
                    ex_type,
                    info_perf,
                    new_default_options,
                    perf,
                )
            end
        end
    end

    println("\n")
    DataFrames.show(perf, allrows = true, allcols = true)
    println("\n")
end;
