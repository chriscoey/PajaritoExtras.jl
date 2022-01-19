#=
run examples tests from the examples folder
=#

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
    oa_solver = gurobi,
    use_extended_form = true,
    use_iterative_method = true,
    # debug_cuts = use_iterative_method,
    # iteration_limit = 30,
    # time_limit = 120.0,
)

# instance sets to run
inst_sets = [
    "minimal",
    # "various",
]

using Test
import DataFrames
include(joinpath(@__DIR__, "Examples.jl"))

perf = Examples.setup_benchmark_dataframe()

@testset "examples tests" begin
    @testset "$ex" for (ex, (ex_type, ex_insts)) in Examples.get_test_instances()
        @testset "$inst_set" for inst_set in inst_sets
            haskey(ex_insts, inst_set) || continue
            inst_subset = ex_insts[inst_set]
            isempty(inst_subset) && continue

            info_perf = (; inst_set, :example => ex)
            str = "$ex $inst_set"
            println("\nstarting $str tests")
            @testset "$str" begin
                Examples.run_instance_set(
                    inst_subset,
                    ex_type,
                    info_perf,
                    default_options,
                    perf,
                )
            end
        end
    end

    println("\n")
    DataFrames.show(perf, allrows = true, allcols = true)
    println("\n")
end;
