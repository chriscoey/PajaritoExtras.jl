#=
utilities for benchmark scripts
=#

function setup_benchmark_dataframe()
    perf = DataFrames.DataFrame(
        example = String[],
        inst_set = String[],
        inst_num = Int[],
        inst_data = Tuple[],
        paj_ext = Bool[],
        paj_iter = Bool[],
        status = String[],
        solve_time = Float64[],
        nodes = Int[],
        iters = Int[],
        num_subp = Int[],
        num_cuts = Int[],
        primal_obj = Float64[],
        rel_gap = Float64[],
        setup_time = Float64[],
        optimize_time = Float64[],
        total_time = Float64[],
    )
    # DataFrames.allowmissing!(perf, 11:DataFrames.ncol(perf))
    return perf
end

function write_perf(
    perf::DataFrames.DataFrame,
    results_path::Union{String, Nothing},
    new_perf::NamedTuple,
)
    push!(perf, new_perf, cols = :subset)
    if !isnothing(results_path)
        CSV.write(
            results_path,
            perf[end:end, :],
            transform = (col, val) -> something(val, missing),
            append = true,
        )
    end
    return
end

const limit_statuses = ["TIME_LIMIT", "ITERATION_LIMIT"]

function run_instance_set(
    inst_subset::Vector{<:Tuple},
    ex_type::Type{<:ExampleInstance},
    info_perf::NamedTuple,
    default_options::NamedTuple,
    perf::DataFrames.DataFrame,
    results_path::Union{String, Nothing},
    skip_limit::Bool,
)
    @testset "inst $inst_num: $(inst[1])" for (inst_num, inst) in enumerate(inst_subset)
        println("inst $inst_num: $(inst[1]) ...")

        total_time = @elapsed run_perf =
            run_instance(ex_type, inst..., default_options = default_options)
        @printf("%8.2e seconds\n", total_time)

        p = (; info_perf..., run_perf..., total_time, inst_num, :inst_data => inst[1])
        write_perf(perf, results_path, p)

        if skip_limit && run_perf.status in limit_statuses
            println("skipping remaining instances due to limit status\n")
            break
        end
        println()
    end
    return
end

function run_examples(
    examples::Vector{String},
    inst_sets::Vector{String},
    default_options::NamedTuple,
    results_path::Union{String, Nothing},
    skip_limit::Bool, # skip larger instances after hitting a solver limit
)
    # setup dataframe
    perf = setup_benchmark_dataframe()
    isnothing(results_path) || CSV.write(results_path, perf)

    # run examples
    @testset "$ex" for ex in examples
        # load instances
        (ex_type, ex_insts) = include(joinpath(@__DIR__, ex, "instances.jl"))

        # run instance sets
        for inst_set in inst_sets
            haskey(ex_insts, inst_set) || continue

            @testset "$inst_set" begin
                inst_subset = ex_insts[inst_set]
                isempty(inst_subset) && continue

                info_perf = (; inst_set, :example => ex)
                str = "$ex $inst_set"
                println("\nstarting $str tests\n")
                @testset "$str" begin
                    run_instance_set(
                        inst_subset,
                        ex_type,
                        info_perf,
                        default_options,
                        perf,
                        results_path,
                        skip_limit,
                    )
                end
            end
        end
    end

    println()
    DataFrames.show(perf, allrows = true, allcols = true)
    println("\n")
    return perf
end
