#=
utilities for benchmark scripts
=#

function setup_benchmark_dataframe()
    perf = DataFrames.DataFrame(
        example = String[],
        inst_set = String[],
        inst_num = Int[],
        inst_data = Tuple[],
        # n = Int[],
        # p = Int[],
        # q = Int[],
        # num_int = Int[],
        # cone_types = Vector{String}[],
        # num_cones = Int[],
        # max_q = Int[],
        paj_ext = Bool[],
        paj_iter = Bool[],
        status = String[],
        solve_time = Float64[],
        iters_nodes = Int[],
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

function run_instance_set(
    inst_subset::Vector,
    ex_type::Type{<:ExampleInstance},
    info_perf::NamedTuple,
    default_options::NamedTuple,
    perf::DataFrames.DataFrame,
    results_path::Union{String, Nothing} = nothing,
)
    for (inst_num, inst) in enumerate(inst_subset)
        test_info = "inst $inst_num: $(inst[1])"

        @testset "$test_info" begin
            println(test_info, " ...")

            total_time = @elapsed run_perf =
                run_instance(ex_type, inst..., default_options = default_options)

            new_perf =
                (; info_perf..., run_perf..., total_time, inst_num, :inst_data => inst[1])
            write_perf(perf, results_path, new_perf)

            @printf("%8.2e seconds\n", total_time)
        end
    end
    return
end

function run_examples(
    inst_sets::Vector{String},
    default_options::NamedTuple,
    results_path::Union{String, Nothing} = nothing,
)
    perf = setup_benchmark_dataframe()
    isnothing(results_path) || CSV.write(results_path, perf)

    @testset "examples tests" begin
        @testset "$ex" for (ex, (ex_type, ex_insts)) in get_test_instances()
            @testset "$inst_set" for inst_set in inst_sets
                haskey(ex_insts, inst_set) || continue
                inst_subset = ex_insts[inst_set]
                isempty(inst_subset) && continue

                info_perf = (; inst_set, :example => ex)
                str = "$ex $inst_set"
                println("\nstarting $str tests")
                @testset "$str" begin
                    run_instance_set(
                        inst_subset,
                        ex_type,
                        info_perf,
                        default_options,
                        perf,
                        results_path,
                    )
                end
            end
        end
    end

    println()
    DataFrames.show(perf, allrows = true, allcols = true)
    println()
    return perf
end
