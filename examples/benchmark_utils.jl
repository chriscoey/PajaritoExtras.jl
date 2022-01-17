#=
utilities for benchmark scripts
=#

function setup_benchmark_dataframe()
    perf = DataFrames.DataFrame(
        example = String[],
        inst_set = String[],
        inst_num = Int[],
        inst_data = Tuple[],
        solver = String[],
        solver_options = Tuple[],
        # nondefaults = Tuple[],
        # script_status = String[],
        # n = Int[],
        # p = Int[],
        # q = Int[],
        # num_int = Int[],
        # cone_types = Vector{String}[],
        # num_cones = Int[],
        # max_q = Int[],
        status = String[],
        solve_time = Float64[],
        # iters = Int[],
        # primal_obj = Float64[],
        # dual_obj = Float64[],
        # rel_obj_diff = Float64[],
        # x_viol = Float64[],
        # y_viol = Float64[],
        # setup_time = Float64[],
        # check_time = Float64[],
        total_time = Float64[],
    )
    # DataFrames.allowmissing!(perf, 11:DataFrames.ncol(perf))
    # DataFrames.allowmissing!(perf, [:solver_options, :nondefaults])
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
    new_default_options::NamedTuple,
    perf::DataFrames.DataFrame,
    results_path::Union{String, Nothing} = nothing,
)
    @show inst_subset
    for (inst_num, inst) in enumerate(inst_subset)
        # nondefaults = get_nondefaults(inst, ex_type_T)
        test_info = "inst $inst_num: $(inst[1])"
        # if !ismissing(nondefaults)
        #     test_info *= " $nondefaults"
        # end

        @testset "$test_info" begin
            println(test_info, " ...")

            total_time = @elapsed run_perf =
                run_instance(ex_type, inst..., default_options = new_default_options)

            new_perf = (;
                info_perf...,
                run_perf...,
                total_time,
                inst_num,
                :solver => "MOIPajarito",
                :inst_data => inst[1],
                # :nondefaults => nondefaults
            )
            write_perf(perf, results_path, new_perf)

            @printf("%8.2e seconds\n", total_time)
        end
    end
    return
end
