#=
analyze benchmark results
~/julia/julia analyze.jl
=#

using Printf
using CSV
using DataFrames

include(joinpath(@__DIR__, "../examples/Examples.jl"))
import Main.Examples

bench_file = joinpath(@__DIR__, "raw", "bench.csv")

output_dir = mkpath(joinpath(@__DIR__, "analysis"))
tex_dir = mkpath(joinpath(output_dir, "tex"))
stats_dir = mkpath(joinpath(output_dir, "stats"))
csv_dir = mkpath(joinpath(output_dir, "csvs"))

# uncomment examples to process
examples_params = [
    # PSD:
    "completablepsd" => ([:d], [1], [:nat, :ext], [:iters]),
    # WSOS:
    "polyfacilitylocation" => ([:n], [1], [:nat, :ext], [:iters]),
    "polyregression" => ([:m], [3], [:nat, :ext], [:iters]),
    "twostagestochastic" => ([:halfdeg], [2], [:nat, :ext], [:iters]),
    # spectral norm:
    "matrixcompletion" => ([:nrow], [3], [:nat, :ext], [:iters]),
    "matrixdecomposition" => ([:n], [2], [:nat, :ext], [:iters]),
    "matrixregression" => ([:p], [2], [:nat, :ext], [:iters]),
    # spectral function:
    "inversecovariance" => ([:d], [1], [:nat, :ext], [:iters]),
    "vectorregression" => ([:n], [1], [:nat, :noext], [:iters]),
    "knapsack" => ([:n], [1], [:cont_geo, :cont_noext_geo], [:iters], "cont_geo"),
    "knapsack" => ([:n], [1], [:cont_log, :cont_noext_log], [:iters], "cont_log"),
    "knapsack" => ([:n], [1], [:cont_inv, :cont_noext_inv], [:iters], "cont_inv"),
    "knapsack" => ([:n], [1], [:nat_geo, :noext_geo, :ext_geo], [:iters], "geo"),
    "experimentdesign" => ([:d], [1], [:nat_rtdet, :ext_rtdet], [:iters], "rtdet"),
    "experimentdesign" => ([:d], [1], [:nat_entr, :ext_entr], [:iters], "entr"),
    # MIP formulations:
    "ballpacking" => ([:num_pts], [4], [:sos2, :logib, :cc], [:nodes]),
    "modulardesign" => ([:n], [3], [:nat, :noext], [:iters], "convex"),
    "modulardesign" => ([:n], [3], [:sos2, :logib, :cc], [:iters], "nonconvex"),
]

function analyze()
    col_types = Dict(:example => String, :inst_set => String, :status => String)
    all_df = CSV.read(bench_file, DataFrame, types = col_types)
    @show all_df

    for (ex_name, ex_params) in examples_params
        println()
        @info("starting $ex_name with params: $ex_params")

        sets = string.(ex_params[3])
        ex_df = filter(
            t -> (t.inst_num > 1) && (t.example == ex_name) && (t.inst_set in sets),
            all_df,
        )
        if isempty(ex_df)
            @info("no data for $ex_name")
            continue
        end

        if length(ex_params) == 5
            ex_name *= "_" * ex_params[5]
        end

        ex_df_wide = make_wide_csv(ex_df, ex_name, ex_params)
        make_table_tex(ex_name, ex_params, ex_df_wide)
        make_plot_csv(ex_name, ex_params, ex_df_wide)

        @info("finished")
    end

    println()
    @info("finished all")
end

status_map = Dict(
    "OPTIMAL" => "co",
    "TIME_LIMIT" => "tl",
    "ITERATION_LIMIT" => "il",
    "OTHER_ERROR" => "er",
)

function relative_tol_satisfied(a::T, b::T, tol::T = 1e-3) where {T <: Real}
    return (abs(a - b) / (1 + max(abs(a), abs(b))) < tol)
end

function make_wide_csv(ex_df, ex_name, ex_params)
    @info("making wide csv for $ex_name")
    inst_keys = ex_params[1]
    str(x) = split(x[2:(end - 1)], ", ")

    replace!(ex_df.status, status_map...)

    for (name, pos) in zip(inst_keys, ex_params[2])
        transform!(
            ex_df,
            :inst_data => ByRow(x -> string(eval(Meta.parse(x))[pos])) => name,
        )
    end

    # check objectives if solver claims optimality
    show_cols = [:inst_set, :inst_num, :status, :primal_obj]
    for group_df in groupby(ex_df, inst_keys)
        # check all pairs of verified converged results
        co_idxs = findall(group_df[:, :status] .== "co")
        (length(co_idxs) < 2) && continue

        first_optval = group_df[co_idxs[1], :primal_obj]
        other_optvals = group_df[co_idxs[2:end], :primal_obj]
        if !all(relative_tol_satisfied.(other_optvals, first_optval))
            dat = group_df[!, :inst_data][1]
            println("$ex_name $dat primal optimal objective values disagree:")
            show(group_df[!, show_cols], eltypes = false)
            println("\n")
        end
    end

    unstacked = [
        unstack(
            ex_df,
            inst_keys,
            :inst_set,
            v,
            renamecols = x -> Symbol(v, :_, x),
            allowduplicates = true,
        ) for v in [:status, :solve_time, :nodes, :iters, :num_subp, :num_cuts]
    ]
    ex_df_wide = outerjoin(unstacked..., on = inst_keys)

    CSV.write(joinpath(stats_dir, ex_name * "_wide.csv"), ex_df_wide)

    return ex_df_wide
end

process_entry(::Missing) = "\$\\ast\$"

function process_entry(x::Float64)
    isnan(x) && return "\$\\ast\$"
    @assert x > 0
    if x < 10
        return @sprintf("%.1f", x)
    else
        return @sprintf("%.0f.", x)
    end
end

process_entry(x) = string(x)

function make_table_tex(ex_name, ex_params, ex_df_wide)
    @info("making table tex for $ex_name")
    inst_sets = string.(ex_params[3])
    inst_keys = ex_params[1]
    col_names = ex_params[4]

    sep = " & "
    ex_tex = open(joinpath(tex_dir, ex_name * "_table.tex"), "w")

    header_str = "% " * prod(string(s) * sep for s in inst_keys)
    for s in inst_sets
        header_str *= s * prod(sep * string(c) for c in col_names) * sep * "time" * sep
    end
    println(ex_tex, header_str)

    for row in eachrow(ex_df_wide)
        row_str = process_entry(row[1])
        for i in 2:length(inst_keys)
            row_str *= sep * process_entry(row[i])
        end
        for inst_set in inst_sets
            row_str *= sep * process_entry(row[Symbol(:status_, inst_set)])
            for col_name in col_names
                row_str *= sep * process_entry(row[Symbol(col_name, :_, inst_set)])
            end
            row_str *= sep * process_entry(row[Symbol(:solve_time_, inst_set)])
        end
        row_str *= " \\\\"
        println(ex_tex, row_str)
    end

    close(ex_tex)

    return
end

function transform_plot_cols(ex_df_wide, inst_set)
    old_cols = Symbol.([:status_, :solve_time_], inst_set)
    # old_cols = Symbol.([:status_, :iters_], inst_set)
    return transform!(
        ex_df_wide,
        old_cols =>
            ByRow((x, y) -> ((!ismissing(x) && x == "co") ? y : missing)) => inst_set,
    )
end

function make_plot_csv(ex_name, ex_params, ex_df_wide)
    @info("making plot csv for $ex_name")
    inst_sets = string.(ex_params[3])
    inst_keys = ex_params[1]

    for inst_set in inst_sets
        transform_plot_cols(ex_df_wide, inst_set)
    end

    plot_file_start = joinpath(csv_dir, ex_name * "_plot")

    axis_name = last(inst_keys)
    if length(inst_keys) == 1
        success_df = select(ex_df_wide, axis_name, inst_sets...)
        CSV.write(plot_file_start * ".csv", success_df)
    else
        group_name = first(inst_keys)
        success_df = select(ex_df_wide, axis_name, group_name, inst_sets...)
        for (group_id, group_df) in pairs(groupby(success_df, group_name))
            CSV.write(
                plot_file_start * "_$(group_id[1]).csv",
                select(group_df, Not(group_name)),
            )
        end
    end

    return
end

analyze();
