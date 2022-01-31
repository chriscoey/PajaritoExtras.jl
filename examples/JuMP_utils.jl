#=
common code for JuMP examples
=#

const sparse_hypatia = MOI.OptimizerWithAttributes(
    Hypatia.Optimizer,
    MOI.Silent() => true,
    "near_factor" => 1000,
    "tol_feas" => 1e-10,
    "tol_rel_opt" => 1e-9,
    "tol_abs_opt" => 1e-8,
    "tol_illposed" => 1e-9,
    "tol_slow" => 2e-2,
    "tol_inconsistent" => 1e-7,
    "syssolver" => Solvers.SymIndefSparseSystemSolver{Float64}(),
    "init_use_indirect" => true,
    "preprocess" => false,
    "reduce" => false,
)

const dense_hypatia = MOI.OptimizerWithAttributes(
    Hypatia.Optimizer,
    MOI.Silent() => true,
    "near_factor" => 1000,
    "tol_feas" => 1e-10,
    "tol_rel_opt" => 1e-9,
    "tol_abs_opt" => 1e-8,
    "tol_illposed" => 1e-9,
    "tol_slow" => 2e-2,
    "tol_inconsistent" => 1e-7,
    "syssolver" => Solvers.QRCholDenseSystemSolver{Float64}(),
    "init_use_indirect" => false,
    "preprocess" => true,
    "reduce" => true,
)

# TODO experiment with best options and optimize path in Hypatia
const sep_hypatia = MOI.OptimizerWithAttributes(
    Hypatia.Optimizer,
    MOI.Silent() => true,
    "near_factor" => 1000,
    "tol_feas" => 1e-10,
    "tol_rel_opt" => 1e-9,
    "tol_abs_opt" => 1e-8,
    "tol_illposed" => 1e-9,
    "tol_slow" => 2e-2,
    "tol_inconsistent" => 1e-7,
    "syssolver" => Solvers.QRCholDenseSystemSolver{Float64}(),
    "init_use_indirect" => false,
    "preprocess" => true,
    "reduce" => true,
    "stepper" => Solvers.PredOrCentStepper{Float64}(),
)

function run_instance(
    ex_type::Type{<:ExampleInstance},
    inst_data::Tuple,
    inst_options::NamedTuple = NamedTuple();
    default_options::NamedTuple = NamedTuple(),
)
    println("setup example instance")
    Random.seed!(1)
    inst = ex_type(inst_data...)

    println("setup JuMP model")
    setup_time = @elapsed model = build(inst)
    # TODO get model_stats

    # set solver options
    println("set solver options")
    JuMP.set_optimizer(model, MOIPajarito.Optimizer)
    options = (; default_options..., conic_solver = dense_hypatia, inst_options...)
    for (option, value) in pairs(options)
        JuMP.set_optimizer_attribute(model, string(option), value)
    end
    JuMP.set_optimizer_attribute(model, "sep_solver", sep_hypatia)

    println("solve")
    flush(stdout)
    flush(stderr)
    optimize_time = @elapsed JuMP.optimize!(model)
    flush(stdout)
    flush(stderr)
    solve_stats = get_solve_stats(model)

    println("check")
    test_extra(inst, model)
    # TODO get solve_stats

    println("done")
    return (;
        # model_stats...,
        solve_stats...,
        setup_time,
        optimize_time,
    )
end

function get_solve_stats(model::JuMP.Model)
    solve_time = JuMP.solve_time(model)
    primal_obj = JuMP.objective_value(model)
    rel_gap = JuMP.relative_gap(model)
    status = string(JuMP.termination_status(model))

    opt = MOI.get(JuMP.backend(model), MOI.RawSolver())
    paj_ext = opt.use_extended_form
    paj_iter = opt.use_iterative_method

    iters_nodes = if paj_iter
        opt.num_iters
    else
        # JuMP.node_count(model) # TODO see https://github.com/jump-dev/Gurobi.jl/issues/444
        # Int(MOI.get(MOI.get(opt.oa_model, MOI.RawSolver()), MOI.NodeCount()))
        Int(MOI.get(opt.oa_opt, MOI.NodeCount()))
    end
    num_subp = length(opt.int_sols_cuts)
    num_cuts = opt.num_cuts

    solve_stats = (;
        paj_ext,
        paj_iter,
        status,
        solve_time,
        iters_nodes,
        num_subp,
        num_cuts,
        primal_obj,
        rel_gap,
    )
    return solve_stats
end

# fallback: just check optimal status
function test_extra(inst::ExampleInstance, model::JuMP.Model)
    @test JuMP.termination_status(model) == MOI.OPTIMAL
end
