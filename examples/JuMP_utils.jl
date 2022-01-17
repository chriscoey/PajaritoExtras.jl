#=
common code for JuMP examples
=#

# fallback: just check optimal status
function test_extra(inst::ExampleInstance, model::JuMP.Model)
    @test JuMP.termination_status(model) == MOI.OPTIMAL
end

function run_instance(
    ex_type::Type{<:ExampleInstance},
    inst_data::Tuple,
    inst_options::NamedTuple = NamedTuple(),
    solver_type = MOIPajarito.Optimizer;
    default_options::NamedTuple = NamedTuple(),
)
    new_options = merge(default_options, inst_options)

    println("setup model")
    setup_time = @elapsed (model, model_stats) =
        setup_model(ex_type, inst_data, new_options, solver_type)

    println("solve and check")
    check_time = @elapsed solve_stats = solve_check(model)

    return (;
        # model_stats...,
        solve_stats...,
        setup_time,
        check_time,
        # :script_status => "Success",
    )
end

function setup_model(
    ex_type::Type{<:ExampleInstance},
    inst_data::Tuple,
    solver_options::NamedTuple,
    solver_type;
    rseed::Int = 1,
)
    # setup example instance and JuMP model
    Random.seed!(rseed)
    inst = ex_type(inst_data...)
    model = build(inst)

    # backend = JuMP.backend(model)
    # MOIU.reset_optimizer(backend, opt)
    # MOIU.attach_optimizer(backend)
    # MOIU.attach_optimizer(backend.optimizer.model)

    JuMP.set_optimizer(model, solver_type)
    for (option, value) in pairs(solver_options)
        JuMP.set_optimizer_attribute(model, string(option), value)
    end

    # return (model, get_model_stats(hyp_model))
    return (model, ())
end

function solve_check(model::JuMP.Model)
    JuMP.optimize!(model)
    flush(stdout)
    flush(stderr)

    # test_extra(inst, model)

    solve_time = JuMP.solve_time(model)
    # iters = MOI.get(model, MOI.BarrierIterations())
    primal_obj = JuMP.objective_value(model)
    # dual_obj = JuMP.dual_objective_value(model)
    status = JuMP.termination_status(model)

    solve_stats = (;
        status,
        solve_time,
        # iters,
        primal_obj,
        # rel_obj_diff,
    )
    return solve_stats
end
