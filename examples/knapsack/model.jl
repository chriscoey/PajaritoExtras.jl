#=
simple integer knapsack problem with convex objective

minₓ f(c ⊙ x) :
b' x ≤ B
x ∈ ℤⁿ

where c > 0, b > 0, B > 0, and f : ℝⁿ₊ → ℝ is a convex spectral function
on the nonnegative domain, which forces x ≥ 0
=#

struct Knapsack <: ExampleInstance
    n::Int
    f::VecSpecExt
    relax_int::Bool # relax integrality
end

function build(inst::Knapsack)
    n = inst.n

    # generate parameters
    c = 2 * rand(n)
    b = rand(n)
    B = 3 * sum(b)

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, x[1:n])
    JuMP.@constraint(model, dot(b, x) <= B)

    if !inst.relax_int
        JuMP.set_integer.(x)
    end

    JuMP.@variable(model, y)
    JuMP.@objective(model, Min, y)
    add_homog_spectral(inst.f, n, vcat(y, c .* x), model)

    # save for use in tests
    model.ext[:x] = x
    model.ext[:c] = c

    return model
end

function test_extra(inst::Knapsack, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat in OPT_OR_LIMIT
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    x = JuMP.value.(model.ext[:x])
    if !inst.relax_int
        @test x ≈ round.(Int, x) atol = tol rtol = tol
    end
    # check conic feasibility
    @test all(>(-tol), x)
    c = JuMP.value.(model.ext[:c])
    obj_result = get_val(pos_only(c .* x), inst.f)
    @test JuMP.objective_value(model) ≈ obj_result atol = tol rtol = tol
    return
end
