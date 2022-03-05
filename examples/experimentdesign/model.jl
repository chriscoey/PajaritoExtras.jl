#=
choose the frequency of experiments to minimize a given convex spectral function
of the information matrix and satisfy an experiment budget constraint

adapted from Boyd and Vandenberghe, "Convex Optimization", section 7.5

min f(V × Diagonal(x) × V')
subject to:
x ≥ 0, Integer
∑x = k

where k = 2d, variable x ∈ ℝᵏ is the frequency of each experiment, k is the
number of experiments to run, the columns of V ∈ ℝ^(d × k) correspond to each
experiment, and f is a convex spectral function
=#

struct ExperimentDesign <: ExampleInstance
    d::Int
    f::MatSpecExt # formulation specifier
end

function build(inst::ExperimentDesign)
    d = inst.d
    @assert d >= 1
    @assert is_domain_pos(inst.f)
    k = 2 * d

    # generate data
    V = randn(d, k)
    V .*= sqrt(d / sum(eigvals(Symmetric(V * V'))))

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, x[1:k] >= 0)
    JuMP.set_integer.(x)
    JuMP.@constraint(model, sum(x) == k)

    # vectorized information matrix
    Q_vec = [
        JuMP.@expression(
            model,
            (i == j ? 1.0 : rt2) * sum(V[i, k] * x[k] * V[j, k] for k in 1:k)
        ) for i in 1:d for j in 1:i
    ]

    # convex objective
    JuMP.@variable(model, epi)
    JuMP.@objective(model, Min, epi)
    add_homog_spectral(inst.f, d, vcat(1.0 * epi, Q_vec), model)

    # save for use in tests
    model.ext[:V] = V
    model.ext[:x] = x

    return model
end

function test_extra(inst::ExperimentDesign, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    x = JuMP.value.(model.ext[:x])
    @test x ≈ round.(Int, x) atol = tol rtol = tol
    # check objective
    V = model.ext[:V]
    λ = eigvals(Symmetric(V * Diagonal(x) * V', :U))
    @test minimum(λ) >= -tol
    obj_result = get_val(pos_only(λ), inst.f)
    @test JuMP.objective_value(model) ≈ obj_result atol = tol rtol = tol
    return
end
