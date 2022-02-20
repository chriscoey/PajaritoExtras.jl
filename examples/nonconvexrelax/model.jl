#=
In the original nonconvex problem, there are disjunctions over nonconvex constraints,
where the nonconvex constraints are like cone boundary or outside-cone conditions, with cones
for which we have extended formulations involving only 3-dimensional cones. We relax this
disjunctive nonconvex problem to a mixed-integer convex problem using efficient `effectively
univariate` piecewise linear formulations of the 3-dimensional cone boundaries in the
extended formulation.

We can think of designing a modular device using parts from a catalog. For each module we
pick one corresponding part from the catalog, satisfying a constraint associated with the
part. These constraints involve certain global design variables (e.g. weight, size,
lifetime, battery power). We minimize the cost of the parts plus a linear function of the
global design variables. We could handle more complex constraints linking parts/modules and
design variables, but we keep this model simple.

sets:
i ∈ 1..m          modules
j ∈ J[i], i ∈ m   parts for module i
k ∈ 1..n          global design variables

variables:
xᵢⱼ ∈ {0,1}   whether part j ∈ J[i] is selected for module i
y   ∈ ℝⁿ      global design characteristics

objective:
min ∑ᵢⱼ cᵢⱼ xᵢⱼ + ∑ₖ dₖ yₖ   where c ≥ 0

constraints:
∑ⱼ xᵢⱼ = 1, ∀i                 disjunction for module i
(xᵢⱼ = 1) ⟹ (y ∈ Sᵢⱼ), ∀i,j   constraint for part i,j
where Sᵢⱼ ⊂ ℝⁿ

use conic copies-of-variables formulation for the disjunctions
=#

struct NonconvexRelax <: ExampleInstance
    m::Int # number of modules
    jmax::Int # number of part options per module
    n::Int # number of design variables
    use_nonconvex::Bool # use relaxation of nonconvex constraints, else convex constraints
    pwl::PWLSOS2 # type of piecewise linear SOS2 formulation to use
    num_pts::Int # number of points per piecewise linearization
end

NonconvexRelax(m::Int, jmax::Int, n::Int) = NonconvexRelax(m, jmax, n, false, SOS2(), 0)

function build(inst::NonconvexRelax)
    (m, jmax, n, num_pts) = (inst.m, inst.jmax, inst.n, inst.num_pts)
    @assert m >= 1
    @assert jmax >= 1
    @assert n >= 2

    # generate random solution to ensure feasibility
    x0 = rand(Bool, m, jmax)
    y0 = 2 * rand(n)
    @show y0

    # generate parameters
    c = 10 * rand(m, jmax)
    # d = vec(sprandn(n, 1, 0.3))
    d = zeros(n)

    # build model
    model = JuMP.Model()
    y = JuMP.@variable(model, [1:n])
    x = JuMP.@variable(model, [1:m, 1:jmax], Bin)
    # JuMP.@objective(model, Min, dot(c, x) + dot(d, y))
    JuMP.@objective(model, Min, dot(c, x) - sum(y))

    # disjunctive copies-of-variables formulation
    # x[i, j] switches on constraint for part j of module i
    z = JuMP.@variable(model, [1:m, 1:jmax, 1:n])
    JuMP.@constraint(model, [i in 1:m], sum(x[i, :]) == 1)
    JuMP.@constraint(model, [i in 1:m], sum(z[i, j, :] for j in 1:jmax) .== y)
    for i in 1:m, j in 1:jmax
        f = VecPower12(2.0) # TODO
        f0 = get_val(y0, f)
        x_ij = x[i, j]
        z_ij = z[i, j, :]
        aff = vcat(f0 * x_ij, x_ij, z_ij)
        if inst.use_nonconvex
            add_nonconvex(f, aff, model, true, false, 0.0, 2.0, inst.pwl, num_pts)
        else
            add_spectral(f, length(z_ij), aff, model)
        end
    end

    # save for use in tests
    model.ext[:x] = x
    model.ext[:y] = y
    model.ext[:z] = z

    # JuMP.set_optimizer(model, Gurobi.Optimizer)
    # JuMP.optimize!(model)

    # @show JuMP.value.(x)
    # @show JuMP.value.(y)
    # @show sum(abs2, JuMP.value.(y)) / 2
    # @show JuMP.value.(z)
    # error()

    return model
end

function test_extra(inst::NonconvexRelax, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    # check integer feasibility
    tol = eps()^0.2
    x = JuMP.value.(model.ext[:x])
    @show x
    @test x ≈ round.(Int, x) atol = tol rtol = tol
    @test all(-tol .<= x .<= 1 + tol)

    y = JuMP.value.(model.ext[:y])
    z = JuMP.value.(model.ext[:z])
    @show y
    @show z
    # TODO

    return
end
