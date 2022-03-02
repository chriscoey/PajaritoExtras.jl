#=
design a modular device using parts from a catalog:
- for each module we pick one corresponding part from the catalog, satisfying a constraint
associated with the part
- these constraints involve certain global design variables (e.g. weight, size,
lifetime, battery power)
- we minimize the cost of the parts plus a linear function of the global design variables

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

    # generate parameters
    c = 10 * rand(m, jmax)
    # d = vec(sprandn(n, 1, 0.3))
    d = zeros(n)

    # build model
    model = JuMP.Model()
    y = JuMP.@variable(model, [1:n])
    x = JuMP.@variable(model, [1:m, 1:jmax], Bin)
    # JuMP.@objective(model, Min, dot(c, x) + dot(d, y)) # TODO
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

    return model
end

function test_extra(inst::NonconvexRelax, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    x = JuMP.value.(model.ext[:x])
    @test x ≈ round.(Int, x) atol = tol rtol = tol
    @test all(-tol .<= x .<= 1 + tol)
    y = JuMP.value.(model.ext[:y])
    z = JuMP.value.(model.ext[:z])
    # TODO

    return
end
