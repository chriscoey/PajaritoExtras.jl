#=
capacitated facility location problem with polynomial flows over continuous time

time is t ∈ [0, 1]
P are polynomials in t of max degree deg ≥ 1, and P₊ ⊂ P are nonnegative on [0, 1]

sets:
i ∈ 1..n   facilities
j ∈ 1..m   customers

parameters:
fᵢ  ∈ ℝ₊  fixed cost of opening i
cᵢⱼ ∈ ℝ₊  cost (per unit of flow) to ship from i to j
dⱼ  ∈ P   demand at j
uᵢ  ∈ P₊  maximum output of i

variables:
xᵢ  ∈ {0,1}  whether i is opened
yᵢⱼ ∈ P₊     flow from i to j

objective:
min ∑ᵢⱼ cᵢⱼ ∫ yᵢⱼ + ∑ᵢ fᵢ xᵢ   minimize total cost

constraints:
∑ⱼ yᵢⱼ ≤ uᵢ xᵢ, ∀i   limit production at i
∑ᵢ yᵢⱼ ≥ dⱼ, ∀j      meet demand at j
=#

struct PolyFacilityLocation <: ExampleInstance
    n::Int
    m::Int
    deg::Int
    use_nat::Bool
end

function build(inst::PolyFacilityLocation)
    # setup polynomial interpolation
    dom = PolyUtils.BoxDomain{Float64}([0.0], [1.0])
    (U, pts, Ps, _, w) = PolyUtils.interpolate(dom, inst.deg, get_quadr = true)

    # generate parameter values
    (n, m) = (inst.n, inst.m)
    f = rand(n)
    c = rand(n, m)
    L = 1 + div(inst.deg, 2)
    P1L = Ps[1][:, 1:L]
    d = [P1L * rand(L) for _ in 1:m]
    # u = [ones(U) for _ in 1:n]
    u = [5 * ones(U) for _ in 1:n]

    # TODO to generate nonnegative polys, generate random poly with largest degree even and positive coef, and add minimum value of poly
    # or generate from random PSD matrix with linear algebra

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, x[1:n], Bin)
    JuMP.@variable(model, y[1:n, 1:m, 1:U])

    JuMP.@objective(
        model,
        Min,
        sum(JuMP.dot(c[i, j] * w, y[i, j, :]) for j in 1:m, i in 1:n) + JuMP.dot(f, x)
    )

    add_wsos(aff) = (inst.use_nat ? add_wsos_nat : add_wsos_ext)(Ps, aff, model)
    for i in 1:n, j in 1:m
        add_wsos(y[i, j, :])
    end
    for i in 1:n
        add_wsos(JuMP.@expression(model, u[i] * x[i] - sum(y[i, j, :] for j in 1:m)))
    end
    for j in 1:m
        add_wsos(JuMP.@expression(model, sum(y[i, j, :] for i in 1:n) - d[j]))
    end

    return model
end

function test_extra(inst::PolyFacilityLocation, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    # TODO data generation doesn't guarantee feasibility
    # @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    # TODO check nonnegativity of wsos constraint functions
    # can just solve the separation problem to check feasibility
    # maybe also plot to check visually
    return
end
