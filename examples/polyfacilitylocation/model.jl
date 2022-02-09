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
dⱼ  ∈ P₊  demand at j
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
    use_nat::Bool # use WSOS cone formulation, else SDP formulation
end

function build(inst::PolyFacilityLocation)
    # setup polynomial interpolation
    dom = PolyUtils.BoxDomain{Float64}([0.0], [1.0])
    (U, _, Ps, _, w) = PolyUtils.interpolate(dom, inst.deg, get_quadr = true)

    # generate parameter values
    (n, m) = (inst.n, inst.m)
    f = rand(n)
    c = rand(n, m)
    d = make_nonneg_polys(m, Ps)
    u = make_nonneg_polys(n, Ps)
    # make feasible by ensuring sum(u) > 2 * sum(d)
    excess = find_poly_min(sum(u) - 2 * sum(d), Ps)
    incr = 0.1 - excess / n
    for i in 1:n
        u[i] .+= incr + 5 * rand() + 1
    end
    excess = find_poly_min(sum(u) - 2 * sum(d), Ps)
    @assert excess > 0

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

    # save for use in tests
    model.ext[:x] = x

    return model
end

function test_extra(inst::PolyFacilityLocation, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    # check integer feasibility
    tol = eps()^0.2
    x = JuMP.value.(model.ext[:x])
    @test x ≈ round.(Int, x) atol = tol rtol = tol
    @test all(-tol .<= x .<= 1 + tol)

    # TODO check nonnegativity of wsos constraint functions
    # can just solve the separation problem to check feasibility
    # maybe also plot to check visually
    return
end
