#=
two stage stochastic farming problem

sets:
i ∈ 1..n  crops

parameters:
̄̄m  ∈ ℝ₊  number of available plots
aᵢ ∈ ℝ₊  fixed cost per farm plot for i
bᵢ ∈ ℝ₊  purchase cost for i
cᵢ ∈ ℝ₊  selling price for i
dᵢ ∈ ℝ₊  contractual demand for i

uncertain parameters:
ξᵢ ∈ ℝ₊  crop yield per farm plot for i
where ξ ~ Ξ

variables:
xᵢ ∈ ℤ₊  number of farm plots for i (1st stage)
yᵢ ∈ P₊  amount of crop i to purchase (2nd stage)
zᵢ ∈ P₊  amount of crop i to sell (2nd stage)
where P₊ are nonnegative polynomials in ξ

objective:
min ∑ᵢ aᵢ xᵢ + 𝔼(ξ ~ Ξ)[∑ᵢ(bᵢ yᵢ - cᵢ zᵢ)]  minimize total expected profit

constraints:
∑ᵢ xᵢ ≤ m  available plots
xᵢ ξᵢ + yᵢ - zᵢ - dᵢ = 0, ∀i  balance inputs and outputs for i

optional linking constraints:
f - ∑ᵢ (yᵢ + zᵢ) ∈ P₊  capacity limit for all crops (f ∈ ℝ₊ is a parameter)

TODO
- CVaR objective or other
- if unlinked, variables can be partitioned when x is fixed, so solve separate subproblems
=#

struct TwoStageStochastic <: ExampleInstance
    n::Int
    deg::Int
    linked::Bool # whether stage two decisions are linked or independent
    use_nat::Bool # use WSOS cone formulation, else SDP formulation
end

function build(inst::TwoStageStochastic)
    n = inst.n
    m = 2n

    # generate parameter values
    a = rand(n)
    b = 1 .+ rand(n)
    c = b - rand(n)
    d = 1 .+ rand(n)
    f = sum(d) + 0.2 * n * rand()

    # setup polynomial interpolation
    n0 = (inst.linked ? n : 1)
    dom = PolyUtils.BoxDomain{Float64}(zeros(n0), ones(n0))
    (U, pts, Ps, _, w) = PolyUtils.interpolate(dom, inst.deg, get_quadr = true)
    w .*= (inst.linked ? 1 : 2^(n - 1))

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, x[1:n], Int)
    JuMP.@variable(model, y[1:n, 1:U])
    JuMP.@variable(model, z[1:n, 1:U])

    JuMP.@objective(
        model,
        Min,
        JuMP.dot(a, x) + JuMP.dot(w, sum(b[i] * y[i, :] - c[i] * z[i, :] for i in 1:n))
    )

    JuMP.@constraint(model, x .>= 0)
    JuMP.@constraint(model, sum(x) <= m)

    add_wsos(aff) = (inst.use_nat ? add_wsos_nat : add_wsos_ext)(Ps, aff, model)

    for i in 1:n
        add_wsos(y[i, :])
        add_wsos(z[i, :])
        pts_i = (inst.linked ? pts[:, i] : vec(pts))
        JuMP.@constraint(model, pts_i * x[i] + y[i, :] - z[i, :] .- d[i] .== 0)
    end

    if inst.linked
        fyz = JuMP.@expression(model, f .- sum(y[i, :] + z[i, :] for i in 1:n))
        add_wsos(fyz)
    end

    # save for use in tests
    model.ext[:x] = x
    model.ext[:y] = y
    model.ext[:z] = z
    if inst.linked
        model.ext[:fyz] = fyz
    end
    model.ext[:Ps] = Ps

    return model
end

function test_extra(inst::TwoStageStochastic, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    x = JuMP.value.(model.ext[:x])
    @test x ≈ round.(Int, x) atol = tol rtol = tol
    # check WSOS feasibility
    Ps = model.ext[:Ps]
    y = JuMP.value.(model.ext[:y])
    z = JuMP.value.(model.ext[:z])
    @test all(check_nonneg(y[i, :] .+ tol, Ps) for i in 1:(inst.n))
    @test all(check_nonneg(z[i, :] .+ tol, Ps) for i in 1:(inst.n))
    if inst.linked
        fyz = JuMP.value.(model.ext[:fyz])
        @test check_nonneg(fyz .+ tol, Ps)
    end
    return
end
