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

# TODO move
import LinearAlgebra
import JuMP
import Hypatia
import Hypatia.PolyUtils
import MOIPajarito
import PajaritoExtras

oa_solver = gurobi
conic_solver = hypatia
opt = JuMP.optimizer_with_attributes(
    MOIPajarito.Optimizer,
    "verbose" => true,
    # "verbose" => false,
    "oa_solver" => oa_solver,
    "conic_solver" => conic_solver,
    # "use_extended_form" => use_extended_form,
    # "use_iterative_method" => use_iterative_method,
    # "debug_cuts" => use_iterative_method,
    # "iteration_limit" => 30,
    # "time_limit" => 120.0,
)

struct PolyFacilityLocation #<: ExampleInstanceJuMP{Float64}
    n::Int
    m::Int
    deg::Int
end

function build(inst::PolyFacilityLocation)
    # setup polynomial interpolation
    dom = PolyUtils.BoxDomain{Float64}([0.0], [1.0])
    (U, pts, Ps, _, w) = PolyUtils.interpolate(dom, inst.deg, get_quadr = true)
    K = Hypatia.WSOSInterpNonnegativeCone{Float64, Float64}(U, Ps)
    @show pts

    # generate parameter values
    (n, m) = (inst.n, inst.m)
    f = rand(n)
    c = rand(n, m)
    L = 1 + div(inst.deg, 2)
    P1L = Ps[1][:, 1:L]
    d = [P1L * rand(L) for _ in 1:m]
    # u = [ones(U) for _ in 1:n]
    u = [5 * ones(U) for _ in 1:n]

    # build model
    model = JuMP.Model(opt)
    # model = JuMP.Model()
    JuMP.@variable(model, x[1:n], Bin)
    JuMP.@variable(model, y[1:n, 1:m, 1:U])
    JuMP.@constraint(model, [i in 1:n, j in 1:m], y[i, j, :] in K)

    JuMP.@objective(
        model,
        Min,
        sum(JuMP.dot(c[i, j] * w, y[i, j, :]) for j in 1:m, i in 1:n) + JuMP.dot(f, x)
    )

    JuMP.@constraint(model, [i in 1:n], u[i] * x[i] - sum(y[i, j, :] for j in 1:m) in K)
    JuMP.@constraint(model, [j in 1:m], sum(y[i, j, :] for i in 1:n) - d[j] in K)

    # TODO
    JuMP.optimize!(model)
    @show JuMP.objective_value(model)
    @show JuMP.value.(x)
    @show JuMP.value.(y)

    return model
end

# TODO
inst = PolyFacilityLocation(2, 2, 2)
model = build(inst)
