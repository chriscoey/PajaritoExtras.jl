#=
vector regression problem with options for loss, regularization, priors

sets:
i ∈ 1..n   observations
j ∈ 1..m   independent variables
l ∈ 1..L   prior constraints

parameters:
Xᵢⱼ ∈ ℝ   jth element of independent variables for observation i
Yᵢ  ∈ ℝ   dependent variable for observation i
f   loss function
g   regularization function
hₗ  lth prior function
Hₗ  lth prior affine transformation
aₗ  ∈ ℝ   lth prior constant

variables:
βⱼ ∈ ℝ   jth parameter

objective:
min ∑ᵢ f(Yᵢ - Xᵢ' β) + g(β)

constraints:
hₗ (Hₗ β) ≤ aₗ, ∀l   lth prior
=#

abstract type RegrFun end


struct VectorRegression <: ExampleInstance
    n::Int
    m::Int
    f::RegrFun
    g::RegrFun
    hs::Vector{<:RegrFun}
end

function build(inst::VectorRegression)



    return model
end

function test_extra(inst::VectorRegression, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return
    # TODO
    return
end
