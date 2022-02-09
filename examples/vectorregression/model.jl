#=
vector regression problem with options for loss, regularization, priors

sets:
i ∈ 1..n   observations
j ∈ 1..m   independent variables
l ∈ 1..L   prior constraints

parameters:
Xᵢⱼ ∈ ℝ   jth element of independent variables for observation i
Yᵢ  ∈ ℝ   dependent variable for observation i
λ   ∈ ℝ₊  regularization parameter
f   loss function (domain ℝ)
hₗ   lth prior function (domain ℝ₊)
aₗ  ∈ ℝ   lth prior constant

variables:
βⱼ ∈ ℝ₊   jth parameter (assumed nonnegative)

objective:
min f(Y - X β) + λ L₀(β)   minimize loss plus L₀ regularization

constraints:
hₗ(β) ≤ aₗ, ∀l   lth prior
=#

struct VectorRegression <: ExampleInstance
    n::Int
    m::Int
    f_l2::Bool # if true, f is L2 norm, else L1 norm
    hs::Vector{<:VecSpecExt}
end

function build(inst::VectorRegression)
    (n, m) = (inst.n, inst.m)
    hs = inst.hs

    # generate parameters
    βmax = 3
    β0 = βmax * rand(m) # random solution used to construct data
    X = randn(n, m)
    Y = X * β0 + 0.1 * randn(n)
    λ = 0.1 # ?
    a = [get_val(β0, h_l) for h_l in hs]

    # TODO maybe need something like specext functions but for all the norm cones
    K_f = if inst.f_l2
        Hypatia.EpiNormEuclCone{Float64}(1 + n)
    else
        Hypatia.EpiNormInfCone{Float64, Float64}(1 + n, true)
    end

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, β[1:m])
    JuMP.@variable(model, z[1:m], Bin)
    JuMP.@constraint(model, β .<= βmax * z)

    JuMP.@variable(model, f_epi)
    JuMP.@constraint(model, vcat(f_epi, Y - X * β) in K_f)

    JuMP.@objective(model, Min, f_epi + λ * sum(z))

    for (h_l, a_l) in zip(hs, a)
        add_homog_spectral(h_l, m, vcat(a_l, β), model)
    end

    # save for use in tests
    model.ext[:z] = z

    return model
end

function test_extra(inst::VectorRegression, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    # check integer feasibility
    tol = eps()^0.2
    z = JuMP.value.(model.ext[:z])
    @test z ≈ round.(Int, z) atol = tol rtol = tol
    @test all(-tol .<= z .<= 1 + tol)
    return
end
