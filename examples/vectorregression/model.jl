#=
vector regression problem with options for loss, regularization, and prior information
expressed through convex constraints

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
hₗ(β) ≤ aₗ, ∀l   lth prior constraint
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
    λ = 0.1
    a = [get_val(β0, h_l) for h_l in hs]

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
    f_aff = JuMP.@expression(model, vcat(f_epi, Y - X * β))
    JuMP.@constraint(model, f_aff in K_f)

    JuMP.@objective(model, Min, f_epi + λ * sum(z))

    for (h_l, a_l) in zip(hs, a)
        add_homog_spectral(h_l, m, vcat(a_l, β), model)
    end

    # save for use in tests
    model.ext[:z] = z
    model.ext[:β] = β
    model.ext[:f_aff] = f_aff
    model.ext[:a] = a

    return model
end

function test_extra(inst::VectorRegression, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    z = JuMP.value.(model.ext[:z])
    @test z ≈ round.(Int, z) atol = tol rtol = tol
    @test all(-tol .<= z .<= 1 + tol)
    # check conic feasibility
    f_aff = JuMP.value.(model.ext[:f_aff])
    f_val = norm(f_aff[2:end], (inst.f_l2 ? 2 : 1))
    @test f_aff[1] ≈ f_val atol = tol rtol = tol
    β = JuMP.value.(model.ext[:β])
    @test all(a_l >= get_val(β, h_l) - tol for (h_l, a_l) in zip(inst.hs, model.ext[:a]))
    return
end
