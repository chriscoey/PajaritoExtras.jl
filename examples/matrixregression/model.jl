#=
regularized matrix regression with row-sparsity

min ‖Y - X A‖ + λₙ ‖A‖⋆ :
Aᵢⱼ ∈ [-1,1], ∀i,j
ℓ₀(‖A₁‖₁, …, ‖Aₙ‖₁) ≤ k
where:
- design matrix X is n x p
- response matrix Y is n x m
- coefficient matrix variable A is p x m
- Aᵢ is the ith row of A
=#

struct MatrixRegression <: ExampleInstance
    n::Int # number of samples
    p::Int # number of predictors
    m::Int # number of responses
    use_nat::Bool # use nuclear norm cone, else PSD cone formulation
end

function build(inst::MatrixRegression)
    (n, m, p) = (inst.n, inst.m, inst.p)
    @assert n >= p >= m > 1
    k = ceil(Int, p * 2 / 3)

    # generate data
    A0 = 2 * rand(p, m) .- 1
    A0[(k + 1):end, :] .= 0
    X = randn(n, p)
    Y = X * A0 + 0.1 * randn(n, m)

    # use dimension reduction via QR
    F = qr(X, ColumnNorm())
    Y = (F.Q' * Y)[1:p, :]
    X = F.R[1:p, 1:p] * F.P'
    λ = 0.2 # TODO?

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, loss)
    JuMP.@variable(model, nuc)
    JuMP.@variable(model, row[1:p])
    JuMP.@variable(model, A[1:p, 1:m])
    JuMP.@variable(model, row_on[1:p], Bin)

    JuMP.@objective(model, Min, loss + λ * nuc)

    JuMP.@constraints(model, begin
        vcat(loss, vec(Y - X * A)) in JuMP.SecondOrderCone()
        [i = 1:p], A[i, :] .<= row_on[i]
        [i = 1:p], A[i, :] .>= -row_on[i]
        sum(row_on) <= k
    end)

    add_nonsymm_specnuc(inst.use_nat, true, nuc, A', model)

    # save for use in tests
    model.ext[:λ] = λ
    model.ext[:X] = X
    model.ext[:Y] = Y
    model.ext[:A] = A
    model.ext[:row_on] = row_on

    return model
end

function test_extra(inst::MatrixRegression, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    # check objective
    tol = eps()^0.2
    A = JuMP.value.(model.ext[:A])
    (Y, X) = (model.ext[:Y], model.ext[:X])
    loss = norm(Y - X * A)
    nuc = sum(svdvals(A))
    obj_result = loss + model.ext[:λ] * nuc
    @test JuMP.objective_value(model) ≈ obj_result atol = tol rtol = tol

    # check feasibility
    row_on = JuMP.value.(model.ext[:row_on])
    @test row_on ≈ round.(Int, row_on) atol = tol rtol = tol
    @test all(-tol .<= row_on .<= 1 + tol)
    for i in 1:size(A, 1)
        if abs(row_on[i]) < tol
            @assert norm(A[i, :]) < tol
        end
    end
    return
end
