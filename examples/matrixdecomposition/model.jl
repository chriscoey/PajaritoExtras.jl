#=
decompose a matrix into a sum of a sparse and a low rank matrix
for the rank condition, use the nuclear norm convex relaxation
assume all entries of the sparse matrix are in [0, 1]

given real or complex matrix C, solve:
min λ ‖A‖₀ + ‖B‖∗
s.t. A + B = C

TODO add continuous variant to Hypatia: relax sparsity to L1
see http://www.mit.edu/~parrilo/pubs/talkfiles/ISMP2009.pdf
=#

struct MatrixDecomposition <: ExampleInstance
    d1::Int
    d2::Int
    A_sparsity::Float64
    B_rank::Int
    use_nat::Bool # use nuclear norm cone, else PSD cone formulation
end

function build(inst::MatrixDecomposition)
    (d1, d2) = (inst.d1, inst.d2)
    B_rank = inst.B_rank
    @assert 1 <= B_rank <= d1 <= d2

    # generate data
    A0 = sprand(d1, d2, inst.A_sparsity)
    B0 = randn(d1, B_rank) * randn(B_rank, d2)
    C = A0 + B0
    λ = 0.02 # TODO?

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, X[1:d1, 1:d2], Bin)
    JuMP.@variable(model, A[1:d1, 1:d2] >= 0)
    JuMP.@constraint(model, A .<= X)

    JuMP.@variable(model, B[1:d1, 1:d2])
    JuMP.@variable(model, u)
    add_nonsymm_specnuc(inst.use_nat, true, u, B, model)

    JuMP.@constraint(model, A + B .== C)
    JuMP.@objective(model, Min, λ * sum(X) + u)

    # save for use in tests
    model.ext[:u] = u
    model.ext[:B] = B

    return model
end

function test_extra(inst::MatrixDecomposition, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    # check nuclear norm value
    tol = eps()^0.2
    u = JuMP.value.(model.ext[:u])
    B = JuMP.value.(model.ext[:B])
    nuc = sum(svdvals(B))
    @test u ≈ nuc atol = tol rtol = tol
    return
end
