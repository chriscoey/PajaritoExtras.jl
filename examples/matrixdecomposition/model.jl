#=
decompose a matrix C into a sum of a sparse binary matrix A and a matrix B with small
spectral norm, constraining the number of nonzeros in A

given matrix C, solve:
min specnorm(C - A) :
Aᵢⱼ ∈ {0,1}
sum(A) = k

TODO add continuous variant to Hypatia: relax sparsity to L1
see http://www.mit.edu/~parrilo/pubs/talkfiles/ISMP2009.pdf
=#

struct MatrixDecomposition <: ExampleInstance
    m::Int
    n::Int
    use_nat::Bool # use nuclear norm cone, else PSD cone formulation
end

function build(inst::MatrixDecomposition)
    (m, n) = (inst.m, inst.n)
    B_rank = 1 + div(m, 2)
    @assert 1 <= B_rank <= m <= n

    # generate data
    A_sparsity = min(1.0, 2 / n)
    A0 = sprand(Bool, m, n, A_sparsity)
    k = nnz(A0)
    B0 = randn(m, B_rank) * randn(B_rank, n)
    C = A0 + B0

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, A[1:m, 1:n], Bin)
    JuMP.@constraint(model, sum(A) == k)

    JuMP.@variable(model, u)
    JuMP.@objective(model, Min, u)

    add_nonsymm_specnuc(inst.use_nat, false, u, C - A, model)

    # save for use in tests
    model.ext[:u] = u
    model.ext[:A] = A
    model.ext[:B] = C - A

    return model
end

function test_extra(inst::MatrixDecomposition, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat in OPT_OR_LIMIT
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    A = JuMP.value.(model.ext[:A])
    @test A ≈ round.(Int, A) atol = tol rtol = tol
    @test all(-tol .<= A .<= 1 + tol)
    # check nuclear norm value
    tol = eps()^0.2
    u = JuMP.value.(model.ext[:u])
    B = JuMP.value.(model.ext[:B])
    @test u ≈ maximum(svdvals(B)) atol = tol rtol = tol
    @test JuMP.objective_value(model) ≈ u
    return
end
