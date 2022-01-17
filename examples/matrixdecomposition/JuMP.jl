#=
decompose a matrix into a sum of a sparse and a low rank matrix
for the rank condition, use the nuclear norm convex relaxation
assume all entries of the sparse matrix are in [0, 1]

given real or complex matrix C, solve:
min γ ‖A‖₀ + ‖B‖∗
s.t. A + B = C

TODO add continuous variant to Hypatia: relax sparsity to L1
see http://www.mit.edu/~parrilo/pubs/talkfiles/ISMP2009.pdf
=#

import SparseArrays

struct MatrixDecomposition <: ExampleInstance
    d1::Int
    d2::Int
    A_sparsity::Float64
    B_rank::Int
    # TODO complex version
end

function build(inst::MatrixDecomposition)
    (d1, d2) = (inst.d1, inst.d2)
    B_rank = inst.B_rank
    @assert 1 <= B_rank <= d1 <= d2
    A0 = SparseArrays.sprand(d1, d2, inst.A_sparsity)
    B0 = randn(d1, B_rank) * randn(B_rank, d2)
    C = A0 + B0
    γ = 0.02 # TODO?

    # build model
    model = JuMP.Model()

    JuMP.@variable(model, X[1:d1, 1:d2], Bin)
    JuMP.@variable(model, A[1:d1, 1:d2] >= 0)
    JuMP.@constraint(model, A .<= X)

    JuMP.@variable(model, B[1:d1, 1:d2])
    JuMP.@variable(model, u)
    K = Hypatia.EpiNormSpectralCone{Float64, Float64}(d1, d2, true)
    JuMP.@constraint(model, vcat(u, vec(B)) in K)

    JuMP.@constraint(model, A + B .== C)
    JuMP.@objective(model, Min, γ * sum(X) + u)

    # # TODO
    # JuMP.optimize!(model)
    # @show JuMP.objective_value(model)
    # A = JuMP.value.(A)
    # @show A
    # @show A0 - A
    # B = JuMP.value.(B)
    # @show LinearAlgebra.svdvals(B)
    # @show B
    # @show B0 - B
    # @show JuMP.value(u)

    return model
end

function test_extra(inst::MatrixDecomposition, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return
    # TODO
    return
end
