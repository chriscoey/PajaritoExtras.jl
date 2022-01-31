#=
given a symmetric sparsity pattern with some sparse entries unknown, find integer values for
the unknowns such that the matrix is PSD-completable, maximimizing the minimum eigenvalue

maximize    y
subject to  S + X - y I ∈ SparsePSD⋆
            Xᵢⱼ ∈ ℤ, ∀ (i,j) unknown
=#

using SparseArrays

struct CompletablePSD <: ExampleInstance
    d::Int # side dimension
    sparsity::Real # fraction of elements in the pattern
    frac_known::Real # fraction of sparse entries known
end

function build(inst::CompletablePSD)
    d = inst.d
    @assert d >= 2

    # generate data
    S0 = inv(sqrt(d)) * rand(-3:3, d, d)
    S0 = S0 * S0'
    pattern = tril!(sprand(Bool, d, d, inst.sparsity) + I)
    (rows, cols, _) = findnz(pattern)
    S0vals = [S0[r, c] for (r, c) in zip(rows, cols)]
    num_sparse = length(S0vals)
    known_pattern = sprand(Bool, num_sparse, inst.frac_known)
    known_pattern[1] = known_pattern[end] = true
    num_unknown = num_sparse - sum(known_pattern)
    @assert num_unknown >= 1
    # @show num_sparse, num_unknown

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, y)
    JuMP.@objective(model, Max, y)
    JuMP.@variable(model, x[1:num_unknown], Int)

    rt2 = sqrt(2.0)
    aff = zeros(JuMP.AffExpr, num_sparse)
    k = 1
    for i in eachindex(known_pattern)
        if known_pattern[i]
            aff[i] = S0vals[i]
        else
            aff[i] = x[k]
            k += 1
        end
        if rows[i] == cols[i]
            aff[i] = aff[i] - y
        else
            aff[i] = rt2 * aff[i]
        end
    end
    @assert k == num_unknown + 1

    KT = Hypatia.PosSemidefTriSparseCone{Hypatia.Cones.PSDSparseDense, Float64, Float64}
    K = KT(d, rows, cols, true)
    JuMP.@constraint(model, aff in K)

    # save for use in tests
    model.ext[:x] = x
    model.ext[:S_eigmin] = eigmin(Hermitian(S0))

    return model
end

function test_extra(inst::CompletablePSD, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    x_opt = JuMP.value.(model.ext[:x])
    @test x_opt ≈ round.(Int, x_opt) atol = tol rtol = tol

    S_eigmin = model.ext[:S_eigmin]
    @test JuMP.objective_value(model) >= S_eigmin - tol
    return
end