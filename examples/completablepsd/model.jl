#=
given a symmetric sparsity pattern with some sparse entries unknown, find integer values
(within some bounds) for the unknowns such that the matrix is PSD-completable,
maximimizing the minimum eigenvalue

max y :
S + X - y I ∈ SparsePSD⋆
Xᵢⱼ ∈ ℤ, -M ≤ Xᵢⱼ ≤ M, ∀ (i,j) unknown

where:
- sparse parameter matrix S contains is the known entries
- sparse variable matrix X contains the unknown entries
- M is a bound on the absolute value of all unknown entries
=#

struct CompletablePSD <: ExampleInstance
    d::Int # side dimension
    sparsity::Real # fraction of elements in the pattern
    frac_known::Real # fraction of sparse entries known
    use_nat::Bool # use sparse PSD cone, else PSD cone formulation
end

function build(inst::CompletablePSD)
    d = inst.d
    @assert d >= 2

    # generate underlying positive definite matrix S0
    S0 = sprand(Bool, d, d, inst.sparsity)
    S0 = Matrix{Float64}(S0 * S0' + I)
    pattern = tril!(sprand(Bool, d, d, inst.sparsity) + I)
    (rows, cols, _) = findnz(pattern)
    S0vals = [S0[r, c] for (r, c) in zip(rows, cols)]
    num_sparse = length(S0vals)

    # select which entries of S0 will be known/unknown
    known_pattern = sprand(Bool, num_sparse, inst.frac_known)
    known_pattern[1] = known_pattern[end] = true
    num_unknown = num_sparse - sum(known_pattern)
    @assert num_unknown >= 1
    M = maximum(abs, S0vals[.!known_pattern])

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, y)
    JuMP.@objective(model, Max, y)

    JuMP.@variable(model, x[1:num_unknown], Int)
    JuMP.@constraint(model, -M .<= x .<= M)

    aff = zeros(JuMP.AffExpr, num_sparse)
    k = 1
    scal = (inst.use_nat ? rt2 : 1.0)
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
            aff[i] = scal * aff[i]
        end
    end
    @assert k == num_unknown + 1

    if inst.use_nat
        KT = Hypatia.PosSemidefTriSparseCone{Hypatia.Cones.PSDSparseDense, Float64, Float64}
        K = KT(d, rows, cols, true)
        JuMP.@constraint(model, aff in K)
    else
        mat = Matrix(sparse(rows, cols, aff))
        for i in 1:d, j in 1:i
            if iszero(mat[i, j])
                mat[i, j] = JuMP.@variable(model)
            end
        end
        LinearAlgebra.copytri!(mat, 'L')
        JuMP.@constraint(model, Symmetric(mat) in JuMP.PSDCone())
    end

    # save for use in tests
    model.ext[:x] = x
    model.ext[:S0_eigmin] = eigmin(Hermitian(S0))

    return model
end

function test_extra(inst::CompletablePSD, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    x = JuMP.value.(model.ext[:x])
    @test x ≈ round.(Int, x) atol = tol rtol = tol
    # check eigmin(S + X) >= eigmin(S0)
    S0_eigmin = model.ext[:S0_eigmin]
    @test JuMP.objective_value(model) >= S0_eigmin - tol
    return
end
