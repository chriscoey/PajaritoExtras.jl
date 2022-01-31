#=
given a nonnegative matrix with missing entries, find nonnegative integer values for the
missing entries to minimize the spectral or nuclear norm, subject to linear constraints
=#

using SparseArrays

struct MatrixCompletion <: ExampleInstance
    use_EF::Bool # in symmetric case, construct EF for spectral/nuclear cone
    nuclear_obj::Bool # use nuclear norm objective, else spectral norm
    symmetric::Bool # use symmetric matrix, else rectangular matrix
    nrow::Int # number of rows in matrix
    ncol::Int # not used in symmetric case
end

function build(inst::MatrixCompletion)
    nrow = inst.nrow
    @assert nrow >= 2
    symmetric = inst.symmetric
    nuclear_obj = inst.nuclear_obj
    ncol = (symmetric ? nrow : inst.ncol)
    @assert ncol >= nrow

    # set some options
    sparsity = min(0.7, 1 - inv(sqrt(nrow * ncol))) # fraction of values that are known
    max_rank = div(nrow, 2) # maximum rank of random matrix
    @show sparsity

    # generate data
    pattern = sprand(Bool, nrow, ncol, sparsity)
    if symmetric
        len = Cones.svec_length(nrow)
        triu!(pattern)
        A0 = rand(0:3, nrow, max_rank)
        A0 = A0 * A0'
    else
        len = nrow * ncol
        A0 = rand(0:3, nrow, max_rank) * rand(0:3, max_rank, ncol)
    end
    (rows, cols, _) = findnz(pattern)
    Avals = [A0[r, c] for (r, c) in zip(rows, cols)]
    num_unknown = len - length(Avals)
    @assert num_unknown >= 1

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, t)
    JuMP.@objective(model, Min, t / nrow)

    JuMP.@variable(model, x[1:num_unknown], Int)
    JuMP.@constraint(model, x .>= 0)
    JuMP.@constraint(model, x .<= 9)
    JuMP.@constraint(model, sum(x) >= length(x))

    X = Matrix{JuMP.AffExpr}(undef, nrow, ncol)
    for (r, c, a) in zip(rows, cols, Avals)
        X[r, c] = a
    end
    k = 1
    for c in 1:ncol, r in 1:(symmetric ? c : nrow)
        if !isassigned(X, r, c)
            X[r, c] = x[k]
            k += 1
        end
    end
    if symmetric
        LinearAlgebra.copytri!(X, 'U')
    end
    @assert k == length(x) + 1
    X ./= nrow # rescale for numerics

    if symmetric
        LinearAlgebra.copytri!(X, 'U')
        if inst.use_EF
            if nuclear_obj
                # EF for symmetric nuclear norm (L1 of eigvals)
                JuMP.@variable(model, X1[1:nrow, 1:nrow], PSD)
                JuMP.@variable(model, X2[1:nrow, 1:nrow], PSD)
                eq = [X[i, j] - X1[i, j] + X2[i, j] for j in 1:nrow for i in 1:j]
                JuMP.@constraint(model, eq .== 0)
                JuMP.@constraint(model, t >= tr(X1) + tr(X2))
            else
                # EF for symmetric spectral norm (Linf of eigvals)
                tI = t * Matrix(I, nrow, nrow)
                JuMP.@constraint(model, Symmetric(tI - X) in JuMP.PSDCone())
                JuMP.@constraint(model, Symmetric(tI + X) in JuMP.PSDCone())
            end
        else
            K = Hypatia.EpiNormSpectralTriCone{Float64, Float64}(1 + len, nuclear_obj)
            Xvec = Vector{JuMP.AffExpr}(undef, len)
            Cones.smat_to_svec!(Xvec, X, sqrt(2.0))
            JuMP.@constraint(model, vcat(t, Xvec) in K)
        end
    else
        if inst.use_EF
            if nuclear_obj
                # EF for nuclear norm
                X1 = JuMP.@variable(model, [1:nrow, 1:nrow], Symmetric)
                X2 = JuMP.@variable(model, [1:ncol, 1:ncol], Symmetric)
                JuMP.@constraint(model, t >= (tr(X1) + tr(X2)) / 2)
            else
                # EF for spectral norm
                X1 = t * Matrix(I, nrow, nrow)
                X2 = t * Matrix(I, ncol, ncol)
            end
            mat = hvcat((2, 2), X1, X, X', X2)
            JuMP.@constraint(model, Symmetric(mat) in JuMP.PSDCone())
        else
            K = Hypatia.EpiNormSpectralCone{Float64, Float64}(nrow, ncol, nuclear_obj)
            JuMP.@constraint(model, vcat(t, vec(X)) in K)
        end
    end

    # save for use in tests
    model.ext[:x_var] = x
    model.ext[:X_var] = X

    return model
end

function test_extra(inst::MatrixCompletion, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    # check objective and feasibility
    tol = eps()^0.2
    x_opt = JuMP.value.(model.ext[:x_var])
    @show x_opt
    @test x_opt ≈ round.(Int, x_opt) atol = tol rtol = tol

    X_opt = JuMP.value.(model.ext[:X_var])
    s = (inst.symmetric ? abs.(eigvals(Symmetric(X_opt, :U))) : svdvals(X_opt))
    snorm = (inst.nuclear_obj ? sum(s) : maximum(s))
    @test JuMP.objective_value(model) * inst.nrow ≈ snorm atol = tol rtol = tol
    return
end
