#=
given a nonnegative matrix with missing entries, find nonnegative integer values for the
missing entries to minimize the spectral or nuclear norm, subject to linear constraints
=#

struct MatrixCompletion <: ExampleInstance
    nuclear_obj::Bool # use nuclear norm objective, else spectral norm
    symmetric::Bool # use symmetric matrix, else rectangular matrix
    nrow::Int # number of rows in matrix
    ncol::Int # not used in symmetric case
    use_nat::Bool # use spectral/nuclear norm cone, else PSD cone formulation
end

function build(inst::MatrixCompletion)
    nrow = inst.nrow
    @assert nrow >= 2
    ncol = (inst.symmetric ? nrow : inst.ncol)
    @assert ncol >= nrow
    # set some options
    sparsity = min(0.7, 1 - inv(sqrt(nrow * ncol))) # fraction of values that are known
    max_rank = div(nrow, 2) # maximum rank of random matrix

    # generate data
    pattern = sprand(Bool, nrow, ncol, sparsity)
    if inst.symmetric
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
    for c in 1:ncol, r in 1:(inst.symmetric ? c : nrow)
        if !isassigned(X, r, c)
            X[r, c] = x[k]
            k += 1
        end
    end
    if inst.symmetric
        LinearAlgebra.copytri!(X, 'U')
        add_specnuc = add_symm_specnuc
    else
        add_specnuc = add_nonsymm_specnuc
    end
    @assert k == length(x) + 1
    X ./= nrow # rescale for numerics

    add_specnuc(inst.use_nat, inst.nuclear_obj, t, X, model)

    # save for use in tests
    model.ext[:x] = x
    model.ext[:X] = X

    return model
end

function test_extra(inst::MatrixCompletion, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat in OPT_OR_LIMIT
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    x = JuMP.value.(model.ext[:x])
    @test x ≈ round.(Int, x) atol = tol rtol = tol
    # check objective
    X = JuMP.value.(model.ext[:X])
    s = (inst.symmetric ? abs.(eigvals(Symmetric(X, :U))) : svdvals(X))
    snorm = (inst.nuclear_obj ? sum(s) : maximum(s))
    @test JuMP.objective_value(model) * inst.nrow ≈ snorm atol = tol rtol = tol
    return
end
