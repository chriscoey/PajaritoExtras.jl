#=
estimate an inverse covariance matrix given an empirical covariance matrix Σ, using a
convex spectral function f, a constraint on the L0 norm of the lower triangle to induce
sparsity, constraints on the sparsity structure, and optionally a convex regularizer g
=#

struct InverseCovariance <: ExampleInstance
    d::Int
    g_l2::Union{Nothing, Bool} # g is L2 norm (true), L1 norm (false), or zero
    f::MatSpecExt # formulation specifier
end

function build(inst::InverseCovariance)
    d = inst.d
    @assert d >= 2
    @assert is_domain_pos(inst.f)
    sparsity = 0.3 # constrain the overall number of nonzeros
    row_sparsity = 0.4 # for each predictor, constrain the number of nonzeros in the row
    vec_dim = Cones.svec_length(d)
    L0_max = max(ceil(vec_dim * sparsity), 1)
    L0_row_max = max(ceil(d * row_sparsity), 1)

    # generate data
    Σ0 = randn(d, d)
    Σ0 = Σ0 * Σ0'
    Σ = Hermitian(inv(tr(Σ0)) * Σ0)
    Σ_vec = zeros(vec_dim)
    Cones.smat_to_svec!(Σ_vec, Σ, rt2)
    λ_f = 1.0 # for MLE, 1.0 with logdet
    λ_g = 0.1 # TODO
    P_max = 10 # TODO depends on reg

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, P_vec[1:vec_dim])

    # sparsity constraints
    JuMP.@variable(model, z[1:(vec_dim - d)], Bin)
    P_tri = P_vec[[i != j for j in 1:d for i in 1:j]]
    for (P_k, z_k) in zip(P_tri, z)
        JuMP.@constraint(model, P_k <= P_max * z_k)
        JuMP.@constraint(model, P_k >= -P_max * z_k)
    end
    JuMP.@constraint(model, sum(z) <= L0_max)

    Z = Matrix{Any}(undef, d, d)
    k = 1
    for j in 1:d
        Z[j, j] = 0.0
        for i in 1:(j - 1)
            Z[i, j] = Z[j, i] = z[k]
            k += 1
        end
    end
    JuMP.@constraint(model, [i in 1:d], sum(Z[i, :]) <= L0_row_max)

    # objective
    JuMP.@variable(model, f_epi)
    add_homog_spectral(inst.f, d, vcat(1.0 * f_epi, P_vec), model)

    if isnothing(inst.g_l2)
        g_epi = 0.0
    else
        g_epi = JuMP.@variable(model)
        K_g = if inst.g_l2
            Hypatia.EpiNormEuclCone{Float64}(1 + vec_dim)
        else
            Hypatia.EpiNormInfCone{Float64, Float64}(1 + vec_dim, true)
        end
        P_tri = scale_svec(Float64, 1.0 * P_vec, inv(rt2))
        JuMP.@constraint(model, vcat(g_epi, P_tri) in K_g)
    end

    JuMP.@objective(model, Min, dot(Σ_vec, P_vec) + λ_f * f_epi + λ_g * g_epi)

    # save for use in tests
    model.ext[:z] = z
    model.ext[:P_vec] = P_vec
    model.ext[:f_epi] = f_epi
    model.ext[:g_epi] = g_epi

    return model
end

function test_extra(inst::InverseCovariance, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat in OPT_OR_LIMIT
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    z = JuMP.value.(model.ext[:z])
    @test z ≈ round.(Int, z) atol = tol rtol = tol
    @test all(-tol .<= z .<= 1 + tol)
    # check conic feasibility
    P_vec = JuMP.value.(model.ext[:P_vec])
    P = Cones.svec_to_smat!(zeros(inst.d, inst.d), P_vec, rt2)
    λ = eigvals(Hermitian(P, :U))
    @test minimum(λ) >= -tol
    f_val = get_val(pos_only(λ), inst.f)
    f_epi = JuMP.value(model.ext[:f_epi])
    @test f_val ≈ f_epi atol = tol rtol = tol
    if !isnothing(inst.g_l2)
        P_tri = scale_svec(Float64, 1.0 * P_vec, inv(rt2))
        g_val = norm(P_tri, (inst.g_l2 ? 2 : 1))
        g_epi = JuMP.value(model.ext[:g_epi])
        @test g_val ≈ g_epi atol = tol rtol = tol
    end
    return
end
