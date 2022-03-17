#=
estimate an inverse covariance matrix given an empirical covariance matrix Σ, using a
convex spectral function f and squared Frobenius norm regularization and a constraint on
the number of nonzero (off diagonal) rows
=#

struct InverseCovariance <: ExampleInstance
    d::Int
    f::MatSpecExt # formulation specifier
end

function build(inst::InverseCovariance)
    d = inst.d
    @assert d >= 3
    @assert is_domain_pos(inst.f)
    vec_dim = Cones.svec_length(d)
    nzrows = d - floor(Int, sqrt(d))

    # generate data
    Σ0 = randn(d, d)
    Σ0 = Σ0 * Σ0'
    Σ = Hermitian(inv(tr(Σ0)) * Σ0)
    Σ_vec = zeros(vec_dim)
    Cones.smat_to_svec!(Σ_vec, Σ, rt2)
    λ_f = 1.0 # for MLE, 1.0 with logdet
    λ_g = 0.1 # TODO

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, P_vec[1:vec_dim])

    JuMP.@variable(model, y[1:d], Bin)
    JuMP.@constraint(model, sum(y) == nzrows)
    P = Symmetric(get_smat(d, 1.0 * P_vec), :U)
    JuMP.@variable(model, z[1:(1 + d)])

    for i in 1:d
        aff = vcat(0.5 * z[i], y[i], P[1:(i - 1), i], P[(i + 1):end, i])
        JuMP.@constraint(model, aff in Hypatia.EpiPerSquareCone{Float64}(1 + d))
    end
    aff = vcat(0.5 * z[1 + d], 1, diag(P))
    JuMP.@constraint(model, aff in Hypatia.EpiPerSquareCone{Float64}(2 + d))

    JuMP.@variable(model, f_epi)
    add_homog_spectral(inst.f, d, vcat(1.0 * f_epi, P_vec), model)

    JuMP.@objective(model, Min, dot(Σ_vec, P_vec) + λ_f * f_epi + λ_g * sum(z))

    # save for use in tests
    model.ext[:y] = y
    model.ext[:z] = z
    model.ext[:P] = P
    model.ext[:f_epi] = f_epi

    return model
end

function test_extra(inst::InverseCovariance, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat in OPT_OR_LIMIT
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    y = JuMP.value.(model.ext[:y])
    @test y ≈ round.(Int, y) atol = tol rtol = tol
    @test all(-tol .<= y .<= 1 + tol)
    # check conic feasibility
    P = JuMP.value.(model.ext[:P])
    λ = eigvals(Hermitian(P, :U))
    @test minimum(λ) >= -tol
    f_val = get_val(pos_only(λ), inst.f)
    f_epi = JuMP.value(model.ext[:f_epi])
    @test f_val ≈ f_epi atol = tol rtol = tol
    z = JuMP.value.(model.ext[:z])
    @test sum(z) ≈ sum(abs2, P) atol = tol rtol = tol
    @test all((y[i] > 0.5 || sum(abs, P[i, :]) - P[i, i] < tol) for i in eachindex(y))
    return
end
