
function possemideftri1(opt)
    TOL = 1e-4
    m = JuMP.Model(opt)

    JuMP.@variable(m, x, Int)
    JuMP.@constraint(m, x >= 0)
    JuMP.@variable(m, y >= 0)
    JuMP.@variable(m, Z[1:2, 1:2], Symmetric)
    JuMP.@constraint(m, svec(1.0 * Z) in Hypatia.PosSemidefTriCone{Float64, Float64}(3))
    JuMP.@objective(m, Max, 3x + y - Z[1, 1])
    JuMP.@constraint(m, 3x + 2y <= 10)
    JuMP.@constraint(m, svec([2 x; x 2]) in Hypatia.PosSemidefTriCone{Float64, Float64}(3))
    c1 = JuMP.@constraint(m, Z[1, 2] >= 1)
    c2 = JuMP.@constraint(m, y >= Z[2, 2])
    JuMP.optimize!(m)

    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), 7.5, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), 7.5, atol = TOL)
    @test isapprox(JuMP.value(x), 2, atol = TOL)
    @test isapprox(JuMP.value(y), 2, atol = TOL)
    @test isapprox(JuMP.value.(vec(Z)), [0.5, 1, 1, 2], atol = TOL)

    JuMP.delete(m, c1)
    JuMP.delete(m, c2)
    JuMP.@constraint(m, x >= 2)
    JuMP.set_lower_bound(Z[1, 2], 2)
    c2 = JuMP.@constraint(m, y >= Z[2, 2] + Z[1, 1])
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.INFEASIBLE
    @test JuMP.primal_status(m) == MOI.NO_SOLUTION
    return
end

function possemideftri2(opt)
    TOL = 1e-4
    d = 3
    for is_complex in (false, true)
        R = (is_complex ? ComplexF64 : Float64)
        dim = svec_length(R, d)
        mat = Hermitian(smat(R, Float64.(1:dim)), :U)
        λ₁ = eigmax(mat)
        m = JuMP.Model(opt)

        x = JuMP.@variable(m, [1:dim])
        JuMP.@constraint(m, x in Hypatia.PosSemidefTriCone{Float64, R}(dim))
        JuMP.@objective(m, Max, dot(svec(mat), x))
        x_diag = [x[svec_idx(R, i, i)] for i in 1:d]
        JuMP.@constraint(m, sum(x_diag) == 1)
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.objective_value(m), λ₁, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), λ₁, atol = TOL)
        x_val = JuMP.value.(x)
        x_mat = Hermitian(smat(R, x_val), :U)
        @test isapprox(tr(x_mat), 1, atol = TOL)
        @test eigmin(x_mat) >= -TOL

        JuMP.set_binary.(x_diag)
        JuMP.optimize!(m)
        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        x_val = JuMP.value.(x)
        x_mat = Hermitian(smat(R, x_val), :U)
        @test isapprox(tr(x_mat), 1, atol = TOL)
        @test eigmin(x_mat) >= -TOL
        @test isapprox(x_mat[d, d], 1, atol = TOL)
    end
    return
end

function possemideftrisparse1(opt)
    # arrow PSD formulation of second-order cone constraint
    TOL = 1e-4
    row_idxs = [1, 2, 3, 4, 2, 3, 4]
    col_idxs = [1, 1, 1, 1, 2, 3, 4]

    for is_complex in (false, true)
        R = (is_complex ? ComplexF64 : Float64)
        KT = Hypatia.PosSemidefTriSparseCone{Hypatia.Cones.PSDSparseDense, Float64, R}
        K = KT(4, row_idxs, col_idxs, false)
        x_dim = vec_length(R, 3)
        x0 = vec_copyto!(zeros(x_dim), -Real.(1:x_dim) .+ 0.5)
        opt_val = norm(Real.(1:x_dim), 2)
        m = JuMP.Model(opt)

        JuMP.@variable(m, y)
        JuMP.@objective(m, Min, y)
        JuMP.@variable(m, x[1:x_dim], Int)
        JuMP.@constraint(m, x .<= x0)
        aff = vcat(y, rt2 * x, y, y, y)
        JuMP.@constraint(m, aff in K)
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        @test isapprox(JuMP.value.(x), -Float64.(1:x_dim), atol = TOL)

        JuMP.@constraint(m, y <= opt_val - 0.2)
        JuMP.optimize!(m)
        @test JuMP.termination_status(m) == MOI.INFEASIBLE
        @test JuMP.primal_status(m) == MOI.NO_SOLUTION
    end
    return
end

function possemideftrisparse2(opt)
    TOL = 1e-4
    row_idxs = [1, 2, 3, 2, 3, 3, 4, 4, 5, 5]
    col_idxs = [1, 1, 1, 2, 2, 3, 3, 4, 4, 5]
    KT = Hypatia.PosSemidefTriSparseCone{Hypatia.Cones.PSDSparseDense, Float64, Float64}
    K = KT(5, row_idxs, col_idxs, true)
    m = JuMP.Model(opt)

    JuMP.@variable(m, y[1:3], Int)
    JuMP.@objective(m, Min, sum(y))
    aff = [8, -2, -4, y[1], 7, y[2], 5, y[3], -6, 5]
    JuMP.@constraint(m, y .>= 3)
    JuMP.@constraint(m, aff in K)

    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), 15, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), 15, atol = TOL)

    JuMP.@constraint(m, y[3] == 3)
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.INFEASIBLE
    @test JuMP.primal_status(m) == MOI.NO_SOLUTION
    return
end

function epinormeucl1(opt)
    TOL = 1e-4
    m = JuMP.Model(opt)

    JuMP.@variable(m, x)
    JuMP.@variable(m, y)
    JuMP.@variable(m, z, Int)
    JuMP.@constraint(m, z <= 2.5)
    JuMP.@objective(m, Min, x + 2y)
    JuMP.@constraint(m, [z, x, y] in Hypatia.EpiNormEuclCone{Float64}(3))
    JuMP.set_integer(x)
    JuMP.optimize!(m)

    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    opt_obj = -1 - 2 * sqrt(3)
    @test isapprox(JuMP.objective_value(m), opt_obj, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), opt_obj, atol = TOL)
    @test isapprox(JuMP.value(x), -1, atol = TOL)
    @test isapprox(JuMP.value(y), -sqrt(3), atol = TOL)
    @test isapprox(JuMP.value(z), 2, atol = TOL)

    JuMP.unset_integer(x)
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    opt_obj = -2 * sqrt(5)
    @test isapprox(JuMP.objective_value(m), opt_obj, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), opt_obj, atol = TOL)
    @test isapprox(abs2(JuMP.value(x)) + abs2(JuMP.value(y)), 4, atol = TOL)
    @test isapprox(JuMP.value(z), 2, atol = TOL)
    return
end

function epinormeucl2(opt)
    TOL = 1e-4
    m = JuMP.Model(opt)

    JuMP.@variable(m, x[1:3], Bin)
    K = Hypatia.EpiNormEuclCone{Float64}(4)
    c1 = JuMP.@constraint(m, vcat(inv(sqrt(2)), x .- 0.5) in K)
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.INFEASIBLE
    @test JuMP.primal_status(m) == MOI.NO_SOLUTION

    JuMP.@objective(m, Min, 1 + sum(x))
    JuMP.delete(m, c1)
    JuMP.@constraint(m, vcat(1.5, x .- 1) in Hypatia.EpiNormEuclCone{Float64}(4))
    JuMP.set_start_value.(x, [1, 0, 0])
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), 2, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), 2, atol = TOL)
    return
end

function epipersquare1(opt)
    m = JuMP.Model(opt)
    JuMP.@variable(m, x[1:3], Int)
    K = Hypatia.EpiPerSquareCone{Float64}(5)
    JuMP.@constraint(m, vcat(0.25, 1, x .- 0.5) in K)
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.INFEASIBLE
    @test JuMP.primal_status(m) == MOI.NO_SOLUTION
    return
end

function epipersquare2(opt)
    TOL = 1e-4
    m = JuMP.Model(opt)

    JuMP.@variable(m, v, Bin)
    JuMP.@variable(m, w[1:2], Int)
    JuMP.@constraint(m, vcat(1, v, w .- 1.5) in Hypatia.EpiPerSquareCone{Float64}(4))
    JuMP.@objective(m, Min, v)
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), 1, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), 1, atol = TOL)
    @test isapprox(JuMP.value(v), 1, atol = TOL)

    JuMP.unset_integer.(w)
    JuMP.@objective(m, Min, v + sum(w) - 1)
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), 1, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), 1, atol = TOL)
    @test isapprox(JuMP.value.(w), [0.5, 0.5], atol = 10TOL)
    return
end

function epinorminf1(opt)
    TOL = 1e-4
    for use_dual in (false, true), is_complex in (false, true)
        R = (is_complex ? ComplexF64 : Float64)
        u_fix = (is_complex ? rt2 : 1.0) * (use_dual ? 0.5 : 0.4)
        w_fix = (is_complex ? [0.3, 0.3, 0.7, 0.7, 1.2, 1.2] : [0.3, 0.7, 1.2])
        w_opt = round.(w_fix)
        dim = vec_length(R, 3)
        m = JuMP.Model(opt)

        JuMP.@variable(m, w[1:dim], Bin)
        K = Hypatia.EpiNormInfCone{Float64, R}(1 + dim, use_dual)
        c1 = JuMP.@constraint(m, vcat(u_fix, w .- 0.5) in K)
        JuMP.optimize!(m)
        @test JuMP.termination_status(m) == MOI.INFEASIBLE
        @test JuMP.primal_status(m) == MOI.NO_SOLUTION

        JuMP.@objective(m, Min, 1 + sum(w))
        JuMP.delete(m, c1)
        K = Hypatia.EpiNormInfCone{Float64, R}(1 + dim, use_dual)
        if use_dual
            u_fix *= 2
        end
        JuMP.@constraint(m, vcat(u_fix, w - w_fix) in K)
        JuMP.set_start_value.(w, w_opt)
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.objective_value(m), 1 + sum(w_opt), atol = TOL)
        @test isapprox(JuMP.objective_bound(m), 1 + sum(w_opt), atol = TOL)
        @test isapprox(JuMP.value.(w), w_opt, atol = TOL)
    end
    return
end

function epinorminf2(opt)
    TOL = 1e-4
    for use_dual in (false, true), is_complex in (false, true)
        R = (is_complex ? ComplexF64 : Float64)
        dim = vec_length(R, 2)
        m = JuMP.Model(opt)

        JuMP.@variable(m, w[1:dim], Int)
        K = Hypatia.EpiNormInfCone{Float64, R}(1 + dim, use_dual)
        aff = [i * w[i] for i in eachindex(w)]
        JuMP.@constraint(m, vcat(10.5, aff) in K)
        JuMP.@objective(m, Min, sum(w))
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        w_opt = JuMP.value.(w)
        @test isapprox(JuMP.objective_value(m), sum(w_opt), atol = TOL)
        @test isapprox(JuMP.objective_bound(m), sum(w_opt), atol = TOL)
        aff_opt = vec_copyto!(zeros(R, 2), JuMP.value.(aff))
        @test norm(aff_opt, (use_dual ? 1 : Inf)) <= 10.5 + TOL
    end
    return
end

function epinormspectral1(opt)
    TOL = 1e-4
    d1 = 2
    d2 = 3
    for use_dual in (false, true), is_complex in (false, true)
        R = (is_complex ? ComplexF64 : Float64)
        dim = vec_length(R, d1 * d2)
        vec = collect(1:dim)
        mat = vec_copyto!(zeros(R, d1, d2), fill(0.1, dim))
        σ = svdvals(mat)
        opt_val = (use_dual ? sum(σ) : σ[1])
        m = JuMP.Model(opt)

        JuMP.@variable(m, y)
        x = JuMP.@variable(m, [1:dim], Int)
        JuMP.@constraint(m, [i in 1:dim], i - 1 <= x[i] <= i + 1)
        K = Hypatia.EpiNormSpectralCone{Float64, R}(d1, d2, use_dual)
        JuMP.@constraint(m, vcat(y, vec - x .+ 0.1) in K)
        JuMP.@objective(m, Min, y)
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        @test isapprox(JuMP.value.(x), vec, atol = TOL)
    end
    return
end

function epinormspectral2(opt)
    TOL = 1e-4
    d1 = d2 = 2
    for use_dual in (false, true), is_complex in (false, true)
        R = (is_complex ? ComplexF64 : Float64)
        dim = vec_length(R, d1 * d2)
        opt_val = (use_dual ? 1 : d1)
        m = JuMP.Model(opt)

        x = JuMP.@variable(m, [1:dim], Bin)
        K = Hypatia.EpiNormSpectralCone{Float64, R}(d1, d2, use_dual)
        JuMP.@constraint(m, vcat(1, x) in K)
        JuMP.@objective(m, Max, sum(x))
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        @test isapprox(sum(JuMP.value.(x)), opt_val, atol = TOL)
    end
    return
end

function hypogeomean1(opt)
    TOL = 1e-4
    m = JuMP.Model(opt)

    JuMP.@variable(m, u)
    JuMP.@variable(m, w[1:3], Int)
    JuMP.@constraint(m, vcat(u, w) in Hypatia.HypoGeoMeanCone{Float64}(4))
    JuMP.@objective(m, Max, u)
    w_max = [1.1, 2.3, 3.5]
    JuMP.@constraint(m, w .<= w_max)
    JuMP.optimize!(m)

    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    opt_obj = PajaritoExtras.geomean(Float64.(1:3))
    @test isapprox(JuMP.objective_value(m), opt_obj, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), opt_obj, atol = TOL)
    @test isapprox(JuMP.value(u), opt_obj, atol = TOL)
    @test isapprox(JuMP.value.(w), 1:3, atol = TOL)

    JuMP.unset_integer.(w)
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    opt_obj = PajaritoExtras.geomean(w_max)
    @test isapprox(JuMP.objective_value(m), opt_obj, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), opt_obj, atol = TOL)
    @test isapprox(JuMP.value.(w), w_max, atol = TOL)
    return
end

function hypogeomean2(opt)
    TOL = 1e-4
    m = JuMP.Model(opt)
    JuMP.@variable(m, x, Bin)
    JuMP.@variable(m, y[1:3], Int)
    JuMP.@constraint(m, vcat(0.5 * x, y) in Hypatia.HypoGeoMeanCone{Float64}(4))
    JuMP.@objective(m, Max, 5x - sum(y))
    JuMP.optimize!(m)

    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), 2, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), 2, atol = TOL)
    @test isapprox(JuMP.value(x), 1, atol = TOL)
    @test isapprox(JuMP.value.(y), ones(3), atol = TOL)
    return
end

function hyporootdettri1(opt)
    TOL = 1e-4
    for d in (2, 3), is_complex in (false, true)
        R = (is_complex ? ComplexF64 : Float64)
        w_dim = svec_length(R, d)
        m = JuMP.Model(opt)

        w = JuMP.@variable(m, [1:w_dim], Bin)
        K = Hypatia.HypoRootdetTriCone{Float64, R}(1 + w_dim)
        JuMP.@constraint(m, vcat(1, w) in K)
        JuMP.@objective(m, Max, sum(w))
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.objective_value(m), d, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), d, atol = TOL)
        w_val = JuMP.value.(w)
        w_mat = Hermitian(smat(R, w_val), :U)
        @test isapprox(w_mat, I, atol = TOL)

        JuMP.@constraint(m, sum(w) <= d - 1)
        JuMP.optimize!(m)
        @test JuMP.termination_status(m) == MOI.INFEASIBLE
        @test JuMP.primal_status(m) == MOI.NO_SOLUTION
    end
    return
end

function hyporootdettri2(opt)
    TOL = 1e-4
    (m, x, Q, opt_x, opt_Q) = _setup_expdesign(opt)
    opt_val = sqrt(det(opt_Q))
    K = Hypatia.HypoRootdetTriCone{Float64, Float64}(4)

    JuMP.@variable(m, y)
    JuMP.@objective(m, Max, y)
    JuMP.@constraint(m, vcat(y, svec(Q)) in K)
    JuMP.optimize!(m)

    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
    @test isapprox(JuMP.value.(x), opt_x, atol = TOL)
    return
end

function vector_epipersepspectral1(opt)
    TOL = 1e-4
    sep_spectral_funs = [
        Hypatia.Cones.NegLogSSF(),
        Hypatia.Cones.NegEntropySSF(),
        Hypatia.Cones.NegSqrtSSF(),
        Hypatia.Cones.NegPower01SSF(0.3),
        Hypatia.Cones.Power12SSF(3 // 2),
    ]

    for h_fun in sep_spectral_funs, use_dual in (false, true)
        D = (use_dual ? Dual : Prim)
        Q = Hypatia.Cones.VectorCSqr{Float64}
        K = Hypatia.EpiPerSepSpectralCone{Float64}(h_fun, Q, 3, use_dual)
        m = JuMP.Model(opt)

        JuMP.@variable(m, u)
        JuMP.@variable(m, w[1:3], Int)
        JuMP.@constraint(m, w .<= [3.2, 5.5, 6.9])
        JuMP.@constraint(m, sum(w) == 8)
        JuMP.@objective(m, Min, u)
        epiper = PajaritoExtras.swap_epiper(D, [u, 1]...)
        JuMP.@constraint(m, vcat(epiper..., w) in K)
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        w_sol = JuMP.value.(w)
        opt_val = PajaritoExtras.val_or_conj(D)(w_sol, h_fun)
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        @test isapprox(JuMP.value(u), opt_val, atol = TOL)
        @test isapprox(w_sol, round.(w_sol), atol = TOL)
    end
    return
end

function vector_epipersepspectral2(opt)
    TOL = 1e-4
    # only functions for which the perspective is zero when the perspective variable is zero
    sep_spectral_funs = [
        Hypatia.Cones.NegLogSSF(),
        Hypatia.Cones.NegSqrtSSF(),
        Hypatia.Cones.NegPower01SSF(3 // 10),
    ]

    for h_fun in sep_spectral_funs, use_dual in (false, true)
        D = (use_dual ? Dual : Prim)
        Q = Hypatia.Cones.VectorCSqr{Float64}
        w = [0.5, 2.3]
        m = JuMP.Model(opt)

        JuMP.@variable(m, x[1:2], Bin)
        JuMP.@constraint(m, sum(x) == 1)
        JuMP.@variable(m, u[1:2])
        JuMP.@objective(m, Min, sum(u))
        for i in 1:2
            K = Hypatia.EpiPerSepSpectralCone{Float64}(h_fun, Q, 2, use_dual)
            epiper = PajaritoExtras.swap_epiper(D, u[i], x[i])
            JuMP.@constraint(m, vcat(epiper..., w .+ i) in K)
        end
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        x_sol = JuMP.value.(x)
        round_x_sol = round.(x_sol)
        @test isapprox(x_sol, round_x_sol, atol = TOL)
        @test round_x_sol in ([0.0, 1.0], [1.0, 0.0])
        i = findfirst(Bool.(round_x_sol))
        u_sol = JuMP.value.(u)
        opt_val = PajaritoExtras.val_or_conj(D)(w .+ i, h_fun)
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        @test isapprox(u_sol[i], opt_val, atol = TOL)
        @test isapprox(sum(u_sol), opt_val, atol = TOL)
    end
    return
end

function matrix_epipersepspectral1(opt)
    TOL = 1e-4
    sep_spectral_funs = [
        Hypatia.Cones.NegLogSSF(),
        Hypatia.Cones.NegEntropySSF(),
        Hypatia.Cones.NegSqrtSSF(),
        Hypatia.Cones.NegPower01SSF(0.8),
        Hypatia.Cones.Power12SSF(1.2),
    ]

    for h_fun in sep_spectral_funs, use_dual in (false, true), is_complex in (false, true)
        R = (is_complex ? ComplexF64 : Float64)
        D = (use_dual ? Dual : Prim)
        Q = Hypatia.Cones.MatrixCSqr{Float64, R}
        K = Hypatia.EpiPerSepSpectralCone{Float64}(h_fun, Q, 3, use_dual)
        w_dim = svec_length(R, 3)
        vec_diag = svec(Matrix{R}(Diagonal([1.2, 2.5, 0.9])))
        m = JuMP.Model(opt)

        JuMP.@variable(m, u)
        JuMP.@variable(m, w[1:w_dim])
        W_diag = [w[svec_idx(R, i, i)] for i in 1:3]
        JuMP.set_integer.(W_diag)
        JuMP.@constraint(m, sum(W_diag) == 8)
        JuMP.@constraint(m, sum(w) >= 9)
        JuMP.@objective(m, Min, u)
        K_vec = w - vec_diag
        epiper = PajaritoExtras.swap_epiper(D, [u, 1]...)
        JuMP.@constraint(m, vcat(epiper..., K_vec) in K)
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        diag_sol = JuMP.value.(W_diag)
        @test isapprox(diag_sol, round.(diag_sol), atol = TOL)
        ω = eigvals(Hermitian(smat(R, JuMP.value.(K_vec)), :U))
        if PajaritoExtras.dom_pos(D, h_fun)
            @test minimum(ω) > -TOL
            ω = max.(ω, 1e-9)
        end
        opt_val = PajaritoExtras.val_or_conj(D)(ω, h_fun)
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        @test isapprox(JuMP.value(u), opt_val, atol = TOL)
    end
    return
end

function matrix_epipersepspectral2(opt)
    TOL = 1e-4
    # only functions that are decreasing
    dual_and_h = [
        (false, Hypatia.Cones.NegLogSSF()),
        (false, Hypatia.Cones.NegSqrtSSF()),
        (false, Hypatia.Cones.NegPower01SSF(3 // 10)),
        (true, Hypatia.Cones.NegLogSSF()),
        (true, Hypatia.Cones.NegEntropySSF()),
        (true, Hypatia.Cones.NegSqrtSSF()),
        (true, Hypatia.Cones.NegPower01SSF(0.5)),
    ]

    for (use_dual, h_fun) in dual_and_h
        D = (use_dual ? Dual : Prim)
        Q = Hypatia.Cones.MatrixCSqr{Float64, Float64}
        K = Hypatia.EpiPerSepSpectralCone{Float64}(h_fun, Q, 2, use_dual)
        (m, x, Q, opt_x, opt_Q) = _setup_expdesign(opt)
        opt_val = PajaritoExtras.val_or_conj(D)(eigvals(Symmetric(opt_Q, :U)), h_fun)

        JuMP.@variable(m, y)
        JuMP.@objective(m, Min, y)
        epiper = PajaritoExtras.swap_epiper(D, [y, 1]...)
        JuMP.@constraint(m, vcat(epiper..., svec(Q)) in K)
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        @test isapprox(JuMP.value.(x), opt_x, atol = TOL)
    end
    return
end

function _setup_expdesign(opt)
    m = JuMP.Model(opt)
    JuMP.@variable(m, x[1:4], Int)
    JuMP.@constraint(m, x[1:2] .>= 0.5) # avoids ill-posedness
    JuMP.@constraint(m, x[3:4] .>= 0)
    JuMP.@constraint(m, sum(x) <= 8.5)

    V = [1 1 -0.2 -0.5; 1 -1 0.5 -0.2]
    Q = V * diagm(x) * V'
    opt_x = [4, 4, 0, 0]
    opt_Q = Symmetric(V * Diagonal(opt_x) * V')
    return (m, x, Q, opt_x, opt_Q)
end

function wsosinterpnonnegative1(opt)
    TOL = 1e-4
    # real case
    dom = PolyUtils.BoxDomain{Float64}(-ones(1), ones(1))
    fs = Function[
        x -> (x - 0.5)^2 + 0.3, # min 0.3
        x -> 2x + 1.2, # min -0.8
        x -> -1 - x, # min -2
        x -> 3 - 2x^2, # min 1
    ]
    (U, pts, Ps) = PolyUtils.interpolate(dom, 2)
    pts = vec(pts)
    f_pts = [f_i.(pts) for f_i in fs]

    for use_dual in (false, true)
        opt_val = (use_dual ? -0.8 : 0.3)
        K = Hypatia.WSOSInterpNonnegativeCone{Float64, Float64}(U, Ps, use_dual)
        (m, x, y) = _setup_polymin(opt, use_dual, K, f_pts, 3.0)
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.value.(x), [1, 0, 0, 1], atol = TOL)
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        @test isapprox(JuMP.value(y), opt_val, atol = TOL)
    end
    return
end

function wsosinterpnonnegative2(opt)
    TOL = 1e-4
    # complex case
    dom = [x -> 1 - abs2(x[1])]
    fs = Function[
        x -> abs2(x - 0.5) - 0.4, # min -0.4
        x -> 2 * real(x) - imag(x) + 1, # min -1.2361
        x -> 2 * abs2(x) + imag(x), # min -0.125
        x -> real(x) + 0.1, # min -0.9
    ]
    (pts, Ps) = PolyUtils.interpolate(ComplexF64, 1, 1, dom, [1])
    pts = first.(pts)
    U = length(pts)
    f_pts = [f_i.(pts) for f_i in fs]

    for use_dual in (false, true)
        opt_val = (use_dual ? -0.9 : -0.4)
        K = Hypatia.WSOSInterpNonnegativeCone{Float64, ComplexF64}(U, Ps, use_dual)
        (m, x, y) = _setup_polymin(opt, use_dual, K, f_pts, 3.0)
        JuMP.optimize!(m)

        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.value.(x), [1, 0, 1, 0], atol = TOL)
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        @test isapprox(JuMP.value(y), opt_val, atol = TOL)
    end
    return
end

function _setup_polymin(opt, use_dual, K, f_pts, bound)
    # pick 2/4 polynomials to either
    # maximize the minimum value (primal) or minimize the maximum value (dual)
    m = JuMP.Model(opt)
    JuMP.@variable(m, x[1:4], Bin)
    JuMP.@constraint(m, sum(x) == 2)
    JuMP.@variable(m, y)

    if use_dual
        JuMP.@objective(m, Min, y)
        U = MOI.dimension(K)
        for i in 1:4
            z_i = JuMP.@variable(m, [1:U])
            JuMP.@constraint(m, z_i in K)
            JuMP.@constraint(m, sum(z_i) == 1)
            JuMP.@constraint(m, y >= JuMP.dot(f_pts[i], z_i) - bound * x[i])
        end
    else
        JuMP.@objective(m, Max, y)
        for i in 1:4
            JuMP.@constraint(m, f_pts[i] .+ (bound * (1 - x[i]) - y) in K)
        end
    end
    return (m, x, y)
end
