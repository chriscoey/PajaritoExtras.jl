
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
        mat = Hermitian(smat(R, x_val), :U)
        @test isapprox(tr(mat), 1, atol = TOL)
        @test eigmin(mat) >= -TOL

        JuMP.set_binary.(x_diag)
        JuMP.optimize!(m)
        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        x_val = JuMP.value.(x)
        mat = Hermitian(smat(R, x_val), :U)
        @test isapprox(tr(mat), 1, atol = TOL)
        @test eigmin(mat) >= -TOL
        @test isapprox(mat[d, d], 1, atol = TOL)
    end
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

function epipersquare1(opt)
    TOL = 1e-4
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
