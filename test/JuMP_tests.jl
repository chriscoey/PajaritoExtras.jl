# JuMP instance tests

module TestJuMP

using Test
using LinearAlgebra
import MathOptInterface
const MOI = MathOptInterface
import JuMP
import MOIPajarito
import PajaritoExtras
import Hypatia
import Hypatia.Cones: svec_length

const rt2 = sqrt(2.0)

function svec(mat::AbstractMatrix{T}) where {T}
    vec = zeros(T, svec_length(size(mat, 1)))
    return Hypatia.Cones.smat_to_svec!(vec, mat, rt2)
end

function runtests(oa_solver, conic_solver)
    @testset "iterative method" begin
        run_jump_tests(true, oa_solver, conic_solver)
    end
    # @testset "one tree method" begin
    #     run_jump_tests(false, oa_solver, conic_solver)
    # end
    return
end

function run_jump_tests(use_iter::Bool, oa_solver, conic_solver)
    opt = JuMP.optimizer_with_attributes(
        MOIPajarito.Optimizer,
        "verbose" => true,
        # "verbose" => false,
        "use_iterative_method" => use_iter,
        "oa_solver" => oa_solver,
        "conic_solver" => conic_solver,
        "iteration_limit" => 30,
        # "time_limit" => 120.0,
    )
    insts = [_psd1, _psd2, ]#_expdesign]
    @testset "$inst" for inst in insts
        println(inst)
        inst(opt)
    end
    return
end

# TODO test complex case
function _psd1(opt)
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

function _psd2(opt)
    TOL = 1e-4
    d = 3
    mat = Symmetric(Matrix{Float64}(reshape(1:(d^2), d, d)), :U)
    λ₁ = eigmax(mat)
    m = JuMP.Model(opt)

    X = JuMP.@variable(m, [1:d, 1:d], Symmetric)
    K = Hypatia.PosSemidefTriCone{Float64, Float64}(svec_length(d))
    JuMP.@constraint(m, svec(1.0 * X) in K)
    JuMP.@objective(m, Max, JuMP.dot(mat, X))
    JuMP.@constraint(m, JuMP.tr(X) == 1)
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), λ₁, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), λ₁, atol = TOL)
    X_val = JuMP.value.(X)
    @test isapprox(tr(X_val), 1, atol = TOL)

    JuMP.set_binary.(X)
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    X_val = JuMP.value.(X)
    @test isapprox(sum(X_val), 1, atol = TOL)
    @test isapprox(X_val[d, d], 1, atol = TOL)
    return
end

# TODO use new cones
function _expdesign(opt)
    TOL = 1e-4
    # experiment design
    V = [1 1 -0.2 -0.5; 1 -1 0.5 -0.2]
    function setup_exp_design()
        m = JuMP.Model(opt)
        JuMP.@variable(m, x[1:4], Int)
        JuMP.@constraint(m, x[1:2] .>= 1) # avoids ill-posedness
        JuMP.@constraint(m, x[3:4] .>= 0)
        JuMP.@constraint(m, sum(x) <= 8)
        Q = V * diagm(x) * V'
        return (m, x, Q)
    end

    # A-optimal
    (m, x, Q) = setup_exp_design()
    JuMP.set_start_value.(x, [4, 4, 0, 0]) # partial warm start
    JuMP.@variable(m, y[1:2])
    JuMP.@objective(m, Min, sum(y))
    for i in 1:2
        ei = zeros(2)
        ei[i] = 1
        Qyi = [Q ei; ei' y[i]]
        JuMP.@constraint(m, Symmetric(Qyi) in JuMP.PSDCone())
    end
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), 1 / 4, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), 1 / 4, atol = TOL)
    x_val = JuMP.value.(x)
    @test isapprox(x_val, [4, 4, 0, 0], atol = TOL)
    @test isapprox(JuMP.value.(y[1]), JuMP.value.(y[2]), atol = TOL)

    # E-optimal
    (m, x, Q) = setup_exp_design()
    JuMP.@variable(m, y)
    JuMP.set_start_value.(vcat(x, y), [4, 4, 0, 0, 8]) # full warm start
    JuMP.@objective(m, Max, y)
    Qy = Q - y * Matrix(I, 2, 2)
    JuMP.@constraint(m, Symmetric(Qy) in JuMP.PSDCone())
    JuMP.optimize!(m)
    @test JuMP.termination_status(m) == MOI.OPTIMAL
    @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
    @test isapprox(JuMP.objective_value(m), 8, atol = TOL)
    @test isapprox(JuMP.objective_bound(m), 8, atol = TOL)
    x_val = JuMP.value.(x)
    @test isapprox(x_val, [4, 4, 0, 0], atol = TOL)

    # D-optimal
    for use_logdet in (true, false)
        opt_x = [4, 4, 0, 0]
        opt_Q = Symmetric(V * Diagonal(opt_x) * V')
        (m, x, Q) = setup_exp_design()
        JuMP.@variable(m, y)
        JuMP.@objective(m, Max, y)
        Qvec = [Q[1, 1], Q[2, 1], Q[2, 2]]
        if use_logdet
            JuMP.@constraint(m, vcat(y, 1.0, Qvec) in MOI.LogDetConeTriangle(2))
            opt_val = logdet(opt_Q)
        else
            JuMP.@constraint(m, vcat(y, Qvec) in MOI.RootDetConeTriangle(2))
            opt_val = sqrt(det(opt_Q))
        end
        JuMP.optimize!(m)
        @test JuMP.termination_status(m) == MOI.OPTIMAL
        @test JuMP.primal_status(m) == MOI.FEASIBLE_POINT
        @test isapprox(JuMP.objective_value(m), opt_val, atol = TOL)
        @test isapprox(JuMP.objective_bound(m), opt_val, atol = TOL)
        x_val = JuMP.value.(x)
        @test isapprox(x_val, opt_x, atol = TOL)
    end
    return
end

end
