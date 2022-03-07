#=
a ball-packing type problem: choose centers and radii of m 2-norm balls into a box [0, 1]ⁿ in
n-dimensional space, maximizing the sum of the radii

constraints are nonconvex, so use SOS2 or PWL techniques with extended formulations
=#

struct BallPacking <: ExampleInstance
    m::Int # number of balls
    n::Int # dimension of space
    pwl::PWLSOS2 # type of piecewise linear SOS2 formulation to use
    num_pts::Int # number of points per piecewise linearization
end

function build(inst::BallPacking)
    (m, n, num_pts) = (inst.m, inst.n, inst.num_pts)
    @assert m >= 2
    @assert n >= 2
    @assert isodd(num_pts) # ensures 0 is included
    pts = collect(range(-1.0, 1.0, length = num_pts))
    @assert 0.0 in pts
    f_pts = abs2.(pts)

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, 0 <= R[1:m] <= 0.5) # radii
    JuMP.@variable(model, 0 <= C[1:m, 1:n] <= 1) # centers
    JuMP.@objective(model, Max, sum(R))

    JuMP.@constraints(model, begin
        # balls inside box
        [i in 1:m, k in 1:n], C[i, k] >= R[i]
        [i in 1:m, k in 1:n], C[i, k] + R[i] <= 1
        # break some symmetry TODO more
        [i in 1:(m - 1)], C[i, 1] <= C[i + 1, 1]
    end)

    # for each pair of balls i < j, need Rᵢ + Rⱼ ≤ ‖Cᵢ - Cⱼ‖
    for i in 1:m, j in (i + 1):m
        # EF variables and linear constraint
        λ = JuMP.@variable(model, [1:n], lower_bound = 0)
        Rij = JuMP.@variable(model)
        Cij = JuMP.@variable(model, [1:n])
        JuMP.@constraints(model, begin
            Rij == R[i] + R[j]
            Rij <= sum(λ)
            Cij .== C[i, :] - C[j, :]
        end)

        # 3-dim cone constraints; note Cij ∈ [-1, 1]ⁿ, ‖Cij‖ ≤ √n, Rij ∈ [0, 1]
        for k in 1:n
            # auxiliary variables and SOS2 formulation
            σ = JuMP.@variable(model, [1:length(pts)], lower_bound = 0)
            PWL_SOS2(inst.pwl, model, σ, 1.0)

            # data constraints
            JuMP.@constraints(model, begin
                dot(σ, f_pts) == λ[k]
                sum(σ) == Rij
                dot(σ, pts) == Cij[k]
            end)
        end
    end

    # save for use in tests
    model.ext[:R] = R
    model.ext[:C] = C

    return model
end

function test_extra(inst::BallPacking, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat in OPT_OR_LIMIT
    (stat == MOI.OPTIMAL) || return

    R = JuMP.value.(model.ext[:R])
    C = JuMP.value.(model.ext[:C])
    # check feasibility for linear constraints
    tol = eps()^0.4
    @test all(-tol .<= R .<= 0.5 + tol)
    @test all(-tol .<= C .<= 1 + tol)
    # check near-feasibility for nonconvex constraints
    tol_loose = eps()^0.1
    approxij(i, j) = (R[i] + R[j] <= norm(C[i, :] - C[j, :]) + tol_loose)
    @test all(approxij(i, j) for i in 1:(inst.m) for j in 1:(i - 1))
    return
end
