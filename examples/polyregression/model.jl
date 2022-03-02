#=
assume two populations and assign observations to the populations and simultaneously fit
polynomial functions to these populations; on a given domain containing the observations,
the functions must satisfy certain shape constraints

given domain D ⊆ ℝⁿ, observations (Xᵢ ∈ D, Yᵢ ∈ ℝ₊), solve for p₁(x), p₂(x):
min  ‖z‖ :                              residual loss
p₂(x) ≥ 0, ∀x ∈ D                       p2 nonnegative
p₁(x) ≥ p₂(x), ∀x ∈ D                   p1 above p2
1 ≥ p₁(x), ∀x ∈ D                       p1 below 1
bᵢ ∈ {0,1}, ∀i                          assignment of points
0.4m <= Σᵢbᵢ <= 0.6m                    approximately balance assignments
|zᵢ - (Yᵢ - p₁(Xᵢ))| ≤ Mᵢ(1 - bᵢ), ∀i   residuals for p1
|zᵢ - (Yᵢ - p₂(Xᵢ))| ≤ Mᵢbᵢ, ∀i         residuals for p2

TODO
- complex case
=#

struct PolyRegression <: ExampleInstance
    n::Int # dimension of independent variables / polynomial domain
    halfdeg::Int # polynomial regressor half-degree
    m::Int # number of observations
    signal_ratio::Float64 # signal-to-noise ratio (zero if no noise)
    use_nat::Bool # use WSOS cone formulation, else SDP formulation
end

function build(inst::PolyRegression)
    (n, m, signal_ratio) = (inst.n, inst.m, inst.signal_ratio)
    @assert m > n >= 1

    # setup interpolation
    D = PolyUtils.BoxDomain{Float64}(zeros(n), ones(n))
    (U, _, Ps, V) = PolyUtils.interpolate(D, inst.halfdeg, calc_V = true)
    @assert m > 2U
    F = qr!(Array(V'), ColumnNorm())

    # generate noisy data in D = [0, 1]ⁿ from two underlying functions
    X = rand(m, n)
    V_X = PolyUtils.make_chebyshev_vandermonde(X, 2 * inst.halfdeg)
    p_X = F \ V_X'

    fs = make_nonneg_polys(2, Ps)
    f1 = fs[1]
    f2 = fs[2]
    @assert minimum(p_X' * f1) >= -1e-4
    @assert minimum(p_X' * f2) >= -1e-4

    f1 .+= f2
    f1_max = -find_poly_min(-f1, Ps)
    @assert f1_max >= -1e-4
    f1 ./= f1_max
    f2 ./= f1_max

    m1 = div(m, 2)
    Y1 = p_X' * f1
    Y2 = p_X' * f2
    @assert minimum(Y1) >= -1e-4
    @assert maximum(Y1) <= 1 + 1e-4
    @assert minimum(Y2) >= -1e-4
    @assert maximum(Y2) <= 1 + 1e-4
    Y = vcat(Y1[1:m1], Y2[(m1 + 1):m])

    if !iszero(signal_ratio)
        noise = randn(m)
        noise .*= norm(Y) / sqrt(signal_ratio) / norm(noise)
        Y .+= noise
        Y = max.(Y, 0.0)
        Y = min.(Y, 1.0)
    end
    @assert all(0 .<= Y .<= 1)

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, p1[1:U])
    JuMP.@variable(model, p2[1:U])
    JuMP.@variable(model, r1[1:m])
    JuMP.@variable(model, r2[1:m])
    JuMP.@variable(model, z[1:m])
    JuMP.@variable(model, b[1:m], Bin)

    JuMP.@variable(model, epi)
    JuMP.@objective(model, Min, epi)

    # TODO tighten big-Ms, inefficient sparsity
    JuMP.@constraints(model, begin
        vcat(epi, z) in MOI.SecondOrderCone(1 + m)
        r1 .== z - Y + p_X' * p1
        r2 .== z - Y + p_X' * p2
        ceil(0.4 * m) <= sum(b)
        sum(b) <= floor(0.6 * m)
        -(1 .- b) .<= r1
        r1 .<= 1 .- b
        -b .<= r2
        r2 .<= b
    end)

    # shape constraints
    add_wsos(aff) = (inst.use_nat ? add_wsos_nat : add_wsos_ext)(Ps, aff, model)
    add_wsos(p2)
    add_wsos(p1 - p2)
    add_wsos(1 .- p1)

    # save for use in tests
    model.ext[:b] = b
    model.ext[:p1] = p1
    model.ext[:p2] = p2
    model.ext[:Ps] = Ps

    return model
end

function test_extra(inst::PolyRegression, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat == MOI.OPTIMAL
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    b = JuMP.value.(model.ext[:b])
    @test b ≈ round.(Int, b) atol = tol rtol = tol
    @test all(-tol .<= b .<= 1 + tol)
    # check WSOS feasibility
    Ps = model.ext[:Ps]
    p1 = JuMP.value.(model.ext[:p1])
    p2 = JuMP.value.(model.ext[:p2])
    @test check_nonneg(p2 .+ tol, Ps)
    @test check_nonneg(p1 - p2 .+ tol, Ps)
    @test check_nonneg((1 + tol) .- p1, Ps)
    # objective is 0 if no noise
    if iszero(inst.signal_ratio)
        @test JuMP.objective_value(model) ≈ 0 atol = tol rtol = tol
    end
    return
end
