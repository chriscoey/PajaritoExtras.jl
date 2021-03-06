#=
design a modular device using parts from a catalog:
- for each module we pick one corresponding part from the catalog, satisfying a constraint
associated with the part
- these constraints involve certain global design variables (e.g. weight, size,
lifetime, battery power)
- we minimize the cost of the parts plus a linear function of the global design variables

sets:
i ∈ 1..m          modules
j ∈ J[i], i ∈ m   parts for module i
k ∈ 1..n          global design variables

variables:
xᵢⱼ ∈ {0,1}   whether part j ∈ J[i] is selected for module i
y   ∈ ℝ₊ⁿ     global design characteristics

objective:
min ∑ᵢⱼ cᵢⱼ xᵢⱼ + ∑ₖ dₖ yₖ   where c ≥ 0

constraints:
∑ⱼ xᵢⱼ = 1, ∀i                 disjunction for module i
ymin ≤ y ≤ ymax                bounds
(xᵢⱼ = 1) ⟹ (y ∈ Sᵢⱼ), ∀i,j   constraint for part i,j where Sᵢⱼ is some set

use the conic copies-of-variables formulation for the disjunctions
=#

struct ModularDesign <: ExampleInstance
    m::Int # number of modules
    p::Int # number of part options per module
    n::Int # number of design variables
    use_nat::Bool # for convex case, use natural formulations
    use_nonconvex::Bool # use nonconvex equality constraints, else convex constraints
    pwl::PWLSOS2 # type of piecewise linear SOS2 formulation to use
    num_pts::Int # maximum number of points per piecewise linearization
end

function ModularDesign(m::Int, p::Int, n::Int, use_nat::Bool)
    return ModularDesign(m, p, n, use_nat, false, SOS2(), 0)
end

function build(inst::ModularDesign)
    (m, p, n) = (inst.m, inst.p, inst.n)
    @assert m >= 1
    @assert p >= 1
    @assert n >= 2
    ymin = 1e-3
    ymax = 2.0
    @assert ymax > ymin > 0
    q = 1 + div(n, 2) # dimension of convex epigraph sets
    @assert q >= 2

    # generate data
    c = 10 * rand(m, p)
    d = vec(sprandn(n, 1, 0.3))
    # random solution to ensure feasibility
    y0 = (ymax - ymin) * rand(n) .+ ymin

    # functions and affine data defining the conic constraints
    # TODO -log has high slope near zero, so need non-evenly-spaced points
    if inst.use_nat
        f_list = VecSepSpecPrim[
            # VecNegLog(),
            VecNegEntropy(),
            VecNegSqrt(),
            VecNegPower01(0.8),
            VecPower12(2.0),
        ]
    else
        f_list = VecSepSpecPrimEF[
            # VecNegLogEF(),
            VecNegEntropyEF(),
            VecNegSqrtEF(),
            VecNegPower01EF(0.8),
            VecPower12EF(2.0),
        ]
    end
    fs = rand(f_list, m, p)
    Gs = Matrix{Matrix{Float64}}(undef, m, p)
    hs = Matrix{Vector{Float64}}(undef, m, p)
    G_sparsity = min(0.7, 4 / (q * n))
    for i in 1:m, j in 1:p
        f_ij = fs[i, j]
        G_ij = sprandn(q, n, G_sparsity)
        for k in 2:q
            if iszero(G_ij[k, :])
                G_ij[k, rand(1:n)] = randn()
            end
        end
        h_ij = G_ij * y0
        w = rand(q - 1) .+ 1e-3
        h_ij[2:end] .+= w
        h_ij[1] += get_val(w, f_ij)
        if !inst.use_nonconvex
            h_ij[1] += rand()
        end
        Gs[i, j] = G_ij
        hs[i, j] = h_ij
    end

    # build model
    model = JuMP.Model()
    JuMP.@variable(model, y[1:n])
    JuMP.@variable(model, x[1:m, 1:p], Bin)
    JuMP.@objective(model, Min, dot(c, x) + dot(d, y))

    # disjunctive copies-of-variables formulation
    # x[i, j] switches on conic constraint for part j of module i
    JuMP.@variable(model, z[1:m, 1:p, 1:n])
    JuMP.@constraint(model, [i in 1:m], sum(x[i, :]) == 1)
    JuMP.@constraint(model, [i in 1:m], sum(z[i, j, :] for j in 1:p) .== y)

    # x[i, j] = 0 ⇒ z[i, j, :] = 0, x[i, j] = 1 ⇒ z[i, j, :] = y
    affs = Matrix{Vector{JuMP.AffExpr}}(undef, m, p)
    for i in 1:m, j in 1:p
        x_ij = x[i, j]
        z_ij = z[i, j, :]
        JuMP.@constraint(model, z_ij .<= ymax * x_ij)
        JuMP.@constraint(model, ymin * x_ij .<= z_ij)

        h_ij = hs[i, j]
        G_ij = Gs[i, j]
        aff = affs[i, j] = h_ij * x_ij - G_ij * z_ij
        u_aff = aff[1]
        w_aff = aff[2:end]
        K = length(w_aff)
        f = fs[i, j]

        if !inst.use_nonconvex
            # add convex constraint
            add_spectral(f, K, vcat(u_aff, x_ij, w_aff), model)
            continue
        end

        # relaxation of nonconvex constraint
        λ = JuMP.@variable(model, [1:K])
        JuMP.@constraint(model, u_aff == sum(λ))

        for k in 1:K
            # auxiliary variable for w_aff[k]
            w_k = JuMP.@variable(model)
            JuMP.@constraint(model, w_k == w_aff[k])

            # convex constraint (above graph)
            add_spectral(f, 1, vcat(λ[k], x_ij, 1.0 * w_k), model)

            # bounds
            G_k = G_ij[1 + k, :]
            min_ij = h_ij[1 + k] - sum((g < 0 ? ymin : ymax) * g for g in G_k)
            max_ij = h_ij[1 + k] - sum((g < 0 ? ymax : ymin) * g for g in G_k)
            @assert min_ij <= max_ij
            max_ij = max(1e-4, max_ij)
            min_ij = max(1e-4, min_ij)

            # interpolation
            pts = collect(range(min_ij, max_ij, length = inst.num_pts))
            f_pts = [get_val([pt], f) for pt in pts]

            # auxiliary variables and SOS2 formulation
            σ = JuMP.@variable(model, [1:length(pts)], lower_bound = 0)
            PWL_SOS2(inst.pwl, model, σ, x_ij)

            # data constraints
            JuMP.@constraints(model, begin
                dot(σ, f_pts) >= λ[k] # below PWL graph
                sum(σ) == x_ij
                dot(σ, pts) == w_k
            end)
        end
    end

    # save for use in tests
    model.ext[:x] = x
    model.ext[:y] = y
    model.ext[:z] = z
    model.ext[:fs] = fs
    model.ext[:affs] = affs

    return model
end

function test_extra(inst::ModularDesign, model::JuMP.Model)
    stat = JuMP.termination_status(model)
    @test stat in OPT_OR_LIMIT
    (stat == MOI.OPTIMAL) || return

    tol = eps()^0.2
    x = JuMP.value.(model.ext[:x])
    x_int = round.(Int, x)
    @test x ≈ x_int atol = tol rtol = tol
    @test all(-tol .<= x .<= 1 + tol)

    # check feasibility for disjunctive linear constraints
    (m, p) = (inst.m, inst.p)
    y = JuMP.value.(model.ext[:y])
    z = JuMP.value.(model.ext[:z])
    @test all(sum(x_int[i, :]) == 1 for i in 1:m)
    z0ij(i, j) = (iszero(x_int[i, j]) ? zero(y) : y)
    approxij(i, j) = isapprox(z[i, j, :], z0ij(i, j), atol = tol, rtol = tol)
    @test all(approxij(i, j) for i in 1:m for j in 1:p)

    # check feasibility for nonlinear constraints
    fs = model.ext[:fs]
    affs = model.ext[:affs]
    viol(fvaff) = begin
        (f, v, aff) = fvaff
        (u, w) = (aff[1], aff[2:end])
        iszero(v) && return -u
        @assert v == 1
        any(<(-tol), w) && return NaN
        return get_val(max.(w, 1e-9), f) - u
    end
    if inst.use_nonconvex
        tol_loose = eps()^0.1
        # check approximate equality (and inside cone)
        conicfeas = (a -> -tol_loose < a < tol)
    else
        # convex constraints: check inside cone
        conicfeas = <(tol)
    end
    ϵs = [viol((fs[i, j], x_int[i, j], JuMP.value.(affs[i, j]))) for i in 1:m, j in 1:p]
    @show maximum(abs, ϵs) # TODO delete
    @test all(conicfeas, ϵs)
    return
end
