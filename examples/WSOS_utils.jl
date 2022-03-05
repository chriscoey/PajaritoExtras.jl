#=
JuMP helpers for constructing formulations for WSOS constraints
=#

# generate random polynomials nonegative on domain (may fail if domain not bounded)
function make_nonneg_polys(num::Int, Ps::Vector{Matrix{Float64}})
    # build JuMP model for finding minima (primal polymin)
    U = size(Ps[1], 1)
    K = Hypatia.WSOSInterpNonnegativeCone{Float64, Float64}(U, Ps)
    model = JuMP.Model(sep_hypatia)
    JuMP.@variable(model, a)
    JuMP.@objective(model, Max, a)
    con = JuMP.@constraint(model, ones(U) .- a in K)

    # generate polys
    return [make_nonneg_poly(U, model, con) for _ in 1:num]
end

function make_nonneg_poly(U::Int, model::JuMP.Model, con::JuMP.ConstraintRef)
    vals = Float64.(rand(0:9, U))
    MOI.modify(JuMP.backend(model), JuMP.index(con), MOI.VectorConstantChange(vals))

    JuMP.optimize!(model)
    stat = JuMP.termination_status(model)

    if stat in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
        # make poly nonnegative if minimum is negative
        min_val = JuMP.dual_objective_value(model)
        if min_val < 1e-5
            vals .+= abs(min_val)
        end
        return vals
    else
        # try again
        @warn("try again")
        return make_nonneg_poly(U, model, con)
    end
end

# find minimum value of a polynomial - solve primal polymin
function find_poly_min(vals::Vector{Float64}, Ps::Vector{Matrix{Float64}})
    U = size(Ps[1], 1)
    K = Hypatia.WSOSInterpNonnegativeCone{Float64, Float64}(U, Ps)
    model = JuMP.Model(sep_hypatia)
    JuMP.@variable(model, a)
    JuMP.@objective(model, Max, a)
    JuMP.@constraint(model, vals .- a in K)

    JuMP.optimize!(model)
    @assert JuMP.termination_status(model) in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL)
    return JuMP.dual_objective_value(model)
end

# check whether polynomial is nonnegative (actually SOS) - solve separation problem
function check_nonneg(vals::Vector{Float64}, Ps::Vector{Matrix{Float64}})
    U = size(Ps[1], 1)
    K = Hypatia.WSOSInterpNonnegativeCone{Float64, Float64}(U, Ps)
    model = JuMP.Model(sep_hypatia)
    JuMP.@constraint(model, vals in K)

    JuMP.optimize!(model)
    return (JuMP.termination_status(model) in (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL))
end

# add a WSOS constraint
function add_wsos(
    use_nat::Bool,
    Ps::Vector{Matrix{Float64}},
    aff::AbstractVector{<:JuMPScalar},
    model::JuMP.Model,
)
    return (use_nat ? add_wsos_nat : add_wsos_ext)(Ps, aff, model)
end

# WSOS cone formulation
function add_wsos_nat(
    Ps::Vector{Matrix{Float64}},
    aff::AbstractVector{<:JuMPScalar},
    model::JuMP.Model,
)
    U = size(Ps[1], 1)
    @assert length(aff) == U
    K = Hypatia.WSOSInterpNonnegativeCone{Float64, Float64}(U, Ps)
    JuMP.@constraint(model, aff in K)
    return
end

# PSD extended formulation
function add_wsos_ext(
    Ps::Vector{Matrix{Float64}},
    aff::AbstractVector{<:JuMPScalar},
    model::JuMP.Model,
)
    U = size(Ps[1], 1)
    @assert length(aff) == U

    psd_vars = []
    for Pr in Ps
        Lr = size(Pr, 2)
        psd_r = JuMP.@variable(model, [1:Lr, 1:Lr], Symmetric)
        if Lr == 1
            JuMP.@constraint(model, psd_r[1, 1] >= 0)
        else
            JuMP.@constraint(model, psd_r in JuMP.PSDCone())
        end
        push!(psd_vars, psd_r)
    end

    coeffs_lhs = JuMP.@expression(
        model,
        [u in 1:U],
        sum(
            sum(
                Pr[u, k] * Pr[u, l] * psd_r[k, l] * (k == l ? 1 : 2) for
                k in 1:size(Pr, 2) for l in 1:k
            ) for (Pr, psd_r) in zip(Ps, psd_vars)
        )
    )
    JuMP.@constraint(model, coeffs_lhs .== aff)
    return
end
