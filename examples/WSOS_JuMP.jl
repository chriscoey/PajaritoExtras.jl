#=
JuMP helpers for constructing formulations for WSOS constraints
=#

const AffVec = Vector{<:Union{JuMP.VariableRef, JuMP.AffExpr}}

function add_wsos(
    use_nat::Bool,
    Ps::Vector{Matrix{Float64}},
    aff::AffVec,
    model::JuMP.Model,
)
    return (use_nat ? add_wsos_nat : add_wsos_ext)(Ps, aff, model)
end

# WSOS cone formulation
function add_wsos_nat(Ps::Vector{Matrix{Float64}}, aff::AffVec, model::JuMP.Model)
    U = size(Ps[1], 1)
    @assert length(aff) == U
    K = Hypatia.WSOSInterpNonnegativeCone{Float64, Float64}(U, Ps)
    JuMP.@constraint(model, aff in K)
    return
end

# PSD extended formulation
function add_wsos_ext(Ps::Vector{Matrix{Float64}}, aff::AffVec, model::JuMP.Model)
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
