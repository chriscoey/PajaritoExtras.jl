#=
JuMP helpers for constructing approximations/relaxations of cone boundary type constraints
for vector function cones, using extended formulations with 3-dimensional cones

TODO only concerned with part of cone boundary for which the nonlinear inequality defining
the cone is at equality

options:
use_equality::Bool - use equality condition in nonconvex constraint, else inequality
add_convex::Bool - add the convex nonlinear constraint (only if use_equality is false)
=#

# for vector separable spectral functions
function add_nonconvex(
    f::VecSpecExt,
    aff::Vector{<:JuMPScalar},
    model::JuMP.Model,
    use_equality::Bool,
    add_convex::Bool,
    pts_min::Float64,
    pts_max::Float64,
    pwl::PWLSOS2,
    num_pts::Int,
)
    if add_convex
        @assert !use_equality
    end

    # epigraph and perspective variables are swapped if dual cone is used
    (epi, per) = swap_if_dual(aff[1], aff[2], f)
    w = aff[3:end]
    d = length(w)

    # TODO could be redundant with other constraints eg binary
    JuMP.@constraint(model, per >= 0)

    # EF variables and linear constraint
    λ = JuMP.@variable(model, [1:d])
    if use_equality
        JuMP.@constraint(model, epi == sum(λ))
    else
        JuMP.@constraint(model, epi <= sum(λ))
    end

    # piecewise linear data
    if is_domain_pos(f)
        @assert pts_min >= 0
    end
    pts = collect(range(pts_min, pts_max, length = num_pts))
    f_pts = [get_val([pt], f) for pt in pts]

    # 3-dim cone constraints
    for i in 1:d
        aff = [λ[i], per, w[i]]

        σ = add_PWL_3D(aff, model, pts, f_pts)
        add_PWL(pwl, model, σ, per)

        if add_convex
            # TODO check why add_spectral function does swap of epi/per again
            add_spectral(f, 1, aff, model)
        end
    end
    return
end

function add_PWL_3D(
    aff::Vector{<:JuMPScalar},
    model::JuMP.Model,
    pts::Vector{Float64},
    f_pts::Vector{Float64},
)
    @assert length(aff) == 3
    σ = JuMP.@variable(model, [1:length(pts)], lower_bound = 0)
    JuMP.@constraints(model, begin
        dot(σ, f_pts) == aff[1]
        sum(σ) == aff[2]
        dot(σ, pts) == aff[3]
    end)
    return σ
end
