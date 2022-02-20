#=
code taken and modified from PiecewiseLinearOpt.jl
https://github.com/joehuchette/PiecewiseLinearOpt.jl
Copyright (c) 2016: Joey Huchette

We homogenize using a given perspective variable v. σ ≥ 0 and sum(σ) = v are already imposed.
In the "bounded" case, we assume v ∈ [0,1] (a MIP formulation can be linear/big-M).
In the "unbounded" case, we assume v ≥ 0 (a MIP formulation needs conic constraints).
=#

import PiecewiseLinearOpt
const PLO = PiecewiseLinearOpt

# SOS2 formulation types
abstract type PWLSOS2 end
struct SOS2 <: PWLSOS2 end
struct CCBounded <: PWLSOS2 end
struct CCUnbounded <: PWLSOS2 end
struct LogIBBounded <: PWLSOS2 end
struct LogIBUnbounded <: PWLSOS2 end
# struct MCBounded <: PWLSOS2 end # TODO invalid?
# struct LogBounded <: PWLSOS2 end # TODO invalid?
# struct ZigZagBounded <: PWLSOS2 end # TODO invalid?

# bounded formulations

function add_PWL(::SOS2, model::JuMP.Model, σ::Vector{JuMP.VariableRef}, ::JuMPScalar)
    JuMP.@constraint(model, σ in JuMP.SOS2())
    return
end

function add_PWL(
    ::CCBounded,
    model::JuMP.Model,
    σ::Vector{JuMP.VariableRef},
    v::JuMPScalar, # binary or 1
)
    n = length(σ)
    y = JuMP.@variable(model, [1:(n - 1)], Bin)
    JuMP.@constraints(model, begin
        sum(y) == v
        σ[1] <= y[1]
        [i in 2:(n - 1)], σ[i] <= y[i - 1] + y[i]
        σ[n] <= y[n - 1]
    end)
    return
end

function add_PWL(
    ::CCUnbounded,
    model::JuMP.Model,
    σ::Vector{JuMP.VariableRef},
    v::JuMPScalar, # binary or 1
)
    n = length(σ)
    y = JuMP.@variable(model, [1:(n - 1)], Bin)
    z = JuMP.@variable(model) # TODO or add n variables?
    RSOC = JuMP.RotatedSecondOrderCone()
    JuMP.@constraints(model, begin
        sum(y) == v
        [z, y[1], σ[1]] in RSOC
        [i in 2:(n - 1)], [z, y[i - 1] + y[i], σ[i]] in RSOC
        [z, y[n - 1], σ[n]] in RSOC
    end)
    return
end

function _logIB_H(k::Int)
    _H = PLO.reflected_gray_codes(k)
    d = length(_H)
    H = Dict(i => _H[i] for i in 1:d)
    H[0] = H[1]
    H[d + 1] = H[d]
    return H
end

function add_PWL(
    ::LogIBBounded,
    model::JuMP.Model,
    σ::Vector{JuMP.VariableRef},
    v::JuMPScalar, # binary or 1
)
    n = length(σ)
    k = ceil(Int, log2(n - 1))
    y = JuMP.@variable(model, [1:k], Bin)
    H = _logIB_H(k)
    for j in 1:k
        H1 = JuMP.@expression(model, sum(σ[i] for i in 1:n if H[i - 1][j] == H[i][j] == 1))
        H0 = JuMP.@expression(model, sum(σ[i] for i in 1:n if H[i - 1][j] == H[i][j] == 0))
        JuMP.@constraint(model, H1 <= y[j])
        JuMP.@constraint(model, H0 <= v - y[j])
    end
    return
end

function add_PWL(
    ::LogIBUnbounded,
    model::JuMP.Model,
    σ::Vector{JuMP.VariableRef},
    v::JuMPScalar,
)
    n = length(σ)
    k = ceil(Int, log2(n - 1))
    y = JuMP.@variable(model, [1:k], Bin)
    z = JuMP.@variable(model) # TODO or add k variables?
    RSOC = JuMP.RotatedSecondOrderCone()
    H = _logIB_H(k)
    for j in 1:k
        H1 = JuMP.@expression(model, sum(σ[i] for i in 1:n if H[i - 1][j] == H[i][j] == 1))
        H0 = JuMP.@expression(model, sum(σ[i] for i in 1:n if H[i - 1][j] == H[i][j] == 0))
        JuMP.@constraint(model, [z, y[j], H1] in RSOC)
        JuMP.@constraint(model, [z, v - y[j], H0] in RSOC)
    end
    return
end

# function add_PWL(
#     ::MCBounded,
#     model::JuMP.Model,
#     σ::Vector{JuMP.VariableRef},
#     v::JuMPScalar,
# )
#     n = length(σ)
#     γ = JuMP.@variable(model, [1:(n - 1), 1:n])
#     y = JuMP.@variable(model, [1:(n - 1)], Bin)
#     JuMP.@constraints(model, begin
#         sum(y) == v
#         sum(γ[i, :] for i in 1:(n - 1)) .== σ
#         [i in 1:(n - 1)], γ[i, i] + γ[i, i + 1] >= y[i]
#     end)
#     return
# end

# function add_PWL(
#     pwl::Union{LogBounded, ZigZagBounded},
#     model::JuMP.Model,
#     σ::Vector{JuMP.VariableRef},
#     ::JuMPScalar,
# )
#     n = length(σ)
#     k = ceil(Int, log2(n - 1))
#     y = JuMP.@variable(model, [1:k], Bin)

#     (codes, planes) = _encoding(pwl)
#     PLO.sos2_encoding_constraints!(model, σ, y, codes(k), planes(k))
#     return
# end

# _encoding(::Log) = (PLO.reflected_gray_codes, PLO.unit_vector_hyperplanes)
# _encoding(::ZigZag) = (PLO.zigzag_codes, PLO.zigzag_hyperplanes)
