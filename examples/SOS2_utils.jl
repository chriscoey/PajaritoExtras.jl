#=
some code taken and modified from PiecewiseLinearOpt.jl
https://github.com/joehuchette/PiecewiseLinearOpt.jl
Copyright (c) 2016: Joey Huchette

We homogenize using a given perspective variable v. σ ≥ 0 and sum(σ) = v are already imposed.
In the "bounded" case, we assume v ∈ [0,1] (MIP formulation can be linear/big-M).
In the "unbounded" case, we assume v ≥ 0 (MIP formulation needs conic constraints).
=#

import PiecewiseLinearOpt
const PLO = PiecewiseLinearOpt

const FloatOrVar = Union{Float64, JuMP.VariableRef}

# SOS2 formulation types
abstract type PWLSOS2 end

struct SOS2 <: PWLSOS2 end
struct CCBounded <: PWLSOS2 end
struct CCUnbounded <: PWLSOS2 end
struct LogIBBounded <: PWLSOS2 end
struct LogIBUnbounded <: PWLSOS2 end

function PWL_SOS2(::SOS2, model::JuMP.Model, σ::Vector{JuMP.VariableRef}, ::FloatOrVar)
    JuMP.@constraint(model, σ in JuMP.SOS2())
    return
end

function PWL_SOS2(
    ::CCBounded,
    model::JuMP.Model,
    σ::Vector{JuMP.VariableRef},
    binper::FloatOrVar, # binary or 1
)
    n = length(σ)
    y = JuMP.@variable(model, [1:(n - 1)], Bin)
    JuMP.@constraints(model, begin
        sum(y) == binper
        σ[1] <= y[1]
        [i in 2:(n - 1)], σ[i] <= y[i - 1] + y[i]
        σ[n] <= y[n - 1]
    end)
    return
end

function PWL_SOS2(
    ::CCUnbounded,
    model::JuMP.Model,
    σ::Vector{JuMP.VariableRef},
    binper::FloatOrVar, # binary or 1
)
    n = length(σ)
    y = JuMP.@variable(model, [1:(n - 1)], Bin)
    z = JuMP.@variable(model) # TODO or add n variables?
    RSOC = JuMP.RotatedSecondOrderCone()
    JuMP.@constraints(model, begin
        sum(y) == binper
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

function PWL_SOS2(
    ::LogIBBounded,
    model::JuMP.Model,
    σ::Vector{JuMP.VariableRef},
    binper::FloatOrVar, # binary or 1
)
    n = length(σ)
    k = ceil(Int, log2(n - 1))
    y = JuMP.@variable(model, [1:k], Bin)
    H = _logIB_H(k)
    for j in 1:k
        H1 = JuMP.@expression(model, sum(σ[i] for i in 1:n if H[i - 1][j] == H[i][j] == 1))
        H0 = JuMP.@expression(model, sum(σ[i] for i in 1:n if H[i - 1][j] == H[i][j] == 0))
        JuMP.@constraint(model, H1 <= y[j])
        JuMP.@constraint(model, H0 <= binper - y[j])
    end
    return
end

function PWL_SOS2(
    ::LogIBUnbounded,
    model::JuMP.Model,
    σ::Vector{JuMP.VariableRef},
    binper::FloatOrVar,
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
        JuMP.@constraint(model, [z, binper - y[j], H0] in RSOC)
    end
    return
end
