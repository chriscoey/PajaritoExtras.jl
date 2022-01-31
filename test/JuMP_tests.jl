# JuMP instance tests

module TestJuMP

using Test
using LinearAlgebra
import MathOptInterface
const MOI = MathOptInterface
import JuMP
import MOIPajarito
import PajaritoExtras
import PajaritoExtras: Prim, Dual, svec_idx
import Hypatia
import Hypatia.Cones: vec_length, svec_length, vec_copyto!
import Hypatia.PolyUtils

# all instances
inst_all = String[
    "possemideftri1",
    "possemideftri2",
    "possemideftrisparse1",
    "possemideftrisparse2",
    "epinormeucl1",
    "epinormeucl2",
    "epipersquare1",
    "epipersquare2",
    "epinorminf1",
    "epinorminf2",
    "epinormspectraltri1",
    "epinormspectraltri2",
    "epinormspectral1",
    "epinormspectral2",
    "hypogeomean1",
    "hypogeomean2",
    "hyporootdettri1",
    "hyporootdettri2",
    "vector_epipersepspectral1",
    "vector_epipersepspectral2",
    "matrix_epipersepspectral1",
    "matrix_epipersepspectral2",
    "wsosinterpnonnegative1",
    "wsosinterpnonnegative2",
]

# instances to test with extended formulations option off
inst_noextend = String[
    "epinormeucl1",
    "epinormeucl2",
    "epipersquare1",
    "epipersquare2",
    "epinorminf1",
    "epinorminf2",
    "hypogeomean1",
    "hypogeomean2",
    "vector_epipersepspectral1",
    "vector_epipersepspectral2",
]

function runtests(oa_solver)
    @testset "iterative" begin
        @info "iterative, all"
        run_jump_tests(inst_all, true, true, oa_solver)
        @info "iterative, not extended"
        run_jump_tests(inst_noextend, false, true, oa_solver)
    end
    @testset "one tree" begin
        @info "one tree, all"
        run_jump_tests(inst_all, true, false, oa_solver)
        @info "one tree, not extended"
        run_jump_tests(inst_noextend, false, false, oa_solver)
    end
    return
end

function run_jump_tests(
    insts::Vector{String},
    use_extended_form::Bool,
    use_iterative_method::Bool,
    oa_solver,
)
    hypatia = MOI.OptimizerWithAttributes(
        Hypatia.Optimizer,
        MOI.Silent() => true,
        "near_factor" => 1000,
        "tol_feas" => 1e-10,
        "tol_rel_opt" => 1e-9,
        "tol_abs_opt" => 1e-8,
        "tol_illposed" => 1e-9,
        "tol_slow" => 2e-2,
        "tol_inconsistent" => 1e-7,
    )

    opt = JuMP.optimizer_with_attributes(
        MOIPajarito.Optimizer,
        # "verbose" => true,
        "oa_solver" => oa_solver,
        "conic_solver" => hypatia,
        "sep_solver" => hypatia,
        "use_extended_form" => use_extended_form,
        "use_iterative_method" => use_iterative_method,
        "debug_cuts" => use_iterative_method,
        "iteration_limit" => 30,
        # "time_limit" => 120.0,
    )

    @testset "$inst" for inst in insts
        @info inst
        eval(Symbol(inst))(opt)
    end
    return
end

# helpers

const rt2 = sqrt(2.0)

function svec(mat::AbstractMatrix{T}) where {T}
    vec = zeros(T, svec_length(size(mat, 1)))
    return Hypatia.Cones.smat_to_svec!(vec, mat, rt2)
end

function svec(mat::AbstractMatrix{Complex{T}}) where {T}
    vec = zeros(T, svec_length(ComplexF64, size(mat, 1)))
    return Hypatia.Cones.smat_to_svec!(vec, mat, rt2)
end

function smat(::Type{Float64}, vec::AbstractVector{T}) where {T}
    side = Hypatia.Cones.svec_side(length(vec))
    mat = zeros(T, side, side)
    return Hypatia.Cones.svec_to_smat!(mat, vec, rt2)
end

function smat(::Type{ComplexF64}, vec::AbstractVector{T}) where {T}
    side = Hypatia.Cones.svec_side(ComplexF64, length(vec))
    mat = zeros(Complex{T}, side, side)
    return Hypatia.Cones.svec_to_smat!(mat, vec, rt2)
end

function scale_svec(
    R::Type{<:Union{Float64, ComplexF64}},
    vec::AbstractVector{T},
    scal::Float64,
) where {T}
    incr = (R == Float64 ? 1 : 2)
    svec = copy(vec)
    k = 1
    for j in 1:Hypatia.Cones.svec_side(R, length(vec))
        for _ in 1:(incr * (j - 1))
            # scale off-diagonal
            svec[k] *= scal
            k += 1
        end
        k += 1
    end
    @assert k == length(vec) + 1
    return svec
end

include("JuMP_instances.jl")

end
