# JuMP instance tests

module TestJuMP

using Test
using LinearAlgebra
import MathOptInterface
const MOI = MathOptInterface
import JuMP
import MOIPajarito
import PajaritoExtras
import PajaritoExtras: svec_idx
import Hypatia
import Hypatia.Cones: svec_length

# all instances
inst_all = [
    # "possemideftri1",
    # "possemideftri2",
    # "epinormeucl1",
    # "epinormeucl2",
    # "hypogeomean1",
    # "hypogeomean2",
    "epipersquare1",
    "epipersquare2",
    # "expdesign",
]

# instances to test with use_extended_form = false
inst_noextend = [
    # "epinormeucl1",
    # "epinormeucl2",
    # "hypogeomean1",
    # "hypogeomean2",
    "epipersquare1",
    "epipersquare2",
]

function runtests(oa_solver, conic_solver)
    @testset "iterative, natural" begin
        @info "iterative, natural"
        run_jump_tests(inst_noextend, false, true, oa_solver, conic_solver)
    end
    @testset "iterative, extended" begin
        @info "iterative, extended"
        run_jump_tests(inst_all, true, true, oa_solver, conic_solver)
    end
    @testset "one tree, natural" begin
        @info "one tree, natural"
        run_jump_tests(inst_noextend, false, false, oa_solver, conic_solver)
    end
    @testset "one tree, extended" begin
        @info "one tree, extended"
        run_jump_tests(inst_all, true, false, oa_solver, conic_solver)
    end
    return
end

function run_jump_tests(
    insts::Vector{String},
    use_extended_form::Bool,
    use_iterative_method::Bool,
    oa_solver,
    conic_solver,
)
    opt = JuMP.optimizer_with_attributes(
        MOIPajarito.Optimizer,
        # "verbose" => true,
        # "verbose" => false,
        "oa_solver" => oa_solver,
        "conic_solver" => conic_solver,
        "use_extended_form" => use_extended_form,
        "use_iterative_method" => use_iterative_method,
        # "debug_cuts" => use_iterative_method,
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

include("JuMP_instances.jl")

end
