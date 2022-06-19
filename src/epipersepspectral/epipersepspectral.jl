#=
epigraph of perspective of a separable spectral function
on a real vector or real symmetric or complex Hermitian (svec scaled triangle) domain
(u, v, w) : v ≥ 0, W ⪰ 0, u ≥ v * h(W / v)

dual cone is
(u, v, w) : u ≥ 0, W ∈ dom(h⋆), v ≥ u * h⋆(W / u)
=#

import Hypatia.Cones: VectorCSqr, MatrixCSqr, SepSpectralFun
import Hypatia.Cones: h_val, h_conj_dom_pos, h_conj, h_der1

include("vector.jl")
include("matrix.jl")

function Pajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiPerSepSpectralCone{RealF},
    opt::Optimizer,
)
    @assert MOI.dimension(cone) == length(oa_s)
    cache = create_sepspectral_cache(cone.Q, cone.use_dual, cone.d, opt)
    cache.oa_s = oa_s
    cache.d = cone.d
    cache.h = cone.h
    return cache
end

function per_sepspec(f::Function, h::SepSpectralFun, v::RealF, w::AbstractVector{RealF})
    v < 1e-12 && return 0.0 # TODO?
    if v == 1
        return f(w, h)
    end
    return v * f(w / v, h)
end

swap_epiper(::Type{Prim}, u::T, v::T) where {T} = (u, v)
swap_epiper(::Type{Dual}, u::T, v::T) where {T} = (v, u)

dom_pos(::Type{Prim}, ::SepSpectralFun) = true
dom_pos(::Type{Dual}, h::SepSpectralFun) = h_conj_dom_pos(h)

val_or_conj(::Type{Prim}) = h_val
val_or_conj(::Type{Dual}) = h_conj
conj_or_val(::Type{Prim}) = h_conj
conj_or_val(::Type{Dual}) = h_val

import Hypatia.Cones: NegLogSSF, NegEntropySSF, NegSqrtSSF, NegPower01SSF, Power12SSF

# r values for initial cuts
init_r_vals(::Type{<:Union{Prim, Dual}}, ::NegLogSSF) = [1e-4, 1.0, 1e4]
init_r_vals(::Type{Prim}, ::NegEntropySSF) = [-9.0, -1.0, 8.0]
init_r_vals(::Type{Dual}, ::NegEntropySSF) = [1e-4, 1.0, 1e3]
init_r_vals(::Type{<:Union{Prim, Dual}}, ::NegSqrtSSF) = [1e-4, 1.0, 1e4]

function init_r_vals(::Type{Prim}, h::NegPower01SSF)
    vals = [1e-4, 1.0, 1e4]
    p = RealF(h.p)
    p <= 0.5 && return vals
    q = p / (p - 1)
    return [p * (v / (1 - p))^inv(q) for v in vals]
end

init_r_vals(::Type{Dual}, h::NegPower01SSF) = [1e-4, 1.0, 1e4]

function init_r_vals(::Type{Prim}, h::Power12SSF)
    p = RealF(h.p)
    p >= 1.5 && return [-1e2, -3.0, 1.0]
    q = p / (p - 1)
    vals = [-1e4, -1.0, 1.0]
    return [sign(v) * p * (abs(v) / (p - 1))^inv(q) for v in vals]
end

init_r_vals(::Type{Dual}, ::Power12SSF) = [1e-3, 1.0, 1e3]

# gradients
h_grad(::Type{Prim}, h::SepSpectralFun, w::Vector{RealF}) = h_der1(similar(w), w, h)

h_grad(::Type{Dual}, h::NegLogSSF, w::Vector{RealF}) = h_der1(similar(w), w, h)
h_grad(::Type{Dual}, ::NegEntropySSF, w::Vector{RealF}) = [-exp(-w_i - 1) for w_i in w]
h_grad(::Type{Dual}, ::NegSqrtSSF, w::Vector{RealF}) = [-0.25 * w_i^-2 for w_i in w]

function h_grad(::Type{Dual}, h::NegPower01SSF, w::Vector{RealF})
    p = RealF(h.p)
    qm1 = inv(p - 1)
    c = -p^-qm1
    return [c * w_i^qm1 for w_i in w]
end

function h_grad(::Type{Dual}, h::Power12SSF, w::Vector{RealF})
    p = RealF(h.p)
    qm1 = inv(p - 1)
    c = -p^-qm1
    g = zero(w)
    for (i, w_i) in enumerate(w)
        if w_i < 0
            g[i] = c * abs(w_i)^qm1
        end
    end
    return g
end
