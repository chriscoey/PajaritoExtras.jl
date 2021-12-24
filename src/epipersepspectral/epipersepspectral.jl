#=
epigraph of perspective of a separable spectral function
on a real vector or real symmetric or complex Hermitian (svec scaled triangle) domain
(u, v, w) : v ≥ 0, W ⪰ 0, u ≥ v * h(W / v)

dual cone is
(u, v, w) : u ≥ 0, W ∈ dom(h⋆), v ≥ u * h⋆(W / u)
=#

import Hypatia.Cones: VectorCSqr, SepSpectralFun, h_val, h_conj_dom_pos, h_conj

include("vector.jl")
# include("matrix.jl")

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiPerSepSpectralCone{Float64},
    extend::Bool,
)
    @assert !cone.use_dual # TODO
    dim = MOI.dimension(cone)
    @assert dim == length(oa_s)
    cache = create_sepspectral_cache(cone.Q, cone.d, extend)
    cache.cone = cone
    cache.oa_s = oa_s
    cache.d = cone.d
    cache.h = cone.h
    return cache
end

function per_sepspec(f::Function, h::SepSpectralFun, v::Float64, w::AbstractVector{Float64})
    v < 1e-12 && return 0.0
    if v == 1
        return f(w, h)
    end
    return v * f(w / v, h)
end

import Hypatia.Cones: NegLogSSF, NegEntropySSF, NegSqrtSSF, NegPower01SSF, Power12SSF

# r values for initial cuts
init_r_vals(::NegLogSSF) = [1e-4, 1.0, 1e4]
init_r_vals(::NegEntropySSF) = [-9.0, -1.0, 8.0]
init_r_vals(::NegSqrtSSF) = [1e-3, 1.0, 1e3]
init_r_vals(::NegPower01SSF) = [1e-3, 1.0, 1e3]
init_r_vals(::Power12SSF) = [-10.0, -1.0, 1.0]
