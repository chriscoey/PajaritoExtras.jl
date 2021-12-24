#=
epigraph of perspective of a separable spectral function
on a real vector or real symmetric or complex Hermitian (svec scaled triangle) domain
(u, v, w) : v ≥ 0, W ⪰ 0, u ≥ v * h(W / v)

dual cone is
(u, v, w) : u ≥ 0, W ∈ dom(h⋆), v ≥ u * h⋆(W / u)
=#

import Hypatia.Cones: SepSpectralFun, h_val, h_conj_dom_pos, h_conj

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
