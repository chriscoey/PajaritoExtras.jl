#=
hypograph of root-determinant of real symmetric or complex Hermitian
positive semidefinite matrices (svec scaled triangle)
(u, w) : W ⪰ 0, u ≤ rtdet(W)

dual cone is
(u, w) : u ≤ 0, W ⪰ 0, u ≥ -d * rtdet(W)
=#

mutable struct HypoRootdetTriCache{C <: RealOrComplex} <: ConeCache
    cone::Hypatia.HypoRootdetTriCone{Float64, C}
    is_complex::Bool
    oa_s::Vector{AE}
    s::Vector{Float64}
    d::Int
    w_temp::Vector{Float64}
    W_temp::Matrix{C}
    HypoRootdetTriCache{C}() where {C <: RealOrComplex} = new{C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.HypoRootdetTriCone{Float64, C},
    ::Bool,
) where {C <: RealOrComplex}
    cache = HypoRootdetTriCache{C}()
    cache.cone = cone
    cache.is_complex = (C == ComplexF64)
    cache.oa_s = oa_s
    dim = MOI.dimension(cone)
    d = cache.d = svec_side(C, dim - 1)
    cache.w_temp = zeros(Float64, dim - 1)
    cache.W_temp = zeros(C, d, d)
    return cache
end

function MOIPajarito.Cones.add_init_cuts(
    cache::HypoRootdetTriCache{C},
    oa_model::JuMP.Model,
) where {C}
    d = cache.d
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    # W_ii ≥ 0
    # linearize at W = I, cut on (u, W) is (-d, I)
    W_diag = [w[svec_idx(C, i, i)] for i in 1:d]
    JuMP.@constraints(oa_model, begin
        W_diag .>= 0
        -d * u + sum(W_diag) >= 0
    end)
    return d + 1
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::HypoRootdetTriCache,
    oa_model::JuMP.Model,
)
    R = cache.W_temp
    @views svec_to_smat!(R, z[2:end], rt2)
    F = eigen!(Hermitian(R, :U))
    cuts = _get_cuts(F.values, F.vectors, cache, oa_model)
    # return AE[]

    # if z[1] > -1e-7
    #     # add eigenvector cuts
    #     append!(cuts, ...)
    # end
    return cuts
end

function MOIPajarito.Cones.get_sep_cuts(cache::HypoRootdetTriCache, oa_model::JuMP.Model)
    # check s ∉ K
    us = cache.s[1]
    Ws = cache.W_temp
    @views svec_to_smat!(Ws, cache.s[2:end], rt2)
    F = eigen!(Hermitian(Ws, :U))
    ω = F.values
    us < 1e-7 && all(>(-1e-7), ω) && return AE[]
    Ws_geom = geomean(ω)
    Ws_geom - us > -1e-7 && return AE[]

    # gradient cut is (-1, V * Diagonal(geom(ω) / d ./ ω) * V')
    ω_pos = max.(ω, 1e-7)
    rω = (geomean(ω_pos) / cache.d) ./ ω_pos
    cuts = _get_cuts(rω, F.vectors, cache, oa_model)

    # if isnan(Ws_geom)
    #     # add eigenvector cuts
    #     # TODO refac with other cones? need to pass in w variables
    #     append!(cuts, _get_psd_cuts(F.vectors, cache, oa_model))
    # end
    return cuts
end

function _get_cuts(
    rω::Vector{Float64},
    V::Matrix{C},
    cache::HypoRootdetTriCache{C},
    oa_model::JuMP.Model,
) where {C}
    # strengthened cut is (-d * geom(rω), V * Diagonal(rω) * V')
    p = -cache.d * geomean(rω)
    R = V * Diagonal(rω) * V'
    r = smat_to_svec!(cache.w_temp, R, rt2)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    cut = JuMP.@expression(oa_model, p * u + JuMP.dot(r, w))
    return [cut]
end
