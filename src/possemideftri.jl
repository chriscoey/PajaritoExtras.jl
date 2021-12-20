
mutable struct PosSemidefTriCache{C <: RealOrComplex} <: ConeCache
    cone::Hypatia.PosSemidefTriCone{Float64, C}
    is_complex::Bool
    oa_s::Vector{AE}
    s::Vector{Float64}
    d::Int
    W::Hermitian
    W_temp::Hermitian{C, Matrix{C}}
    PosSemidefTriCache{C}() where {C <: RealOrComplex} = new{C}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.PosSemidefTriCone{Float64, C},
    ::Bool,
) where {C <: RealOrComplex}
    cache = PosSemidefTriCache{C}()
    cache.cone = cone
    cache.is_complex = (C == ComplexF64)
    cache.oa_s = oa_s
    d = cache.d = svec_side(C, MOI.dimension(cone))
    AEC = (cache.is_complex ? Complex{AE} : AE)
    cache.W = Hermitian(zeros(AEC, d, d), :U)
    svec_to_smat!(cache.W.data, oa_s, rt2)
    cache.W_temp = Hermitian(zeros(C, d, d), :U)
    return cache
end

# initial OA polyhedron is the dual cone of diagonally dominant matrices
function MOIPajarito.Cones.add_init_cuts(
    cache::PosSemidefTriCache{C},
    oa_model::JuMP.Model,
) where {C}
    d = cache.d
    W = cache.W
    # W_ii ≥ 0
    JuMP.@constraint(oa_model, [i in 1:d], W[i, i] >= 0)
    # a linearization of W_ii * W_jj ≥ |W_ij|^2
    for j in 1:d, i in 1:(j - 1)
        _add_init_cuts(i, j, cache, oa_model)
    end
    return (cache.is_complex ? 2d^2 - d : d^2)
end

# real: cuts on (W_ii, W_jj, W_ij) are (1, 1, ±2), ∀i != j
function _add_init_cuts(
    i::Int,
    j::Int,
    cache::PosSemidefTriCache{Float64},
    oa_model::JuMP.Model,
)
    W = cache.W
    w_ii = W[i, i]
    w_jj = W[j, j]
    w_ij = W[i, j]
    JuMP.@constraints(oa_model, begin
        w_ii + w_jj + 2w_ij >= 0
        w_ii + w_jj - 2w_ij >= 0
    end)
    return
end

# complex: cuts on (W_ii, W_jj, Wr_ij, Wi_ij) are (1, 1, ±2, ±2), ∀i != j
function _add_init_cuts(
    i::Int,
    j::Int,
    cache::PosSemidefTriCache{ComplexF64},
    oa_model::JuMP.Model,
)
    W = cache.W
    w_ii = W[i, i]
    w_jj = W[j, j]
    (wr_ij, wi_ij) = reim(W[i, j])
    JuMP.@constraints(oa_model, begin
        w_ii + w_jj + 2wr_ij + 2wi_ij >= 0
        w_ii + w_jj + 2wr_ij - 2wi_ij >= 0
        w_ii + w_jj - 2wr_ij + 2wi_ij >= 0
        w_ii + w_jj - 2wr_ij - 2wi_ij >= 0
    end)
    return
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::PosSemidefTriCache,
    oa_model::JuMP.Model,
)
    # strengthened cuts from eigendecomposition are λᵢ * rᵢ * rᵢ'
    R = cache.W_temp
    svec_to_smat!(R.data, z, rt2)
    F = eigen!(R, 1e-10, Inf) # TODO tune
    isempty(F.values) && return AE[]
    R_eig = F.vectors * Diagonal(sqrt.(F.values))
    return _get_cuts(R_eig, cache, oa_model)
end

function MOIPajarito.Cones.get_sep_cuts(cache::PosSemidefTriCache, oa_model::JuMP.Model)
    # check s ∉ K
    sW = cache.W_temp
    svec_to_smat!(sW.data, cache.s, rt2)
    F = eigen!(sW, -Inf, -1e-7)
    isempty(F.values) && return AE[]
    return _get_cuts(F.vectors, cache, oa_model)
end

function _get_cuts(R_eig::Matrix{Float64}, cache::PosSemidefTriCache, oa_model::JuMP.Model)
    # cuts from eigendecomposition are rᵢ * rᵢ'
    cuts = AE[]
    for i in size(R_eig, 2)
        @views r_i = R_eig[:, i]
        R_i = r_i * r_i'
        clean_array!(R_i) && continue
        cut = dot_expr(Hermitian(R_i, :U), cache.W, oa_model)
        push!(cuts, cut)
    end
    return cuts
end
