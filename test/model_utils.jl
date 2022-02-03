# tools for model building

import Hypatia
import Hypatia.Cones: smat_to_svec!, svec_to_smat!, svec_side, svec_length

const rt2 = sqrt(2)

function svec(mat::AbstractMatrix{T}) where {T}
    vec = zeros(T, svec_length(size(mat, 1)))
    return smat_to_svec!(vec, mat, rt2)
end

function svec(mat::AbstractMatrix{Complex{T}}) where {T}
    vec = zeros(T, svec_length(ComplexF64, size(mat, 1)))
    return smat_to_svec!(vec, mat, rt2)
end

function smat(::Type{Float64}, vec::AbstractVector{T}) where {T}
    side = svec_side(length(vec))
    mat = zeros(T, side, side)
    return svec_to_smat!(mat, vec, rt2)
end

function smat(::Type{ComplexF64}, vec::AbstractVector{T}) where {T}
    side = svec_side(ComplexF64, length(vec))
    mat = zeros(Complex{T}, side, side)
    return svec_to_smat!(mat, vec, rt2)
end

function scale_svec(
    R::Type{<:Union{Float64, ComplexF64}},
    vec::AbstractVector{T},
    scal::Float64,
) where {T}
    incr = (R == Float64 ? 1 : 2)
    svec = copy(vec)
    k = 1
    for j in 1:svec_side(R, length(vec))
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
