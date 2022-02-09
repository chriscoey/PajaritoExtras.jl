#=
JuMP helpers for constructing formulations for spectral/nuclear norm cone constraints
=#

function add_symm_specnuc(
    use_nat::Bool,
    is_nuc::Bool,
    u_aff::JuMPScalar,
    W_aff::AbstractMatrix{<:JuMPScalar},
    model::JuMP.Model,
)
    nrow = size(W_aff, 1)
    @assert issymmetric(W_aff)
    if use_nat
        len = Cones.svec_length(nrow)
        K = Hypatia.EpiNormSpectralTriCone{Float64, Float64}(1 + len, is_nuc)
        W_vec = Vector{JuMP.AffExpr}(undef, len)
        Cones.smat_to_svec!(W_vec, W_aff, rt2)
        JuMP.@constraint(model, vcat(u_aff, W_vec) in K)
    else
        if is_nuc
            # EF for symmetric nuclear norm (L1 of eigvals)
            JuMP.@variable(model, X1[1:nrow, 1:nrow], PSD)
            JuMP.@variable(model, X2[1:nrow, 1:nrow], PSD)
            eq = [W_aff[i, j] - X1[i, j] + X2[i, j] for j in 1:nrow for i in 1:j]
            JuMP.@constraint(model, eq .== 0)
            JuMP.@constraint(model, u_aff >= tr(X1) + tr(X2))
        else
            # EF for symmetric spectral norm (Linf of eigvals)
            uI = u_aff * Matrix(I, nrow, nrow)
            JuMP.@constraint(model, Symmetric(uI - W_aff) in JuMP.PSDCone())
            JuMP.@constraint(model, Symmetric(uI + W_aff) in JuMP.PSDCone())
        end
    end
    return
end

function add_nonsymm_specnuc(
    use_nat::Bool,
    is_nuc::Bool,
    u_aff::JuMPScalar,
    W_aff::AbstractMatrix{<:JuMPScalar},
    model::JuMP.Model,
)
    (nrow, ncol) = size(W_aff)
    if use_nat
        K = Hypatia.EpiNormSpectralCone{Float64, Float64}(nrow, ncol, is_nuc)
        JuMP.@constraint(model, vcat(u_aff, vec(W_aff)) in K)
    else
        if is_nuc
            # EF for nuclear norm
            X1 = JuMP.@variable(model, [1:nrow, 1:nrow], Symmetric)
            X2 = JuMP.@variable(model, [1:ncol, 1:ncol], Symmetric)
            JuMP.@constraint(model, u_aff >= (tr(X1) + tr(X2)) / 2)
        else
            # EF for spectral norm
            X1 = u_aff * Matrix(I, nrow, nrow)
            X2 = u_aff * Matrix(I, ncol, ncol)
        end
        mat = hvcat((2, 2), X1, W_aff, W_aff', X2)
        JuMP.@constraint(model, Symmetric(mat) in JuMP.PSDCone())
    end
    return
end
