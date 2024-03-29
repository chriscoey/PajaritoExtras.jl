#=
interpolant-basis weighted sum-of-squares polynomials (real or real-valued complex)
parametrized by vector of matrices `Ps` derived from interpolant basis and domain constraints

primal cone
w : ∃ Sₗ ⪰ 0 ∀ l, w = ∑ₗ diag(Pₗ Sₗ Pₗ')

dual cone
w : Pₗ' Diagonal(w) Pₗ ⪰ 0 ∀ l
=#

mutable struct WSOSInterpNonnegative{D <: PrimDual, C <: RealCompF} <: Cache
    oa_s::Vector{AE}
    Ps::Vector{Matrix{C}}
    d::Int
    sep_constr::CR
    WSOSInterpNonnegative{D, C}() where {D <: PrimDual, C <: RealCompF} = new{D, C}()
end

function Pajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.WSOSInterpNonnegativeCone{RealF, C},
    opt::Optimizer,
) where {C <: RealCompF}
    D = primal_or_dual(cone.use_dual)
    cache = WSOSInterpNonnegative{D, C}()
    cache.oa_s = oa_s
    cache.Ps = cone.Ps
    d = cache.d = MOI.dimension(cone)
    @assert length(oa_s) == d
    @assert size(cache.Ps[1], 1) == d
    if D == Prim
        cache.sep_constr = get_sep_constr(cone, opt)
    end
    return cache
end

# TODO maybe define Hypatia MOI cone equal/hash in Hypatia
function hash_cone(cone::Hypatia.WSOSInterpNonnegativeCone{RealF})
    @assert !cone.use_dual
    return hash(cone.Ps)
end

# primal cone functions

function Pajarito.Cones.add_init_cuts(cache::WSOSInterpNonnegative{Prim}, opt::Optimizer)
    # add variable bounds wᵢ ≥ 0
    JuMP.@constraint(opt.oa_model, cache.oa_s .>= 0)
    return
end

function Pajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::WSOSInterpNonnegative{Prim},
    opt::Optimizer,
)
    cut = dot_expr(z, cache.oa_s, opt)
    return [cut]
end

# dual cone functions

function Pajarito.Cones.add_init_cuts(cache::WSOSInterpNonnegative{Dual}, opt::Optimizer)
    opt.use_init_fixed_oa || return

    # add cuts to make diagonal of each Pₗ' Diagonal(w) Pₗ nonnegative
    # TODO could add SOC cuts or linearizations for 2x2 principal matrix PSD condition
    w = cache.oa_s
    r = zeros(cache.d)
    for P in cache.Ps, P_i in eachcol(P)
        @. r = abs2(P_i)
        JuMP.@constraint(opt.oa_model, dot(r, w) >= 0)
    end
    return
end

function Pajarito.Cones.get_subp_cuts(
    z::Vector{RealF},
    cache::WSOSInterpNonnegative{Dual},
    opt::Optimizer,
)
    # TODO to get an extreme ray decomposition, use the technique from PY
    # that uses the primal-dual solution from Hypatia to get values for each Sₗ ⪰ 0
    cut = dot_expr(z, cache.oa_s, opt)
    return [cut]
end

function Pajarito.Cones.get_sep_cuts(
    s::Vector{RealF},
    cache::WSOSInterpNonnegative{Dual},
    opt::Optimizer,
)
    # extreme ray decomposition adds eigenvector cuts to make
    # each Pₗ' Diagonal(w) Pₗ PSD
    r = zeros(cache.d)
    cuts = AE[]
    for P in cache.Ps
        Λ = P' * Diagonal(s) * P
        F = eigen!(Hermitian(Λ, :U), -Inf, -opt.tol_feas)
        isempty(F.values) && continue

        PV = P * F.vectors
        for PV_i in eachcol(PV)
            @. r = abs2(PV_i)
            cut = dot_expr(r, cache.oa_s, opt)
            push!(cuts, cut)
        end
    end
    return cuts
end
