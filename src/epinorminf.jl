#=
epigraph of ∞-norm of real or complex vectors
(u, w) : u ≥ ‖Wᵢ‖, ∀i
no extended formulation needed
real case (polyhedral): u ≥ Wᵢ, u ≥ -Wᵢ
complex case: (u, re(Wᵢ), im(Wᵢ)) in EpiNormEucl

dual cone is epigraph of 1-norm
(u, w) : u ≥ ∑ᵢ ‖Wᵢ‖
extended formulation: ∃ λ, u ≥ Σᵢ λᵢ
real case (polyhedral): λᵢ ≥ Wᵢ, λᵢ ≥ -Wᵢ
complex case: (λᵢ, re(Wᵢ), im(Wᵢ)) in EpiNormEucl
=#

mutable struct EpiNormInfCache{C <: RealOrComplex} <: ConeCache
    cone::Hypatia.EpiNormInfCone{Float64, C}
    oa_s::Vector{AE}
    s::Vector{Float64}
    d::Int
    w_temp::Vector{Float64}
    W_temp::Vector{C}
    EpiNormInfCache{C}() where {C <: RealOrComplex} = new{C}()
end

mutable struct DualEpiNormInfCache{C <: RealOrComplex, E <: Extender} <: ConeCache
    cone::Hypatia.EpiNormInfCone{Float64, C}
    oa_s::Vector{AE}
    s::Vector{Float64}
    d::Int
    λ::Vector{VR}
    w_temp::Vector{Float64}
    W_temp::Vector{C}
    DualEpiNormInfCache{C, E}() where {C <: RealOrComplex, E <: Extender} = new{C, E}()
end

function MOIPajarito.Cones.create_cache(
    oa_s::Vector{AE},
    cone::Hypatia.EpiNormInfCone{Float64, C},
    extend::Bool,
) where {C <: RealOrComplex}
    is_complex = (C == ComplexF64)
    dim = MOI.dimension(cone)
    @assert dim == length(oa_s)
    d = (is_complex ? div(dim - 1, 2) : dim - 1)
    cache = if cone.use_dual
        E = extender(extend, d)
        cache = DualEpiNormInfCache{C, E}()
    else
        EpiNormInfCache{C}()
    end
    cache.cone = cone
    # cache.is_complex = is_complex
    cache.oa_s = oa_s
    cache.d = d
    cache.w_temp = zeros(Float64, dim - 1)
    cache.W_temp = zeros(C, d)
    return cache
end

function MOIPajarito.Cones.add_init_cuts(
    cache::EpiNormInfCache{Float64},
    oa_model::JuMP.Model,
)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    # polyhedral (no OA)
    # JuMP.@constraints(oa_model, begin
    #     u .>= w
    #     u .>= -w
    # end)
    return 2 * cache.d
end

function MOIPajarito.Cones.add_init_cuts(
    cache::EpiNormInfCache{ComplexF64},
    oa_model::JuMP.Model,
)
    u = cache.oa_s[1]
    @views w = cache.oa_s[2:end]
    d = cache.d
    # cuts on (u, re(Wᵢ), im(Wᵢ)) are (rt2, ±1, ±1) ∈ EpiNormEucl
    # for i in 1:d
    #     re_w = w[2i - 1]
    #     im_w = w[2i]
    #     JuMP.@constraints(oa_model, begin
    #         rt2 * u - re_w - im_w >= 0
    #         rt2 * u - re_w + im_w >= 0
    #         rt2 * u + re_w - im_w >= 0
    #         rt2 * u + re_w + im_w >= 0
    #     end)
    # end
    return 4d
end

function MOIPajarito.Cones.get_subp_cuts(
    ::Vector{Float64},
    ::EpiNormInfCache{Float64},
    ::JuMP.Model,
)
    return AE[]
end

function MOIPajarito.Cones.get_subp_cuts(
    z::Vector{Float64},
    cache::EpiNormInfCache{ComplexF64},
    oa_model::JuMP.Model,
)
    # strengthened cut on (u, re(Wᵢ), im(Wᵢ)) is (‖Rᵢ‖, re(Rᵢ), im(Rᵢ))
    R = cache.W_temp
    @views vec_copyto!(R, z[2:end])
    cuts = AE[]
    for (i, R_i) in enumerate(R)
        abs_i = abs(R_i)
        abs_i < 1e-7 && continue
        cut = _get_cut(R_i, i, cache, oa_model)
        push!(cuts, cut)
    end
    return cuts
end

function MOIPajarito.Cones.get_sep_cuts(::EpiNormInfCache{Float64}, ::JuMP.Model)
    return AE[]
end

function MOIPajarito.Cones.get_sep_cuts(
    cache::EpiNormInfCache{ComplexF64},
    oa_model::JuMP.Model,
)
    us = cache.s[1]
    Ws = cache.W_temp
    @views vec_copyto!(Ws, cache.s[2:end])
    # check s ∉ K
    if us - norm(Ws, Inf) > -1e-7
        return AE[]
    end

    # gradient cut is (1, (-Wsᵢ / ‖Wsᵢ‖)ᵢ)
    cuts = AE[]
    for (i, Ws_i) in enumerate(Ws)
        abs_i = abs(Ws_i)
        abs_i < 1e-7 && continue
        R_i = Ws_i / -abs(Ws_i)
        cut = _get_cut(R_i, i, cache, oa_model)
        push!(cuts, cut)
    end
    return cuts
end

function _get_cut(
    R_i::ComplexF64,
    i::Int,
    cache::EpiNormInfCache{ComplexF64},
    oa_model::JuMP.Model,
)
    u = cache.oa_s[1]
    re_w = cache.oa_s[2i]
    im_w = cache.oa_s[2i + 1]
    return JuMP.@expression(oa_model, u + real(R_i) * re_w + imag(R_i) * im_w)
end

# function MOIPajarito.Cones.add_init_cuts(
#     cache::DualEpiNormInfCache{Float64, Unextended},
#     oa_model::JuMP.Model,
# )
#     u = cache.oa_s[1]
#     @views w = cache.oa_s[2:end]
#     # polyhedral but needs exponentially many constraints
#     JuMP.@constraints(oa_model, begin
#         # u .>= w
#         # u .>= -w
#     end)
#     # return 2d
# end

# function MOIPajarito.Cones.add_init_cuts(
#     cache::DualEpiNormInfCache{ComplexF64, Unextended},
#     oa_model::JuMP.Model,
# )
#     u = cache.oa_s[1]
#     @views w = cache.oa_s[2:end]
#     d = cache.d
#     # cuts on (u, re(Wᵢ), im(Wᵢ)) are (rt2, ±1, ±1) ∈ EpiNormEucl
#     for i in 1:d
#         re_w = w[2i - 1]
#         im_w = w[2i]
#         JuMP.@constraints(oa_model, begin
#             # rt2 * u - re_w - im_w >= 0
#             # rt2 * u - re_w + im_w >= 0
#             # rt2 * u + re_w - im_w >= 0
#             # rt2 * u + re_w + im_w >= 0
#         end)
#     end
#     return 4d
# end

# mutable struct EpiNormInfCache{E <: Extender} <: ConeCache
#     cone::Hypatia.EpiNormInfCone{Float64}
#     oa_s::Vector{AE}
#     s::Vector{Float64}
#     d::Int
#     λ::Vector{VR}
#     EpiNormInfCache{E}() where {E <: Extender} = new{E}()
# end

# function MOIPajarito.Cones.create_cache(
#     oa_s::Vector{AE},
#     cone::Hypatia.EpiNormInfCone{Float64},
#     extend::Bool,
# )
#     @assert !cone.use_dual # TODO
#     dim = MOI.dimension(cone)
#     @assert dim == length(oa_s)
#     d = dim - 1
#     E = extender(extend, d)
#     cache = EpiNormInfCache{E}()
#     cache.cone = cone
#     cache.oa_s = oa_s
#     cache.d = d
#     return cache
# end

# function MOIPajarito.Cones.get_subp_cuts(
#     z::Vector{Float64},
#     cache::EpiNormInfCache,
#     oa_model::JuMP.Model,
# )
#     return _get_cuts(z[2:end], cache, oa_model)
# end

# function MOIPajarito.Cones.get_sep_cuts(cache::EpiNormInfCache, oa_model::JuMP.Model)
#     s = cache.s
#     us = s[1]
#     @views ws = s[2:end]
#     @assert all(>(-1e-7), ws)
#     # check s ∉ K
#     if us < 1e-7 || geomean(ws) - us > -1e-7
#         return AE[]
#     end

#     # gradient cut is (-1, geom(w) / d ./ w)
#     w_pos = max.(ws, 1e-7)
#     c1 = geomean(w_pos) / cache.d
#     r = c1 ./ w_pos
#     return _get_cuts(r, cache, oa_model)
# end

# # unextended formulation

# function MOIPajarito.Cones.add_init_cuts(
#     cache::EpiNormInfCache{Unextended},
#     oa_model::JuMP.Model,
# )
#     u = cache.oa_s[1]
#     @views w = cache.oa_s[2:end] # TODO cache?
#     d = cache.d
#     # variable bounds wᵢ ≥ 0 and cut (-d, e)
#     JuMP.@constraints(oa_model, begin
#         [i in 1:d], w[i] >= 0
#         -d * u + sum(w) >= 0
#     end)
#     return d + 1
# end

# function _get_cuts(
#     r::Vector{Float64},
#     cache::EpiNormInfCache{Unextended},
#     oa_model::JuMP.Model,
# )
#     # strengthened cut is (-d * geom(r), r)
#     clean_array!(r) && return AE[]
#     p = -cache.d * geomean(r)
#     u = cache.oa_s[1]
#     @views w = cache.oa_s[2:end]
#     cut = JuMP.@expression(oa_model, p * u + JuMP.dot(r, w))
#     return [cut]
# end

# # extended formulation

# MOIPajarito.Cones.num_ext_variables(cache::EpiNormInfCache{Extended}) = 1 + cache.d

# function MOIPajarito.Cones.extend_start(
#     cache::EpiNormInfCache{Extended},
#     s_start::Vector{Float64},
# )
#     u_start = s_start[1]
#     w_start = s_start[2:end]
#     w_geom = geomean(w_start)
#     @assert w_geom - u_start >= -1e-7 # TODO
#     if u_start < 1e-8
#         return zeros(1 + cache.d)
#     end
#     λ_start = [u_start * log(w_i / u_start) for w_i in w_start]
#     @assert sum(λ_start) >= -eps()
#     return vcat(u_start, λ_start)
# end

# function MOIPajarito.Cones.setup_auxiliary(
#     cache::EpiNormInfCache{Extended},
#     oa_model::JuMP.Model,
# )
#     @assert cache.d >= 2
#     θ = cache.θ = JuMP.@variable(oa_model, lower_bound = 0)
#     u = cache.oa_s[1]
#     JuMP.@constraint(oa_model, θ >= u)
#     λ = cache.λ = JuMP.@variable(oa_model, [1:(cache.d)])
#     JuMP.@constraint(oa_model, sum(λ) >= 0)
#     return vcat(θ, λ)
# end

# function MOIPajarito.Cones.add_init_cuts(
#     cache::EpiNormInfCache{Extended},
#     oa_model::JuMP.Model,
# )
#     @views w = cache.oa_s[2:end]
#     d = cache.d
#     θ = cache.θ
#     λ = cache.λ
#     # variable bounds wᵢ ≥ 0 and cut (-d, e)
#     # disaggregated cut on (λᵢ, θ, wᵢ) is (-1, -1, 1)
#     # TODO check implication in math
#     JuMP.@constraints(oa_model, begin
#         [i in 1:d], w[i] >= 0
#         [i in 1:d], -λ[i] - θ + w[i] >= 0
#     end)
#     return 2d
# end

# function _get_cuts(
#     r::Vector{Float64},
#     cache::EpiNormInfCache{Extended},
#     oa_model::JuMP.Model,
# )
#     clean_array!(r) && return AE[]
#     p = -geomean(r)
#     @views w = cache.oa_s[2:end]
#     θ = cache.θ
#     λ = cache.λ
#     cuts = AE[]
#     for i in 1:(cache.d)
#         r_i = r[i]
#         iszero(r_i) && continue # TODO ?
#         # strengthened disaggregated cut on (λᵢ, θ, wᵢ) is
#         # p = -geom(r), (p, p * (log(-rᵢ / p) + 1), rᵢ)
#         # TODO check math
#         q = p * (log(-r_i / p) + 1)
#         cut = JuMP.@expression(oa_model, p * λ[i] + q * θ + r_i * w[i])
#         push!(cuts, cut)
#     end
#     return cuts
# end
