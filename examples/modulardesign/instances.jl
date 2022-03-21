
sparse_options = (; conic_solver = sparse_hypatia)
noext_options = (sparse_options..., use_extended_form = false)
nosubp_options = (sparse_options..., solve_subproblems = false) # for SOS2
noconic_options = (nosubp_options..., solve_relaxation = false)

insts = OrderedDict()

insts["test"] = [
    ((1, 1, 3, true),),
    ((1, 1, 3, false),),
    ((2, 3, 4, true),),
    ((2, 3, 4, false),),
    ((3, 2, 5, true),),
    ((3, 2, 5, false), noext_options),
    ((3, 2, 5, true), nosubp_options),
    ((3, 2, 5, false), noconic_options),
    ((1, 2, 3, true, true, SOS2(), 16), nosubp_options),
    ((2, 2, 4, true, true, SOS2(), 16), noconic_options),
    ((2, 2, 4, true, true, CCBounded(), 16), sparse_options),
    ((2, 2, 4, true, true, CCBounded(), 16), noconic_options),
    ((3, 3, 6, true, true, LogIBBounded(), 32), nosubp_options),
    ((3, 3, 6, true, true, LogIBBounded(), 32), noconic_options),
]

# convex instances
function modulardesign_insts(use_nat::Bool, max_n::Int, options::NamedTuple)
    return [((8, 4, n, use_nat), options) for n in vcat(20, 20:20:max_n)]
end

insts["nat"] = modulardesign_insts(true, 240, sparse_options)
insts["nat_noext"] = modulardesign_insts(true, 240, noext_options)
insts["ext"] = modulardesign_insts(false, 240, sparse_options)

# nonconvex instances
function modulardesign_insts(pwl::PWLSOS2, max_n::Int)
    return [((3, 3, n, true, true, pwl, 512), noconic_options) for n in vcat(4, 4:2:max_n)]
end

insts["sos2"] = modulardesign_insts(SOS2(), 30)
insts["logib"] = modulardesign_insts(LogIBBounded(), 30)
insts["cc"] = modulardesign_insts(CCBounded(), 18)

return (ModularDesign, insts)
