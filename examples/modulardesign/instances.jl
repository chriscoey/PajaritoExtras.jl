
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
    ((2, 2, 4, true, true, SOS2(), 16), noext_options),
    ((2, 2, 4, true, true, CCBounded(), 16), sparse_options),
    ((2, 2, 4, true, true, CCBounded(), 16), noconic_options),
    ((3, 3, 6, true, true, LogIBBounded(), 32), nosubp_options),
    ((3, 3, 6, true, true, LogIBBounded(), 32), noconic_options),
]

# convex instances
function modulardesign_insts(use_nat::Bool, options::NamedTuple)
    return [((10, 5, n, use_nat), options) for n in vcat(10, 10:10:100)]
end

insts["nat"] = modulardesign_insts(true, sparse_options)
insts["nat_noext"] = modulardesign_insts(true, noext_options)
insts["ext"] = modulardesign_insts(false, sparse_options)

# nonconvex instances
function modulardesign_insts(pwl::PWLSOS2, options::NamedTuple)
    return [((4, 4, n, true, true, pwl, 128), options) for n in vcat(4, 4:2:20)]
end

insts["sos2"] = modulardesign_insts(SOS2(), noconic_options)
insts["logib"] = modulardesign_insts(LogIBBounded(), noconic_options)
insts["cc"] = modulardesign_insts(CCBounded(), noconic_options)

return (ModularDesign, insts)
