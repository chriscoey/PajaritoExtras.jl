
sparse_options = (; conic_solver = sparse_hypatia)
nosubp_options = (; conic_solver = sparse_hypatia, solve_subproblems = false) # for SOS2

insts = OrderedDict()

insts["test"] = [
    ((1, 1, 3, true),),
    ((1, 1, 3, false),),
    ((2, 3, 4, true),),
    ((2, 3, 4, false),),
    ((3, 2, 5, true),),
    ((3, 2, 5, true), sparse_options),
    ((3, 2, 5, true), nosubp_options),
    ((1, 6, 10, true), (; solve_relaxation = false)),
    ((1, 3, 6, true), (; solve_relaxation = false)),
    ((1, 2, 3, true, true, SOS2(), 10), nosubp_options),
    ((1, 3, 6, true, true, CCBounded(), 10),),
    ((2, 2, 4, true, true, SOS2(), 20), nosubp_options),
    ((2, 2, 4, true, true, CCBounded(), 20), sparse_options),
    ((2, 2, 4, true, true, CCBounded(), 20), nosubp_options),
    ((3, 3, 6, true, true, LogIBBounded(), 20), sparse_options),
    ((3, 3, 6, true, true, LogIBBounded(), 20), nosubp_options),
    ((1, 6, 5, true, true, SOS2(), 20), nosubp_options),
]

# convex instances
function modulardesign_insts(use_nat::Bool)
    return [((5, 5, n, use_nat), sparse_options) for n in vcat(20, 25:25:250)]
end

insts["nat"] = modulardesign_insts(true)
insts["ext"] = modulardesign_insts(false)

# nonconvex instances
function modulardesign_insts(pwl::PWLSOS2, options::NamedTuple)
    return [((3, 3, n, true, true, pwl, 128), options) for n in vcat(20, 20:20:160)]
end

insts["logib"] = modulardesign_insts(LogIBBounded(), nosubp_options)
insts["cc"] = modulardesign_insts(CCBounded(), nosubp_options)
insts["sos2"] = modulardesign_insts(SOS2(), nosubp_options)

return (ModularDesign, insts)
