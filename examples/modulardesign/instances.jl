
sparse_options = (; conic_solver = sparse_hypatia)
nosubp_options = (; conic_solver = sparse_hypatia, solve_subproblems = false) # for SOS2

insts = OrderedDict()

insts["test"] = [
    ((1, 1, 3),),
    ((2, 3, 4),),
    ((3, 2, 5),),
    ((3, 2, 5), sparse_options),
    ((3, 2, 5), nosubp_options),
    ((1, 6, 10), (; solve_relaxation = false)),
    ((1, 3, 6), (; solve_relaxation = false)),
    ((1, 2, 3, true, false, SOS2(), 10), nosubp_options),
    ((1, 2, 3, true, true, SOS2(), 10), nosubp_options),
    ((1, 3, 6, true, false, CCBounded(), 10),),
    ((1, 3, 6, true, true, CCBounded(), 40),),
    ((2, 2, 4, true, false, SOS2(), 20), nosubp_options),
    ((2, 2, 4, true, false, CCBounded(), 20), sparse_options),
    ((2, 2, 4, true, false, CCBounded(), 20), nosubp_options),
    ((3, 3, 6, true, true, LogIBBounded(), 20), sparse_options),
    ((3, 3, 6, true, true, LogIBBounded(), 20), nosubp_options),
    ((1, 6, 5, true, false, SOS2(), 20), nosubp_options),
]

# function modulardesign_insts(specs::Vector)
#     ts = Tuple[]
#     for (max_d, f) in specs
#         t = [((d, nothing, f),) for d in vcat(3, 3:3:max_d)]
#         append!(ts, t)
#     end
#     return ts
# end

# insts["nat"] = modulardesign_insts([(12, MatNegLog())])
# insts["ext"] = modulardesign_insts([(12, MatNegLogDirect())])

return (ModularDesign, insts)
