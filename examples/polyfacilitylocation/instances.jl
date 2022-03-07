
sparse_options = (; conic_solver = sparse_hypatia)

insts = OrderedDict()

insts["test"] = [
    ((1, 2, 2, true),),
    ((1, 2, 2, false),),
    ((2, 3, 3, true), sparse_options),
    ((2, 3, 3, false), sparse_options),
    ((6, 6, 3, true), sparse_options),
    ((6, 6, 3, false), sparse_options),
    ((20, 20, 3, true), sparse_options),
]

# function polyfacilitylocation_insts(specs::Vector)
#     ts = Tuple[]
#     for (max_d, f) in specs
#         t = [((d, nothing, f),) for d in vcat(3, 3:3:max_d)]
#         append!(ts, t)
#     end
#     return ts
# end

# insts["nat"] = polyfacilitylocation_insts([(12, MatNegLog())])
# insts["ext"] = polyfacilitylocation_insts([(12, MatNegLogDirect())])

return (PolyFacilityLocation, insts)
