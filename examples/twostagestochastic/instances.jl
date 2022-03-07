
sparse_options = (; conic_solver = sparse_hypatia)

insts = OrderedDict()

insts["test"] = [
    ((2, 2, false, true), sparse_options),
    ((2, 2, false, false), sparse_options),
    ((2, 2, true, true), sparse_options),
    ((2, 2, true, false), sparse_options),
    ((3, 3, false, true), sparse_options),
    ((3, 3, false, false), sparse_options),
    ((3, 3, true, true), sparse_options),
    ((10, 3, false, true), sparse_options),
    ((5, 2, true, true), sparse_options),
]

# function twostagestochastic_insts(specs::Vector)
#     ts = Tuple[]
#     for (max_d, f) in specs
#         t = [((d, nothing, f),) for d in vcat(3, 3:3:max_d)]
#         append!(ts, t)
#     end
#     return ts
# end

# insts["nat"] = twostagestochastic_insts([(12, MatNegLog())])
# insts["ext"] = twostagestochastic_insts([(12, MatNegLogDirect())])

return (TwoStageStochastic, insts)
