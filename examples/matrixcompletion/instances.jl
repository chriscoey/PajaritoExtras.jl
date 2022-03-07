
insts = OrderedDict()

insts["test"] = [
    ((true, false, 3, 4, true),),
    ((true, false, 3, 4, false),),
    ((false, false, 3, 4, true),),
    ((false, false, 3, 4, false),),
    ((true, true, 6, 6, true),),
    ((true, true, 6, 6, false),),
    ((false, true, 6, 6, true),),
    ((false, true, 6, 6, false),),
]

# function matrixcompletion_insts(specs::Vector)
#     ts = Tuple[]
#     for (max_d, f) in specs
#         t = [((d, nothing, f),) for d in vcat(3, 3:3:max_d)]
#         append!(ts, t)
#     end
#     return ts
# end

# insts["nat"] = matrixcompletion_insts([(12, MatNegLog())])
# insts["ext"] = matrixcompletion_insts([(12, MatNegLogDirect())])

return (MatrixCompletion, insts)
