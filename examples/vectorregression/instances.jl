
insts = OrderedDict()

insts["test"] = [
    ((2, 3, false, VecSpecExt[]),),
    ((3, 3, true, VecSpecExt[]),),
    ((5, 5, false, [VecNegRtdet()]),),
    ((6, 4, true, [VecNegEntropy(), VecNegSqrtConj()]),),
]

# function vectorregression_insts(specs::Vector)
#     ts = Tuple[]
#     for (max_d, f) in specs
#         t = [((d, nothing, f),) for d in vcat(3, 3:3:max_d)]
#         append!(ts, t)
#     end
#     return ts
# end

# insts["nat"] = vectorregression_insts([(12, MatNegLog())])
# insts["ext"] = vectorregression_insts([(12, MatNegLogDirect())])

return (VectorRegression, insts)
