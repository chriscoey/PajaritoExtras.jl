
insts = OrderedDict()

insts["test"] = [
    ((3, 4, 0.3, 2, true),)
    ((3, 4, 0.3, 2, false),)
    ((5, 8, 0.2, 4, true),)
    ((5, 8, 0.2, 4, false),)
]

# function matrixdecomposition_insts(specs::Vector)
#     ts = Tuple[]
#     for (max_d, f) in specs
#         t = [((d, nothing, f),) for d in vcat(3, 3:3:max_d)]
#         append!(ts, t)
#     end
#     return ts
# end

# insts["nat"] = matrixdecomposition_insts([(12, MatNegLog())])
# insts["ext"] = matrixdecomposition_insts([(12, MatNegLogDirect())])

return (MatrixDecomposition, insts)
