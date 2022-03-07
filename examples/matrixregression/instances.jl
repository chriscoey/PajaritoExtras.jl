
sparse_options = (; conic_solver = sparse_hypatia)

insts = OrderedDict()

insts["test"] = [
    ((5, 3, 3, true),),
    ((5, 3, 3, false),),
    ((15, 7, 5, true), sparse_options),
    ((15, 7, 5, false), sparse_options),
    ((50, 15, 10, true), sparse_options),
    ((50, 15, 10, false), sparse_options),
]

# function matrixregression_insts(specs::Vector)
#     ts = Tuple[]
#     for (max_d, f) in specs
#         t = [((d, nothing, f),) for d in vcat(3, 3:3:max_d)]
#         append!(ts, t)
#     end
#     return ts
# end

# insts["nat"] = matrixregression_insts([(12, MatNegLog())])
# insts["ext"] = matrixregression_insts([(12, MatNegLogDirect())])

return (MatrixRegression, insts)
