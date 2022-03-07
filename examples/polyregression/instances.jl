
insts = OrderedDict()

insts["test"] = [
    ((1, 2, 15, 50.0, true),),
    ((1, 2, 15, 50.0, false),),
    ((3, 1, 40, 0.0, true),),
    ((3, 1, 40, 0.0, false),),
    ((2, 2, 40, 100.0, true),),
]

# function polyregression_insts(specs::Vector)
#     ts = Tuple[]
#     for (max_d, f) in specs
#         t = [((d, nothing, f),) for d in vcat(3, 3:3:max_d)]
#         append!(ts, t)
#     end
#     return ts
# end

# insts["nat"] = polyregression_insts([(12, MatNegLog())])
# insts["ext"] = polyregression_insts([(12, MatNegLogDirect())])

return (PolyRegression, insts)
