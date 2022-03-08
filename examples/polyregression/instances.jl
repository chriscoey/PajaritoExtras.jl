
insts = OrderedDict()

insts["test"] = [
    ((1, 2, 15, 50.0, true),),
    ((1, 2, 15, 50.0, false),),
    ((3, 1, 40, 0.0, true),),
    ((3, 1, 40, 0.0, false),),
    ((2, 2, 40, 100.0, true),),
]

function polyregression_insts(use_nat::Bool)
    return [((n, 2, 50, 10.0, use_nat), sparse_options) for n in vcat(1, 1:6)]
end

insts["nat"] = polyregression_insts(true)
insts["ext"] = polyregression_insts(false)

return (PolyRegression, insts)
