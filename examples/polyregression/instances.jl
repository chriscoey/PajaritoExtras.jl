
insts = OrderedDict()

insts["test"] = [
    ((1, 2, 15, 50.0, true),),
    ((1, 2, 15, 50.0, false),),
    ((3, 1, 40, 0.0, true),),
    ((3, 1, 40, 0.0, false),),
    ((2, 2, 40, 100.0, true),),
]

function polyregression_insts(use_nat::Bool)
    return [((2, 2, m, 25.0, use_nat),) for m in vcat(30, 30:5:70)]
end

insts["nat"] = polyregression_insts(true)
insts["ext"] = polyregression_insts(false)

return (PolyRegression, insts)
