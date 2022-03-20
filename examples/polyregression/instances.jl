
sparse_options = (; conic_solver = sparse_hypatia)

insts = OrderedDict()

insts["test"] = [
    ((1, 2, 15, 50.0, true),),
    ((1, 2, 15, 50.0, false),),
    ((3, 1, 40, 0.0, true), sparse_options),
    ((3, 1, 40, 0.0, false), sparse_options),
    ((2, 2, 40, 100.0, true), sparse_options),
]

function polyregression_insts(use_nat::Bool)
    return [((2, 2, m, 25.0, use_nat), sparse_options) for m in vcat(30, 30:5:85)]
end

insts["nat"] = polyregression_insts(true)
insts["ext"] = polyregression_insts(false)

return (PolyRegression, insts)
