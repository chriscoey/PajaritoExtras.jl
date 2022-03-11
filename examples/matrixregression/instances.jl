
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

function matrixregression_insts(use_nat::Bool)
    return [((2 * p, p, 5, use_nat), sparse_options) for p in vcat(10, 10:10:100)]
end

insts["nat"] = matrixregression_insts(true)
insts["ext"] = matrixregression_insts(false)

return (MatrixRegression, insts)
