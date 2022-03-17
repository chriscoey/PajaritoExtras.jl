
noconic_options = (; solve_subproblems = false, solve_relaxation = false)

insts = OrderedDict()

insts["test"] = [
    ((3, 4, true),)
    ((3, 4, false),)
    ((3, 4, true), noconic_options)
    ((3, 4, false), noconic_options)
    ((5, 8, true),)
    ((5, 8, false),)
]

function matrixdecomposition_insts(use_nat::Bool)
    return [((10, d2, use_nat),) for d2 in vcat(30, 30:30:360)]
end

insts["nat"] = matrixdecomposition_insts(true)
insts["ext"] = matrixdecomposition_insts(false)

return (MatrixDecomposition, insts)
