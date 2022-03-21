
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

function matrixdecomposition_insts(use_nat::Bool, d2_max::Int)
    return [((15, d2, use_nat),) for d2 in vcat(50, 50:50:d2_max)]
end

insts["nat"] = matrixdecomposition_insts(true, 700)
insts["ext"] = matrixdecomposition_insts(false, 400)

return (MatrixDecomposition, insts)
