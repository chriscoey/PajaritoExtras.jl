
sparse_options = (; conic_solver = sparse_hypatia)

insts = OrderedDict()

insts["test"] = [
    ((2, 2, false, true), sparse_options),
    ((2, 2, false, false), sparse_options),
    ((2, 2, true, true), sparse_options),
    ((2, 2, true, false), sparse_options),
    ((3, 3, false, true), sparse_options),
    ((3, 3, false, false), sparse_options),
    ((3, 3, true, true), sparse_options),
    ((10, 3, false, true), sparse_options),
    ((5, 2, true, true), sparse_options),
]

function twostagestochastic_insts(use_nat::Bool)
    return [((3, d, false, use_nat), sparse_options) for d in [3, 10, 30, 100, 300, 1000]]
end

insts["nat"] = twostagestochastic_insts(true)
insts["ext"] = twostagestochastic_insts(false)

return (TwoStageStochastic, insts)
