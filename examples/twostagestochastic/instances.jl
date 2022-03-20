
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
    p_max = (use_nat ? 11 : 7) # ext runs out of memory early
    return [((3, 2^p, false, use_nat), sparse_options) for p in vcat(3, 3:p_max)]
end

insts["nat"] = twostagestochastic_insts(true)
insts["ext"] = twostagestochastic_insts(false)

return (TwoStageStochastic, insts)
