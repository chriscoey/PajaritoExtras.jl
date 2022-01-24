
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
return (TwoStageStochastic, insts)
