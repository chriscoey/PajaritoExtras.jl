
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

return (MatrixRegression, insts)
