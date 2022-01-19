
sparse_options = (; conic_solver = sparse_hypatia)

insts = OrderedDict()
insts["test"] = [
    ((1, 2, 2, true),),
    ((1, 2, 2, false),),
    ((2, 3, 3, true), sparse_options),
    ((2, 3, 3, false), sparse_options),
    ((6, 6, 3, true), sparse_options),
    ((6, 6, 3, false), sparse_options),
    ((20, 20, 3, true), sparse_options),
]
return (PolyFacilityLocation, insts)
