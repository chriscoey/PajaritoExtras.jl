
sparse_options = (; conic_solver = sparse_hypatia)

insts = OrderedDict()
insts["minimal"] = [((1, 2, 2),), ((2, 3, 3),), ((20, 20, 3), sparse_options)]
return (PolyFacilityLocation, insts)
