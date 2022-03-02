
insts = OrderedDict()

insts["test"] = [
    ((1, 2, 15, 50.0, true),),
    ((1, 2, 15, 50.0, false),),
    ((3, 1, 40, 0.0, true),),
    ((3, 1, 40, 0.0, false),),
    ((2, 2, 40, 100.0, true),),
]

return (PolyRegression, insts)
