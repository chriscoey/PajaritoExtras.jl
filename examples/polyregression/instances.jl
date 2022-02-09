
insts = OrderedDict()

insts["test"] = [
    ((1, 2, 15, true),),
    ((1, 2, 15, false),),
    ((3, 1, 40, true),),
    ((3, 1, 40, false),),
    ((2, 2, 40, true),),
]

return (PolyRegression, insts)
