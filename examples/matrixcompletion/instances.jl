
insts = OrderedDict()
insts["test"] = [
    ((true, false, 3, 4, true),),
    ((true, false, 3, 4, false),),
    ((false, false, 3, 4, true),),
    ((false, false, 3, 4, false),),
    ((true, true, 6, 6, true),),
    ((true, true, 6, 6, false),),
    ((false, true, 6, 6, true),),
    ((false, true, 6, 6, false),),
]
return (MatrixCompletion, insts)
