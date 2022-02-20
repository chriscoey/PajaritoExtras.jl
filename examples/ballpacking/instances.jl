
insts = OrderedDict()

insts["test"] = [
    # ((3, 3, SOS2(), 21),), # TODO MOIPajarito can't handle actual SOS2 constraints
    ((2, 2, CCBounded(), 11),),
    ((2, 2, CCUnbounded(), 11),),
    ((2, 2, LogIBBounded(), 11),),
    ((2, 2, LogIBUnbounded(), 11),),
    # ((3, 3, MC(), 21),),
    # ((3, 3, Log(), 21),),
    # ((3, 3, ZigZag(), 21),),
]

return (BallPacking, insts)
