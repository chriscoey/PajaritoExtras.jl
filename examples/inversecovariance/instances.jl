
insts = OrderedDict()
insts["test"] = [
    ((3, MatNegLog(), nothing),),
    ((3, MatNegLogDirect(), nothing),),
    ((5, MatNegRtdet(), true),),
    ((10, MatNegSqrt(), false),),
]
return (InverseCovariance, insts)
