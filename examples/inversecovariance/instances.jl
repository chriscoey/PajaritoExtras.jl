
insts = OrderedDict()

insts["test"] = [
    ((4, MatNegLog()),),
    ((4, MatNegLogDirect()),),
    ((6, MatNegRtdet()),),
    ((6, MatNegRtdetEFExp()),),
]

inversecovariance_insts(f::MatSpecExt, max_d::Int) = [((d, f),) for d in vcat(5, 5:max_d)]

insts["nat"] = inversecovariance_insts(MatNegLog(), 16)
insts["ext"] = inversecovariance_insts(MatNegLogDirect(), 16)

return (InverseCovariance, insts)
