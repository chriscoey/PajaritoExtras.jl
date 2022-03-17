
insts = OrderedDict()

insts["test"] = [
    ((4, MatNegLog()),),
    ((4, MatNegLogDirect()),),
    ((6, MatNegRtdet()),),
    ((6, MatNegRtdetEFExp()),),
]

inversecovariance_insts(f::MatSpecExt) = [((d, f),) for d in vcat(5, 5:2:25)]

insts["nat"] = inversecovariance_insts(MatNegLog())
insts["ext"] = inversecovariance_insts(MatNegLogDirect())

return (InverseCovariance, insts)
