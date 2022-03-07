
insts = OrderedDict()

insts["test"] = [
    ((3, nothing, MatNegLog()),),
    ((3, nothing, MatNegLogDirect()),),
    ((5, false, MatNegSqrt()),),
    ((5, false, MatNegSqrtEigOrd()),),
    ((10, true, MatNegRtdet()),),
    ((10, true, MatNegRtdetEFExp()),),
]

inversecovariance_insts(f::MatSpecExt) = [((d, nothing, f),) for d in vcat(3, 4:4:16)]

insts["nat"] = inversecovariance_insts(MatNegLog())
insts["ext"] = inversecovariance_insts(MatNegLogDirect())

return (InverseCovariance, insts)
