
insts = OrderedDict()

insts["test"] = [
    ((3, nothing, MatNegLog()),),
    ((3, nothing, MatNegLogDirect()),),
    ((5, false, MatNegSqrt()),),
    ((5, false, MatNegSqrtEigOrd()),),
    ((10, true, MatNegRtdet()),),
    ((10, true, MatNegRtdetEFExp()),),
]

return (InverseCovariance, insts)
