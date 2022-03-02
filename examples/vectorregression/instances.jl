
insts = OrderedDict()

insts["test"] = [
    ((2, 3, false, VecSpecExt[]),),
    ((3, 3, true, VecSpecExt[]),),
    ((5, 5, false, [VecNegRtdet()]),),
    ((6, 4, true, [VecNegEntropy(), VecNegSqrtConj()]),),
]

return (VectorRegression, insts)
