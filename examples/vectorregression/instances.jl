
insts = OrderedDict()
insts["test"] = [
    ((2, 3, false, VecSpecExt[]),),
    ((3, 3, true, VecSpecExt[]),),
    ((3, 3, false, [VecNegRtdet()]),),
    ((3, 2, true, [VecNegEntropy(), VecNegSqrtConj()]),),
]
return (VectorRegression, insts)
