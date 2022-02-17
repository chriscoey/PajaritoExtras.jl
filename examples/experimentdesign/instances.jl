
insts = OrderedDict()

insts["test"] = [
    # rootdet
    ((4, MatNegRtdet()),),
    ((4, MatNegRtdetEFExp()),),
    ((4, MatNegRtdetEFPow()),),
    # tr neglog
    ((5, MatNegLog()),),
    ((5, MatNegLogDirect()),),
    # tr negentropy
    ((3, MatNegEntropy()),),
    ((3, MatNegEntropyEigOrd()),),
    # tr negsqrt
    ((3, MatNegSqrt()),),
    ((3, MatNegSqrtEigOrd()),),
    # tr negpower01
    ((3, MatNegPower01(0.7)),),
    ((3, MatNegPower01EigOrd(0.7)),),
    # tr power12
    ((3, MatPower12(1.3)),),
    ((3, MatPower12EigOrd(1.3)),),
]

experimentdesign_insts(ext::MatSpecExt) = [((d, ext),) for d in vcat(3, 10:10:10)]

insts["nat"] = vcat(experimentdesign_insts.([
    # MatNegRtdet(),
    MatNegLog(), # TODO or MatLogdetCone if supported
    # MatNegSqrtConj(),
    # MatNegPower01(1/3),
])...)

insts["ext"] = vcat(experimentdesign_insts.([
    # MatNegRtdetEFExp(),
    MatNegLogDirect(),
    # MatNegSqrtConjDirect(),
    # MatNegPower01EigOrd(1/3),
])...)

return (ExperimentDesign, insts)
