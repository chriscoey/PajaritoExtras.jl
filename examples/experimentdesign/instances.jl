
insts = OrderedDict()

insts["test"] = [
    # rootdet
    ((4, MatNegRtdet()),),
    ((4, MatNegRtdetEFExp()),),
    ((4, MatNegRtdetEFPow()),),
    # tr neglog
    ((4, MatNegLog()),),
    ((4, MatNegLogDirect()),),
    # tr negentropy
    ((3, MatNegEntropy()),),
    ((3, MatNegEntropyEigOrd()),),
    # tr negsqrt
    ((2, MatNegSqrt()),),
    ((2, MatNegSqrtEigOrd()),),
    # tr negpower01
    ((3, MatNegPower01(0.7)),),
    ((3, MatNegPower01EigOrd(0.7)),),
    # tr power12
    ((2, MatPower12(1.3)),),
    ((2, MatPower12EigOrd(1.3)),),
]

experimentdesign_insts(f::MatSpecExt) = [((d, f),) for d in vcat(3, 4:4:16)]

insts["nat_rtdet"] = experimentdesign_insts(MatNegRtdet())
insts["ext_rtdet"] = experimentdesign_insts(MatNegRtdetEFExp())
insts["nat_entr"] = experimentdesign_insts(MatNegEntropy())
insts["ext_entr"] = experimentdesign_insts(MatNegEntropyEigOrd())

return (ExperimentDesign, insts)
