
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

function experimentdesign_insts(specs::Vector)
    ts = Tuple[]
    for (max_d, f) in specs
        t = [((d, f),) for d in vcat(2, 3:3:max_d)]
        append!(ts, t)
    end
    return ts
end

insts["nat"] = experimentdesign_insts([
    (15, MatNegRtdet()),
    (15, MatNegLog()),
    (15, MatNegSqrtConj()),
    (15, MatNegPower01(1 / 3)),
])

insts["ext"] = experimentdesign_insts([
    (15, MatNegRtdetEFExp()),
    (15, MatNegLogDirect()),
    (15, MatNegSqrtConjDirect()),
    (15, MatNegPower01EigOrd(1 / 3)),
])

return (ExperimentDesign, insts)
