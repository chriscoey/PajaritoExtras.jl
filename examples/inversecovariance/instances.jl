
insts = OrderedDict()

insts["test"] = [
    ((3, nothing, MatNegLog()),),
    ((3, nothing, MatNegLogDirect()),),
    ((5, false, MatNegSqrt()),),
    ((5, false, MatNegSqrtEigOrd()),),
    ((10, true, MatNegRtdet()),),
    ((10, true, MatNegRtdetEFExp()),),
]

function inversecovariance_insts(specs::Vector)
    ts = Tuple[]
    for (max_d, f) in specs
        t = [((d, nothing, f),) for d in vcat(3, 3:3:max_d)]
        append!(ts, t)
    end
    return ts
end

insts["nat"] = inversecovariance_insts([(12, MatNegLog())])
insts["ext"] = inversecovariance_insts([(12, MatNegLogDirect())])

return (InverseCovariance, insts)
