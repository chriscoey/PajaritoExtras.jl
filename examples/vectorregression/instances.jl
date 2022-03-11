
insts = OrderedDict()

insts["test"] = [
    ((2, 3, false, VecSpecExt[]),),
    ((3, 3, true, VecSpecExt[]),),
    ((5, 5, false, [VecNegRtdet()]),),
    ((6, 4, true, [VecNegEntropy(), VecNegSqrtConj()]),),
]

function vectorregression_insts(use_extended_form::Bool)
    options = (; use_extended_form = use_extended_form)
    hs = [VecNegEntropy(), VecNegSqrtConj()]
    return [((n, n, false, hs), options) for n in vcat(10, 15:15:150)]
end

insts["nat"] = vectorregression_insts(true)
insts["nat_noext"] = vectorregression_insts(false)

return (VectorRegression, insts)
