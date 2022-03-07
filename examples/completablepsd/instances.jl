
insts = OrderedDict()

insts["test"] = [
    ((3, 0.7, 0.5, true),),
    ((3, 0.7, 0.5, false),),
    ((6, 0.3, 0.7, true),),
    ((6, 0.3, 0.7, false),),
    ((10, 0.7, 0.2, true),),
    ((10, 0.7, 0.2, false),),
    ((15, 0.2, 0.3, true),),
    ((15, 0.2, 0.3, false),),
]

function completablepsd_insts(use_nat::Bool)
    return [((d, inv(sqrt(d)), 0.7, use_nat),) for d in vcat(10, 10:10:100)]
end

insts["nat"] = completablepsd_insts(true)
insts["ext"] = completablepsd_insts(false)

return (CompletablePSD, insts)
