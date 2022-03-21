
insts = OrderedDict()

insts["test"] = [
    ((true, false, 3, 4, true),),
    ((true, false, 3, 4, false),),
    ((false, false, 3, 4, true),),
    ((false, false, 3, 4, false),),
    ((true, true, 6, 6, true),),
    ((true, true, 6, 6, false),),
    ((false, true, 6, 6, true),),
    ((false, true, 6, 6, false),),
]

function matrixcompletion_insts(use_nat::Bool, nrow_max::Int)
    return [((true, true, nrow, nrow, use_nat),) for nrow in vcat(20, 20:20:nrow_max)]
end

insts["nat"] = matrixcompletion_insts(true, 220)
insts["ext"] = matrixcompletion_insts(false, 120)

return (MatrixCompletion, insts)
