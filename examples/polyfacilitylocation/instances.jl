
sparse_options = (; conic_solver = sparse_hypatia)

insts = OrderedDict()

insts["test"] = [
    ((1, 2, 2, true),),
    ((1, 2, 2, false),),
    ((2, 3, 3, true), sparse_options),
    ((2, 3, 3, false), sparse_options),
    ((6, 6, 3, true), sparse_options),
    ((6, 6, 3, false), sparse_options),
    ((20, 20, 3, true), sparse_options),
]

function polyfacilitylocation_insts(use_nat::Bool, n_max::Int)
    return [((n, 2n, 3, use_nat), sparse_options) for n in vcat(5, 5:5:n_max)]
end

insts["nat"] = polyfacilitylocation_insts(true, 50)
insts["ext"] = polyfacilitylocation_insts(false, 25)

return (PolyFacilityLocation, insts)
