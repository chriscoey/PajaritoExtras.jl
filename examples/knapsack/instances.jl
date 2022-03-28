
noext_options = (; use_extended_form = false)
noconic_options = (; solve_subproblems = false, solve_relaxation = false)
noext_noconic_options = merge(noext_options, noconic_options)

insts = OrderedDict()

insts["test"] = [
    ((3, VecNegRtdet(), false),),
    ((3, VecNegRtdet(), true),),
    ((3, VecNegRtdet(), true), noconic_options),
    ((3, VecNegRtdet(), false), noext_options),
    ((3, VecNegRtdet(), false), noext_noconic_options),
    ((3, VecNegRtdetEFExp(), false),),
    ((5, VecNegSqrtConj(), false),),
    ((5, VecNegSqrtConj(), true),),
    ((5, VecNegSqrtConjEF(), false),),
]

function knapsack_insts(f::VecSpecExt, relax_int::Bool, max_n::Int, options::NamedTuple)
    return [((n, f, relax_int), options) for n in vcat(2, 2:2:max_n)]
end

# continuous, separation only:
insts["cont_geo"] = knapsack_insts(VecNegRtdet(), true, 50, noconic_options)
insts["cont_noext_geo"] = knapsack_insts(VecNegRtdet(), true, 50, noext_noconic_options)

# integer, use Hypatia:
insts["nat_geo"] = knapsack_insts(VecNegRtdet(), false, 30, (;))
insts["noext_geo"] = knapsack_insts(VecNegRtdet(), false, 30, noext_options)
insts["ext_geo"] = knapsack_insts(VecNegRtdetEFExp(), false, 30, (;))
# insts["log_nat"] = knapsack_insts(VecNegLog(), true, 10)
# insts["log_noext"] = knapsack_insts(VecNegLog(), false, 10)
# insts["log_ext"] = knapsack_insts(VecNegLogEF(), false, 10)

return (Knapsack, insts)
