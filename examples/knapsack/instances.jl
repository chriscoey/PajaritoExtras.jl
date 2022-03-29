
noext_options = (; use_extended_form = false)
noconic_options =
    (; solve_subproblems = false, solve_relaxation = false, use_init_fixed_oa = true)
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

# continuous, separation only:
function cont_knapsack_insts(f::VecSpecExt, max_n::Int, options::NamedTuple)
    return [((n, f, true), options) for n in vcat(5, 5:5:max_n)]
end

insts["cont_geo"] = cont_knapsack_insts(VecNegRtdet(), 200, noconic_options)
insts["cont_noext_geo"] = cont_knapsack_insts(VecNegRtdet(), 70, noext_noconic_options)
insts["cont_log"] = cont_knapsack_insts(VecNegLog(), 200, noconic_options)
insts["cont_noext_log"] = cont_knapsack_insts(VecNegLog(), 70, noext_noconic_options)
insts["cont_inv"] = cont_knapsack_insts(VecNegSqrtConj(), 200, noconic_options)
insts["cont_noext_inv"] = cont_knapsack_insts(VecNegSqrtConj(), 70, noext_noconic_options)

# integer:
function int_knapsack_insts(f::VecSpecExt, max_n::Int, options::NamedTuple)
    return [((n, f, false), options) for n in vcat(3, 3:3:15, 20:5:35, 40:10:max_n)]
end

insts["nat_geo"] = int_knapsack_insts(VecNegRtdet(), 140, (;))
insts["noext_geo"] = int_knapsack_insts(VecNegRtdet(), 25, noext_options)
insts["ext_geo"] = int_knapsack_insts(VecNegRtdetEFExp(), 140, (;))
# insts["nat_log"] = int_knapsack_insts(VecNegLog(), 75, (;))
# insts["noext_log"] = int_knapsack_insts(VecNegLog(), 25, noext_options)
# insts["ext_log"] = int_knapsack_insts(VecNegLogEF(), 75, (;))

return (Knapsack, insts)
