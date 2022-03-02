
nosubp_options = (; solve_subproblems = false) # needed for SOS2 in MOIPajarito

insts = OrderedDict()

insts["test"] = [
    ((1, 2, 3),),
    ((2, 2, 2),),
    ((1, 3, 3, true, CCBounded(), 10),),
    ((1, 3, 3, true, LogIBBounded(), 10),),
    ((1, 3, 3, true, SOS2(), 10), nosubp_options),
    ((1, 3, 3, true, CCBounded(), 10), nosubp_options),
    ((1, 3, 3, true, LogIBBounded(), 10), nosubp_options),
    # ((1, 3, 3, true, MC(), 10),),
    # ((1, 3, 3, true, Log(), 10),),
    # ((1, 3, 3, true, ZigZag(), 10),),
]

return (NonconvexRelax, insts)
