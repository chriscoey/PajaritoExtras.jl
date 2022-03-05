
nosubp_options = (; solve_subproblems = false) # needed for SOS2 in MOIPajarito

insts = OrderedDict()

insts["test"] = [
    ((1, 1, 3),),
    ((2, 3, 4),),
    ((3, 2, 5),),
    ((1, 6, 10), (; solve_relaxation = false)),
    ((1, 3, 6), (; solve_relaxation = false)),
    ((1, 2, 3, true, false, SOS2(), 10), nosubp_options),
    ((1, 2, 3, true, true, SOS2(), 10), nosubp_options),
    ((1, 3, 6, true, false, CCBounded(), 10),),
    ((1, 3, 6, true, true, CCBounded(), 40),),
    ((2, 2, 4, true, false, SOS2(), 20), nosubp_options),
    ((2, 2, 4, true, false, CCBounded(), 20),),
    ((2, 2, 4, true, false, CCBounded(), 20), nosubp_options),
    ((2, 2, 6, true, true, LogIBBounded(), 40),),
    ((2, 2, 6, true, true, LogIBBounded(), 40), nosubp_options),
    ((1, 6, 5, true, false, SOS2(), 20), nosubp_options),
]

return (ModularDesign, insts)
