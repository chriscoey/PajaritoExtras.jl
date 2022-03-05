
nosubp_options = (; solve_subproblems = false) # needed for SOS2 in MOIPajarito

insts = OrderedDict()

insts["test"] = [
    ((2, 2, SOS2(), 21), nosubp_options),
    ((2, 2, CCBounded(), 21), nosubp_options),
    ((2, 2, LogIBBounded(), 21), nosubp_options),
    ((3, 3, SOS2(), 11), nosubp_options),
    ((3, 3, CCBounded(), 11), nosubp_options),
    ((3, 3, LogIBBounded(), 11), nosubp_options),
    # ((2, 2, CCUnbounded(), 11), nosubp_options),
    # ((2, 2, LogIBUnbounded(), 11), nosubp_options),
]

return (BallPacking, insts)
