
nosubp_options = (; solve_subproblems = false) # needed for SOS2 in Pajarito

insts = OrderedDict()

insts["test"] = [
    ((2, 2, SOS2(), 21), nosubp_options),
    ((2, 2, CCBounded(), 21), nosubp_options),
    ((2, 2, LogIBBounded(), 21), nosubp_options),
    ((3, 3, SOS2(), 11), nosubp_options),
    ((3, 3, CCBounded(), 11), nosubp_options),
    ((3, 3, LogIBBounded(), 11), nosubp_options),
    # ((1, 2, CCNonneg(), 11), nosubp_options),
    # ((1, 2, LogIBNonneg(), 11), nosubp_options),
]

function ballpacking_insts(pwl::PWLSOS2, max_p::Int)
    ts = [((4, 3, pwl, 7), nosubp_options)] # compile
    for p in 3:max_p
        push!(ts, ((4, 3, pwl, 2^p - 1), nosubp_options))
    end
    return ts
end

insts["sos2"] = ballpacking_insts(SOS2(), 10)
insts["logib"] = ballpacking_insts(LogIBBounded(), 13)
insts["cc"] = ballpacking_insts(CCBounded(), 7)

return (BallPacking, insts)
