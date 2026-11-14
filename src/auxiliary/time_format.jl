"""
    time_unit(t) -> (value, unit)

Model time `t` (seconds) in the largest unit that keeps it ≥ 1, for the run
log: `(90.0, "d")` for `t = 7.776e6`, `(1.5, "h")` for `5400`, `(2.5, "min")`
for `150`, `(12.5, "s")` for `12.5`. Days are the largest unit on purpose: the
length of a model year differs between cases (e.g. 360 d in the ocean decks).

Allocation-free — the unit is a string literal and nothing is formatted here —
so a log line can carry it as one more `@printf` field at no extra cost:

    @printf(" # ... at t=%.6g s (%.2f %s)\\n", t, time_unit(t)...)

(checked in test/test_time_format.jl).
"""
@inline function time_unit(t::Real)
    x = Float64(t)
    a = abs(x)
    a >= 86400.0 && return (x / 86400.0, "d")
    a >= 3600.0  && return (x / 3600.0,  "h")
    a >= 60.0    && return (x / 60.0,    "min")
    return (x, "s")
end
