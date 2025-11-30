"""
    ERA5 LES - Primitive Variables

    Convert between conservative and primitive variables.

    Author: Jexpresso Development Team
    Date: 2025-11-30
"""

function user_primitives!(qp, q, qe, ip, ::TOTAL)
    # Total formulation
    ρ = q[ip,1]

    qp[ip,1] = q[ip,1]           # ρ
    qp[ip,2] = q[ip,2]/ρ         # u
    qp[ip,3] = q[ip,3]/ρ         # v
    qp[ip,4] = q[ip,4]/ρ         # w
    qp[ip,5] = q[ip,5]/ρ         # θ
    qp[ip,6] = q[ip,6]/ρ         # qt
    qp[ip,7] = q[ip,7]/ρ         # ql
    qp[ip,8] = q[ip,end]         # P
end

function user_primitives!(qp, q, qe, ip, ::PERT)
    # Perturbation formulation
    ρ = q[ip,1] + qe[ip,1]

    qp[ip,1] = ρ                 # ρ (total)
    qp[ip,2] = (q[ip,2] + qe[ip,2])/ρ  # u
    qp[ip,3] = (q[ip,3] + qe[ip,3])/ρ  # v
    qp[ip,4] = (q[ip,4] + qe[ip,4])/ρ  # w
    qp[ip,5] = (q[ip,5] + qe[ip,5])/ρ  # θ
    qp[ip,6] = (q[ip,6] + qe[ip,6])/ρ  # qt
    qp[ip,7] = (q[ip,7] + qe[ip,7])/ρ  # ql
    qp[ip,8] = q[ip,end]               # P
end
