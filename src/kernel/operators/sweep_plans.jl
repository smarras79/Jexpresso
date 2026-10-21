const USE_FACE_PLANS = Ref(true)

struct FinePlan
    remote_gel::TInt
    gids::Matrix{Int}                      # (ngl,ngl) parent-grid gids, 0 = no entry
    L1::Matrix{Float64}                    # (ngl^2, ngl)
    L2::Matrix{Float64}
end

struct CoarseContrib
    rg::TInt                       # remote child gel
    exact::Bool                    # true: node coincides with a child node, copy its value
    gid::Int                       # gid of that child node (exact case)
    gids::Vector{Int}              # child grid gids, index q1 + (q2-1)*ngl (interpolated case)
    L1::Vector{Float64}            # Lagrange row along axis 1 (interpolated case)
    L2::Vector{Float64}            # Lagrange row along axis 2 (interpolated case)
end

struct CoarsePlan
    nodes::Vector{Vector{CoarseContrib}}    # per coarse face node, contributions in global-gel order
end

const FINE_PLANS   = IdDict{Any,Dict{Tuple{Int,Int},FinePlan}}()      # keyed by the remote_gid_coords object's identity
const COARSE_PLANS = IdDict{Any,Dict{Tuple{Int,Int},CoarsePlan}}()

const USE_LOCAL_NCF_PLANS = Ref(true)
const USE_FLAT_FACES = Ref(true)

# ---------- per-branch closure-evaluation timers ----------
const BRANCH_EVAL_LABELS = vcat(BRANCH_LABELS, [:none])      # :none = faces without an upwind closure
const BRANCH_EVAL = Dict(b => Symbol("fe_$(b)") for b in BRANCH_EVAL_LABELS)
for b in BRANCH_EVAL_LABELS; push!(PROF_SECTIONS, BRANCH_EVAL[b]); end

# ---------- local (same-rank) NCF plans ----------
struct LocalFinePlan
    plidx::Vector{Int32}          # parent_grid entry idx -> local index in the parent vector (0 = unset)
    skip::BitVector               # child node skipped (value 0.0), as in the original
    W1::Matrix{Float64}           # (n_child_nodes, ngl): L1[cj, :]
    W2::Matrix{Float64}           # (n_child_nodes, ngl): L2[ci, :]
end
struct LocalChildSlot
    child_iel::Int
    clidx::Vector{Int32}          # child_grid entry idx -> local index in the child vector (0 = unset)
    covered::BitVector            # parent face node covered by this child
    R1::Matrix{Float64}           # (ngl^2, ngl) Lagrange rows of the covered nodes
    R2::Matrix{Float64}
end
struct LocalCoarsePlan
    children::Vector{Int}         # child elements in global-element-id order
    slots::Vector{LocalChildSlot}
end

const LOCAL_FINE_PLANS   = IdDict{Any,Dict{Int,LocalFinePlan}}()                  # keyed by the identity of ncf_L1L2
const LOCAL_COARSE_PLANS = IdDict{Any,Dict{Tuple{Int,Int},LocalCoarsePlan}}()
const FLAT_CORR = IdDict{Any,Dict{Tuple{Int,Int},Array{Int32,3}}}()               # keyed by the identity of face_correspondence

const USE_ANGULAR_PLANS = Ref(true)
const ANG_STATS = zeros(Int, 5)   # same-rank built, same-rank evaluated, cross-rank built, cross-rank evaluated, slot tables built

struct AngularSamePlan
    jel::Int
    lids::Vector{Int}          # local index of each face node (face_ip_range loop order)
    nidx::Vector{Int}          # matching local index in the neighbour element's vectors
    src_ips::Vector{Int}       # nsrc^2 source directions, a outer / b inner
    w::Vector{Float64}         # Lθ[a]*Lϕ[b], same order
end
mutable struct AngularRemotePlan
    remote_gel::TInt
    lids::Vector{Int}
    gids::Vector{Int}          # global node id of each face node
    src_ips::Vector{Int}
    w::Vector{Float64}
    slots::Matrix{Int32}       # (nsrc^2, nface) receive-buffer slots; filled on first use with a GidTable
    slot_owner::Any
end

const ANGULAR_SAME_PLANS   = IdDict{Any,Dict{Tuple{Int,Int,Int},AngularSamePlan}}()                      # keyed by identity of angular_deps
const ANGULAR_REMOTE_PLANS = IdDict{Any,Dict{Tuple{Int,Int,Int},Union{Nothing,AngularRemotePlan}}}()     # keyed by identity of angular_deps_cr

function clear_face_plans!()      # REPLACES the existing definition
    empty!(FINE_PLANS); empty!(COARSE_PLANS)
    empty!(LOCAL_FINE_PLANS); empty!(LOCAL_COARSE_PLANS); empty!(FLAT_CORR)
    empty!(ANGULAR_SAME_PLANS); empty!(ANGULAR_REMOTE_PLANS)
    fill!(ANG_STATS, 0)
    return nothing
end

function make_local_fine_plan(mesh, parent_iel::Int, ncf_idx::Int, L1, L2, ngl::Int)
    IPc = mesh.IPc_list[:,ncf_idx]; IPp = mesh.IPp_list[:,ncf_idx]
    pidx = Dict{Int,Int}()
    for k=1:ngl, j=1:ngl, i=1:ngl; pidx[mesh.connijk[parent_iel,i,j,k]] = lidx(i,j,k,ngl); end
    plidx = zeros(Int32, length(IPp))
    for idx in 1:length(IPp)
        node = Int(IPp[idx]); node <= 0 && continue
        plidx[idx] = pidx[node]
    end
    nc = length(IPc); skip = trues(nc); W1 = zeros(nc, ngl); W2 = zeros(nc, ngl)
    for n in 1:nc
        Int(IPc[n]) <= 0 && continue
        ci = div(n-1, ngl) + 1; cj = mod(n-1, ngl) + 1
        (ci < 1 || ci > size(L1,1) || cj < 1 || cj > size(L2,1)) && continue
        skip[n] = false; W1[n,:] .= L1[cj,:]; W2[n,:] .= L2[ci,:]
    end
    return LocalFinePlan(plidx, skip, W1, W2)
end

function eval_local_fine_plan!(out, pl::LocalFinePlan, parent_vals, ngl::Int)
    pg = zeros(ngl, ngl)
    for idx in eachindex(pl.plidx)
        l = pl.plidx[idx]; l == 0 && continue
        pg[((idx-1)%ngl)+1, ((idx-1)÷ngl)+1] = parent_vals[l]
    end
    for n in eachindex(pl.skip)
        if pl.skip[n]; out[n] = 0.0; continue; end
        val = 0.0
        for p1 = 1:ngl, p2 = 1:ngl
            val += pl.W1[n,p1]*pl.W2[n,p2]*pg[p1,p2]
        end
        out[n] = val
    end
    return out
end

function make_local_coarse_plan(mesh, ncf_group, ngl::Int)
    order = sort(collect(ncf_group); by = idx -> mesh.el2gel[Int(mesh.cip[idx])])     # same order as the same-rank builder
    slots = LocalChildSlot[]
    for ncf_idx in order
        child_iel = Int(mesh.cip[ncf_idx]); IPc = mesh.IPc_list[:,ncf_idx]; IPp = mesh.IPp_list[:,ncf_idx]
        child_ijk = Dict{Int,NTuple{3,Int}}()
        for k=1:ngl, j=1:ngl, i=1:ngl; child_ijk[mesh.connijk[child_iel,i,j,k]] = (i,j,k); end
        line1 = Int[]; line2 = Int[]
        for i = 1:ngl; idx = (i-1)*ngl + 1; idx <= length(IPc) || continue; push!(line1, Int(IPc[idx])); end
        for j = 1:ngl; idx = j;             idx <= length(IPc) || continue; push!(line2, Int(IPc[idx])); end
        axis_c1 = axis_of_max_spread(([mesh.x[n] for n in line2], [mesh.y[n] for n in line2], [mesh.z[n] for n in line2]))
        axis_c2 = axis_of_max_spread(([mesh.x[n] for n in line1], [mesh.y[n] for n in line1], [mesh.z[n] for n in line1]))
        coord(node, axis) = axis==1 ? mesh.x[node] : (axis==2 ? mesh.y[node] : mesh.z[node])
        c1_coords = [coord(n,axis_c1) for n in line2]; c2_coords = [coord(n,axis_c2) for n in line1]
        c1_min,c1_max = extrema(c1_coords); c2_min,c2_max = extrema(c2_coords)
        tol1 = 1e-6*max(c1_max-c1_min,1.0); tol2 = 1e-6*max(c2_max-c2_min,1.0)
        clidx = zeros(Int32, ngl*ngl)
        for idx in 1:length(IPc)
            node = Int(IPc[idx]); node <= 0 && continue
            clidx[idx] = lidx(child_ijk[node]..., ngl)
        end
        covered = falses(ngl*ngl); R1 = zeros(ngl*ngl, ngl); R2 = zeros(ngl*ngl, ngl)
        for idx in 1:length(IPp)
            node = Int(IPp[idx]); node <= 0 && continue
            pc1 = coord(node,axis_c1); pc2 = coord(node,axis_c2)
            (pc1 < c1_min-tol1 || pc1 > c1_max+tol1 || pc2 < c2_min-tol2 || pc2 > c2_max+tol2) && continue
            covered[idx] = true
            R1[idx,:] .= vec(build_1d_lagrange_matrix(c1_coords, [pc1])[1,:])
            R2[idx,:] .= vec(build_1d_lagrange_matrix(c2_coords, [pc2])[1,:])
        end
        push!(slots, LocalChildSlot(child_iel, clidx, covered, R1, R2))
    end
    return LocalCoarsePlan([s.child_iel for s in slots], slots)
end

function eval_local_coarse_plan(pl::LocalCoarsePlan, x_dg, d::Int, ngl::Int)
    psum = zeros(ngl, ngl); pcnt = zeros(Int, ngl, ngl); cg = zeros(ngl, ngl)
    for s in pl.slots
        cv = x_dg[(s.child_iel, d)]
        fill!(cg, 0.0)
        for idx in 1:ngl*ngl
            l = s.clidx[idx]; l == 0 && continue
            cg[((idx-1)%ngl)+1, ((idx-1)÷ngl)+1] = cv[l]
        end
        for idx in 1:ngl*ngl
            s.covered[idx] || continue
            val = 0.0
            for a = 1:ngl, b = 1:ngl
                val += s.R1[idx,a]*s.R2[idx,b]*cg[a,b]
            end
            p1 = ((idx-1)%ngl)+1; p2 = ((idx-1)÷ngl)+1
            psum[p1,p2] += val; pcnt[p1,p2] += 1
        end
    end
    return [pcnt[a,b] > 0 ? psum[a,b]/pcnt[a,b] : NaN for a in 1:ngl, b in 1:ngl]
end

function local_fine_plan!(mesh, ncf_L1L2, ncf_idx::Int, parent_iel::Int, ngl::Int)
    dct = get!(() -> Dict{Int,LocalFinePlan}(), LOCAL_FINE_PLANS, ncf_L1L2)
    return get!(dct, ncf_idx) do
        L1, L2 = ncf_L1L2[ncf_idx]
        make_local_fine_plan(mesh, parent_iel, ncf_idx, L1, L2, ngl)
    end
end
function local_coarse_plan!(mesh, ncf_L1L2, iel::Int, fid::Int, ncf_group, ngl::Int)
    dct = get!(() -> Dict{Tuple{Int,Int},LocalCoarsePlan}(), LOCAL_COARSE_PLANS, ncf_L1L2)
    return get!(() -> make_local_coarse_plan(mesh, ncf_group, ngl), dct, (iel,fid))
end

# ---------- flat index arrays for same-rank conforming faces ----------
function flat_corr_index!(face_correspondence, iel::Int, fid::Int, ngl::Int)
    dct = get!(() -> Dict{Tuple{Int,Int},Array{Int32,3}}(), FLAT_CORR, face_correspondence)
    return get!(dct, (iel, fid)) do
        inner = face_correspondence[(iel, fid)]
        arr = zeros(Int32, ngl, ngl, ngl)
        irange, jrange, krange = face_ip_range(fid, ngl)
        for k in krange, j in jrange, i in irange
            c = inner[(i, j, k)]
            arr[i, j, k] = lidx(c[1], c[2], c[3], ngl)
        end
        arr
    end
end

# ---------- verification against the original builders (synthetic field; any rank count) ----------
function verify_local_ncf_plans(setup, comm; n_dirs_sample::Int=3)
    mesh = setup.mesh; ngl = mesh.ngl
    ds = unique(round.(Int, range(1, minimum(setup.n_dirs); length=n_dirs_sample)))
    need = Set{Int}()
    for (_, idx) in setup.ncf_as_child; push!(need, Int(mesh.pip[idx])); end
    for (_, grp) in setup.ncf_as_parent, idx in grp; push!(need, Int(mesh.cip[idx])); end
    x_dg = Dict{Tuple{Int,Int},Vector{Float64}}()
    for iel in need, d in ds
        gel = Float64(mesh.el2gel[iel])
        x_dg[(iel,d)] = [sin(0.37*gel + 1.1*d + 0.013*p) for p in 1:ngl^3]
    end
    nf = 0; nc = 0; bad = 0; maxd = 0.0
    clear_face_plans!()
    for ((c_iel, c_fid), idx) in setup.ncf_as_child, d in ds
        p_iel = Int(mesh.pip[idx]); L1, L2 = setup.ncf_L1L2[idx]
        ref = build_dg_ncf_parent_to_child_local(mesh, x_dg, c_iel, p_iel, c_fid, d, mesh.IPc_list[:,idx], mesh.IPp_list[:,idx], L1, L2, ngl)
        pl = local_fine_plan!(mesh, setup.ncf_L1L2, idx, p_iel, ngl)
        out = Vector{Float64}(undef, length(pl.skip)); eval_local_fine_plan!(out, pl, x_dg[(p_iel,d)], ngl)
        nf += 1
        reinterpret(UInt64, ref) == reinterpret(UInt64, out) || (bad += 1)
        dd = ref .- out; maxd = max(maxd, maximum(abs.(dd[isfinite.(dd)]); init=0.0))
    end
    for ((p_iel, p_fid), grp) in setup.ncf_as_parent, d in ds
        entries = [(idx, Int(mesh.cip[idx]), mesh.IPc_list[:,idx], mesh.IPp_list[:,idx]) for idx in grp]
        ref = build_dg_ncf_child_to_parent_v5_local(mesh, x_dg, p_iel, p_fid, d, entries, ngl)
        pl = local_coarse_plan!(mesh, setup.ncf_L1L2, p_iel, p_fid, grp, ngl)
        out = eval_local_coarse_plan(pl, x_dg, d, ngl)
        nc += 1
        reinterpret(UInt64, vec(ref)) == reinterpret(UInt64, vec(out)) || (bad += 1)
        dd = ref .- out; maxd = max(maxd, maximum(abs.(dd[isfinite.(dd)]); init=0.0))
    end
    cnt = MPI.Allreduce(Float64[nf, nc, bad], +, comm); mx = MPI.Allreduce(maxd, MPI.MAX, comm)
    MPI.Comm_rank(comm) == 0 || return String[]
    mxs = @sprintf("%.2e", mx)
    return String["local NCF plan verification: parent->child $(Int(cnt[1])) evaluations, child->parent $(Int(cnt[2])) evaluations, " *
                  "$(Int(cnt[3])) bitwise mismatches, max |Δ| $mxs (expected 0)"]
end


plan_getc(c, ax) = ax == 1 ? c[1] : (ax == 2 ? c[2] : c[3])

function face_plane_info(mesh, iel::Int, fid::Int, ngl::Int)
    irange, jrange, krange = face_ip_range(fid, ngl)
    coords = NTuple{3,Float64}[]
    for k in krange, j in jrange, i in irange
        ip = mesh.connijk[iel,i,j,k]
        push!(coords, (mesh.x[ip], mesh.y[ip], mesh.z[ip]))
    end
    xs = [c[1] for c in coords]; ys = [c[2] for c in coords]; zs = [c[3] for c in coords]
    spread = (maximum(xs)-minimum(xs), maximum(ys)-minimum(ys), maximum(zs)-minimum(zs))
    excl = argmin(spread)
    normal = excl == 1 ? sum(xs)/length(xs) : (excl == 2 ? sum(ys)/length(ys) : sum(zs)/length(zs))
    axs = setdiff(1:3, [excl])
    return coords, excl, normal, axs[1], axs[2]
end

function make_fine_plan(mesh, child_iel::Int, fid::Int, remote_gel::TInt, remote_gid_coords, ngl::Int;
                        tol::Float64=1e-4, cluster_tol::Float64=1e-6)
    coords, excl, normal, a1, a2 = face_plane_info(mesh, child_iel, fid, ngl)
    cand = [(gid,c) for (gid,c) in remote_gid_coords[remote_gel] if isapprox(plan_getc(c,excl), normal; atol=tol)]
    length(cand) == ngl*ngl || error("fine plan: expected $(ngl*ngl) remote nodes for remote_gel=$remote_gel (iel=$child_iel, fid=$fid), found $(length(cand))")
    p1v, c1 = cluster_values_with_index([plan_getc(c,a1) for (_,c) in cand]; tol=cluster_tol)
    p2v, c2 = cluster_values_with_index([plan_getc(c,a2) for (_,c) in cand]; tol=cluster_tol)
    (length(p1v) == ngl && length(p2v) == ngl) || error("fine plan: degenerate parent grid $(length(p1v))x$(length(p2v)) for remote_gel=$remote_gel")
    gids = zeros(Int, ngl, ngl)
    for (idx,(gid,_)) in enumerate(cand); gids[c1[idx], c2[idx]] = gid; end
    L1 = build_1d_lagrange_matrix(p1v, [plan_getc(c,a1) for c in coords])
    L2 = build_1d_lagrange_matrix(p2v, [plan_getc(c,a2) for c in coords])
    return FinePlan(remote_gel, gids, L1, L2)
end

function eval_fine_plan!(out, pl::FinePlan, d::Int, gr, ngl::Int)          # false = blocked
    vals = Matrix{Float64}(undef, ngl, ngl)
    for i2 in 1:ngl, i1 in 1:ngl
        g = pl.gids[i1,i2]
        if g == 0
            vals[i1,i2] = 0.0
        else
            v = get(gr, (pl.remote_gel, d, g), NaN)
            isnan(v) && return false
            vals[i1,i2] = v
        end
    end
    for node in 1:ngl*ngl
        acc = 0.0
        for p1 in 1:ngl, p2 in 1:ngl
            acc += pl.L1[node,p1] * pl.L2[node,p2] * vals[p1,p2]
        end
        out[node] = acc
    end
    return true
end

function make_coarse_plan(mesh, coarse_iel::Int, fid::Int, remote_gels::Vector{TInt}, remote_gid_coords, ngl::Int;
                          tol::Float64=1e-4, cluster_tol::Float64=1e-6)
    coords, excl, normal, a1, a2 = face_plane_info(mesh, coarse_iel, fid, ngl)

    # ---- the grid of each remote child on this face's plane (unchanged) ----
    child_grids = Dict{TInt,NamedTuple}()
    for rg in remote_gels
        on_plane = [(gid,c) for (gid,c) in remote_gid_coords[rg] if isapprox(plan_getc(c,excl), normal; atol=tol)]
        length(on_plane) == ngl*ngl ||
            error("coarse plan: remote_gel=$rg has $(length(on_plane)) on-plane nodes, expected $(ngl*ngl)")
        v1, cl1 = cluster_values_with_index([plan_getc(c,a1) for (_,c) in on_plane]; tol=cluster_tol)
        v2, cl2 = cluster_values_with_index([plan_getc(c,a2) for (_,c) in on_plane]; tol=cluster_tol)
        (length(v1) == ngl && length(v2) == ngl) || error("coarse plan: degenerate child grid for remote_gel=$rg")
        grid_gid = Dict{Tuple{Int,Int},Int}()
        for (idx,(gid,_)) in enumerate(on_plane); grid_gid[(cl1[idx], cl2[idx])] = gid; end
        child_grids[rg] = (a1_vals=v1, a2_vals=v2, grid_gid=grid_gid,
                           bbox=(minimum(v1)-tol, maximum(v1)+tol, minimum(v2)-tol, maximum(v2)+tol))
    end

    # ---- NEW: for every coarse node, collect ALL children covering it, in canonical (global-gel) order ----
    sorted_rgs = sort(collect(keys(child_grids)))
    nodes = Vector{Vector{CoarseContrib}}(undef, ngl*ngl)
    for idx in 1:ngl*ngl
        cv1 = plan_getc(coords[idx], a1); cv2 = plan_getc(coords[idx], a2)
        list = CoarseContrib[]
        for rg in sorted_rgs
            cg = child_grids[rg]
            (cv1 >= cg.bbox[1] && cv1 <= cg.bbox[2] && cv2 >= cg.bbox[3] && cv2 <= cg.bbox[4]) || continue
            e1 = findfirst(v -> isapprox(v, cv1; atol=tol), cg.a1_vals)
            e2 = findfirst(v -> isapprox(v, cv2; atol=tol), cg.a2_vals)
            if e1 !== nothing && e2 !== nothing
                push!(list, CoarseContrib(rg, true, get(cg.grid_gid, (e1,e2), 0), Int[], Float64[], Float64[]))
            else
                L1 = build_1d_lagrange_matrix(cg.a1_vals, [cv1])
                L2 = build_1d_lagrange_matrix(cg.a2_vals, [cv2])
                gids = zeros(Int, ngl*ngl)
                for q1 in 1:ngl, q2 in 1:ngl
                    gids[q1 + (q2-1)*ngl] = get(cg.grid_gid, (q1,q2), 0)
                end
                push!(list, CoarseContrib(rg, false, 0, gids, vec(L1[1,:]), vec(L2[1,:])))
            end
        end
        isempty(list) && error("coarse plan: node (a1=$cv1, a2=$cv2) not in any child's bounding box")
        nodes[idx] = list
    end
    return CoarsePlan(nodes)
end

@inline function coarse_contrib_value(c::CoarseContrib, d::Int, gr, ngl::Int)
    if c.exact
        c.gid == 0 && return 0.0
        return get(gr, (c.rg, d, c.gid), NaN)
    end
    acc = 0.0
    for q1 in 1:ngl, q2 in 1:ngl
        g = c.gids[q1 + (q2-1)*ngl]; v = 0.0
        if g != 0
            v = get(gr, (c.rg, d, g), NaN)
            isnan(v) && return NaN
        end
        acc += c.L1[q1] * c.L2[q2] * v
    end
    return acc
end

function eval_coarse_plan!(out, pl::CoarsePlan, d::Int, gr, ngl::Int)       # false = blocked (a value has not arrived)
    for idx in 1:ngl*ngl
        lst = pl.nodes[idx]; acc = 0.0
        for c in lst
            v = coarse_contrib_value(c, d, gr, ngl)
            isnan(v) && return false
            acc += v
        end
        out[idx] = acc / length(lst)
    end
    return true
end

function fine_plan!(mesh, rgc, iel::Int, fid::Int, remote_gel::TInt, ngl::Int)
    dct = get!(() -> Dict{Tuple{Int,Int},FinePlan}(), FINE_PLANS, rgc)
    return get!(() -> make_fine_plan(mesh, iel, fid, remote_gel, rgc, ngl), dct, (iel,fid))
end
function coarse_plan!(mesh, rgc, iel::Int, fid::Int, remote_gels::Vector{TInt}, ngl::Int)
    dct = get!(() -> Dict{Tuple{Int,Int},CoarsePlan}(), COARSE_PLANS, rgc)
    return get!(() -> make_coarse_plan(mesh, iel, fid, remote_gels, rgc, ngl), dct, (iel,fid))
end

function verify_face_plans(setup, I_cg_all, comm; max_dirs::Int=6)
    rank = MPI.Comm_rank(comm); mesh = setup.mesh; ngl = mesh.ngl
    gr  = exchange_ncf_gid_values!(setup.gid_exch, I_cg_all, ngl)           # collective; old Dict path supplies the reference data
    rgc = setup.remote_gid_coords

    items = Tuple{Symbol,Int,Int,Vector{TInt},Int}[]
    for ((iel, fid), kind) in setup.neighbor_kind
        (kind == :i_am_fine || kind == :i_am_coarse) || continue
        rgs = setup.neighbor_gel_group[(iel, fid)]
        all(g -> haskey(rgc, g), rgs) || continue
        nd = setup.n_dirs[Int(iel)]
        for d in unique(round.(Int, range(1, nd; length=min(max_dirs, nd))))
            push!(items, (kind, Int(iel), Int(fid), rgs, d))
        end
    end
    plan_of(kind, iel, fid, rgs) = kind == :i_am_fine ? fine_plan!(mesh, rgc, iel, fid, rgs[1], ngl) :
                                                          coarse_plan!(mesh, rgc, iel, fid, rgs, ngl)

    # coverage: how many coarse nodes of the verified faces have several covering children
    n_fine = 0; n_coarse = 0; n_multi = 0; seen = Set{Tuple{Int,Int}}()
    for (kind, iel, fid, rgs, d) in items
        pl = plan_of(kind, iel, fid, rgs)
        if kind == :i_am_fine
            n_fine += 1
        else
            n_coarse += 1
            if !((iel, fid) in seen)
                push!(seen, (iel, fid)); n_multi += count(l -> length(l) > 1, pl.nodes)
            end
        end
    end

    function run_ref()
        out = Vector{Vector{Float64}}(undef, length(items))
        for (n, (kind, iel, fid, rgs, d)) in enumerate(items)
            out[n] = kind == :i_am_fine ?
                build_dg_ncf_parent_to_child_remote_by_coords(mesh, iel, fid, d, rgs[1], rgc, gr, ngl) :
                vec(build_dg_ncf_child_to_parent_remote_by_coords(mesh, iel, fid, d, rgs, rgc, gr, ngl))
        end
        return out
    end
    function run_plan()
        out = Vector{Union{Nothing,Vector{Float64}}}(undef, length(items))
        for (n, (kind, iel, fid, rgs, d)) in enumerate(items)
            buf = Vector{Float64}(undef, ngl*ngl)
            pl = plan_of(kind, iel, fid, rgs)
            ok = kind == :i_am_fine ? eval_fine_plan!(buf, pl, d, gr, ngl) : eval_coarse_plan!(buf, pl, d, gr, ngl)
            out[n] = ok ? buf : nothing
        end
        return out
    end

    run_ref(); run_plan()                                     # warm-up: compile both paths before timing
    t0 = time_ns(); ref = run_ref();  t_ref = (time_ns() - t0) * 1e-9
    t0 = time_ns(); new = run_plan(); t_new = (time_ns() - t0) * 1e-9

    n_bad = 0; n_blk = 0; max_abs = 0.0
    for n in eachindex(items)
        r = ref[n]; p = new[n]
        if p === nothing
            n_blk += 1
            any(isnan, r) || (n_bad += 1)                     # plan blocked although the reference is complete
        elseif any(isnan, r)
            n_bad += 1                                        # plan returned a value where the reference is blocked
        else
            reinterpret(UInt64, r) == reinterpret(UInt64, p) || (n_bad += 1)
            max_abs = max(max_abs, maximum(abs.(r .- p)))
        end
    end

    tot = MPI.Allreduce(Float64[length(items), n_bad, n_blk, t_ref, t_new, n_fine, n_coarse, n_multi], +, comm)
    mx  = MPI.Allreduce(max_abs, MPI.MAX, comm)
    rank == 0 || return String[]

    nev = Int(tot[1]); nfine = Int(tot[6]); ncoarse = Int(tot[7]); nmulti = Int(tot[8])
    nbad = Int(tot[2]); nblk = Int(tot[3])
    us_ref = round(1e6*tot[4]/max(tot[1],1); digits=1)
    us_new = round(1e6*tot[5]/max(tot[1],1); digits=2)
    mxs = @sprintf("%.2e", mx)                       # single literal, so the macro accepts it
    return String["face-plan verification: $nev evaluations (fine $nfine, coarse $ncoarse; " *
                  "coarse faces contain $nmulti multi-child nodes), $nbad bitwise mismatches, " *
                  "$nblk plan-blocked with complete data, max |Δ| $mxs " *
                  "(counts expected 0; max |Δ| expected 0 or about 1e-13); " *
                  "reference $us_ref us/eval, plan $us_new us/eval (second pass, summed over ranks)"]

end


function report_neighbors(setup, comm, ranks_of_interest)
    rank = MPI.Comm_rank(comm)
    nb = Float64.(sort(unique(Int.(collect(values(setup.neighbor_rank))))))
    counts = MPI.Allgather(Int32(length(nb)), comm)
    all = MPI.Allgatherv(nb, Cint.(counts), comm)
    if rank == 0
        offs = cumsum(vcat(0, Int.(counts)))
        for r in ranks_of_interest
            println("rank $r neighbors: ", Int.(all[offs[r+1]+1:offs[r+2]]))
        end
    end
end

function make_angular_same_plan(mesh, iel::Int, fid::Int, jel::Int, d::Int, dep_nodes, face_correspondence, ngl::Int)
    em_jel = mesh.extra_mesh[jel]
    src_e_ext = find_containing_angular_element(em_jel, Set(dep_nodes))
    src_e_ext == 0 && error("angular plan: no source element contains the dependency nodes (iel=$iel, jel=$jel, d=$d)")
    em_iel = mesh.extra_mesh[iel]
    θd = 0.0; ϕd = 0.0
    for e = 1:em_iel.extra_nelem
        n = em_iel.extra_nop[e] + 1
        for a = 1:n, b = 1:n
            if em_iel.extra_connijk[e,a,b] == d
                θd = em_iel.extra_coords[1, d]; ϕd = em_iel.extra_coords[2, d]
            end
        end
    end
    nsrc = em_jel.extra_nop[src_e_ext] + 1
    θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds(em_jel.extra_coords, em_jel.extra_connijk, em_jel.extra_nop, src_e_ext)
    near_seam = isapprox(ϕmax, 2π; atol=1e-6)
    ϕd_fix = (near_seam && ϕd < 1e-8) ? 2π : ϕd
    θ_src = Vector{Float64}(undef, nsrc)
    for a = 1:nsrc; θ_src[a] = em_jel.extra_coords[1, em_jel.extra_connijk[src_e_ext,a,1]]; end
    ϕ_src = Vector{Float64}(undef, nsrc)
    for b = 1:nsrc
        raw = em_jel.extra_coords[2, em_jel.extra_connijk[src_e_ext,1,b]]
        ϕ_src[b] = (near_seam && raw < 1e-8) ? 2π : raw
    end
    Lθ = build_1d_lagrange_matrix(θ_src, [θd]); Lϕ = build_1d_lagrange_matrix(ϕ_src, [ϕd_fix])
    src_ips = Vector{Int}(undef, nsrc*nsrc); w = Vector{Float64}(undef, nsrc*nsrc)
    for a = 1:nsrc, b = 1:nsrc
        q = (a-1)*nsrc + b
        src_ips[q] = em_jel.extra_connijk[src_e_ext,a,b]
        w[q] = Lθ[1,a]*Lϕ[1,b]
    end
    inner = face_correspondence[(iel,fid)]
    irange, jrange, krange = face_ip_range(fid, ngl)
    lids = Int[]; nidx = Int[]
    for k in krange, j in jrange, i in irange
        c = inner[(i,j,k)]
        push!(lids, lidx(i,j,k,ngl)); push!(nidx, lidx(c[1],c[2],c[3],ngl))
    end
    ANG_STATS[1] += 1
    return AngularSamePlan(jel, lids, nidx, src_ips, w)
end

function eval_angular_same_plan!(out::Vector{Float64}, pl::AngularSamePlan, x_dg)
    n2 = length(pl.src_ips)
    vecs = Vector{Vector{Float64}}(undef, n2)
    for q = 1:n2; vecs[q] = x_dg[(pl.jel, pl.src_ips[q])]; end      # 49 lookups, once per call
    w = pl.w
    for n in eachindex(pl.lids)
        l = pl.nidx[n]; val = 0.0
        for q = 1:n2
            val += w[q]*vecs[q][l]
        end
        out[pl.lids[n]] = val
    end
    ANG_STATS[2] += 1
    return out
end

function make_angular_remote_plan(mesh, iel::Int, fid::Int, d::Int, remote_gel, remote_nodes, remote_angular_mesh, ngl::Int)
    em_remote = remote_angular_mesh[remote_gel]
    src_e_ext = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
    src_e_ext == 0 && return nothing                       # the original returns "not ready" in this case
    em_iel = mesh.extra_mesh[iel]
    θd = em_iel.extra_coords[1,d]; ϕd = em_iel.extra_coords[2,d]
    nsrc = em_remote.extra_nop[src_e_ext] + 1
    θmin,θmax,ϕmin,ϕmax = get_angular_element_bounds_ragged(em_remote.extra_coords, em_remote.extra_connijk, em_remote.extra_nop, src_e_ext)
    near_seam = isapprox(ϕmax, 2π; atol=1e-6)
    ϕd_fix = (near_seam && ϕd < 1e-8) ? 2π : ϕd
    θ_src = [em_remote.extra_coords[1, em_remote.extra_connijk[src_e_ext][a,1]] for a=1:nsrc]
    ϕ_src = Float64[]
    for b = 1:nsrc
        raw = em_remote.extra_coords[2, em_remote.extra_connijk[src_e_ext][1,b]]
        push!(ϕ_src, (near_seam && raw < 1e-8) ? 2π : raw)
    end
    Lθ = build_1d_lagrange_matrix(θ_src, [θd]); Lϕ = build_1d_lagrange_matrix(ϕ_src, [ϕd_fix])
    src_ips = Vector{Int}(undef, nsrc*nsrc); w = Vector{Float64}(undef, nsrc*nsrc)
    for a = 1:nsrc, b = 1:nsrc
        q = (a-1)*nsrc + b
        src_ips[q] = em_remote.extra_connijk[src_e_ext][a,b]
        w[q] = Lθ[1,a]*Lϕ[1,b]
    end
    irange, jrange, krange = face_ip_range(fid, ngl)
    lids = Int[]; gids = Int[]
    for k in krange, j in jrange, i in irange
        push!(lids, lidx(i,j,k,ngl)); push!(gids, Int(mesh.ip2gip[mesh.connijk[iel,i,j,k]]))
    end
    ANG_STATS[3] += 1
    return AngularRemotePlan(TInt(remote_gel), lids, gids, src_ips, w, Matrix{Int32}(undef, 0, 0), nothing)
end

function remote_plan_slots!(pl::AngularRemotePlan, t::GidTable)
    n2 = length(pl.src_ips); nf = length(pl.lids)
    S = Matrix{Int32}(undef, n2, nf)
    for n = 1:nf, q = 1:n2
        S[q,n] = get(t.slot, (pl.remote_gel, pl.src_ips[q], pl.gids[n]), Int32(0))
    end
    pl.slots = S; pl.slot_owner = t.slot
    ANG_STATS[5] += 1
    return S
end

function eval_angular_remote_plan!(out::Vector{Float64}, pl::AngularRemotePlan, gr)
    n2 = length(pl.src_ips); nf = length(pl.lids); w = pl.w
    if gr isa GidTable
        pl.slot_owner === gr.slot || remote_plan_slots!(pl, gr)
        S = pl.slots; vals = gr.vals
        for n = 1:nf
            val = 0.0
            for q = 1:n2
                s = S[q,n]
                v = s == 0 ? NaN : vals[s]                 # missing key or not yet arrived -> NaN, as in the original
                val += w[q]*v
            end
            out[pl.lids[n]] = val
        end
    else                                                    # plain Dict (round 1, old-path verification)
        for n = 1:nf
            val = 0.0; g = pl.gids[n]
            for q = 1:n2
                v = get(gr, (pl.remote_gel, pl.src_ips[q], g), NaN)
                val += w[q]*v
            end
            out[pl.lids[n]] = val
        end
    end
    ANG_STATS[4] += 1
    return out
end

function angular_same_plan!(mesh, angular_deps, iel::Int, fid::Int, jel::Int, d::Int, dep_nodes, face_correspondence, ngl::Int)
    dct = get!(() -> Dict{Tuple{Int,Int,Int},AngularSamePlan}(), ANGULAR_SAME_PLANS, angular_deps)
    return get!(dct, (iel,fid,d)) do
        make_angular_same_plan(mesh, iel, fid, jel, d, dep_nodes, face_correspondence, ngl)
    end
end
function angular_remote_plan!(mesh, angular_deps_cr, iel::Int, fid::Int, d::Int, remote_gel, remote_nodes, remote_angular_mesh, ngl::Int)
    dct = get!(() -> Dict{Tuple{Int,Int,Int},Union{Nothing,AngularRemotePlan}}(), ANGULAR_REMOTE_PLANS, angular_deps_cr)
    return get!(dct, (iel,fid,d)) do
        make_angular_remote_plan(mesh, iel, fid, d, remote_gel, remote_nodes, remote_angular_mesh, ngl)
    end
end

# ---------------- verification: plans vs your original builders, on the real final solution ----------------
function verify_angular_plans(setup, I_cg_all, comm; max_items::Int=300)
    mesh = setup.mesh; ngl = mesh.ngl
    gr_tab  = exchange_ncf_gid_table!(setup.gid_exch, I_cg_all, ngl)          # collective
    gr_dict = exchange_ncf_gid_values!(setup.gid_exch, I_cg_all, ngl)         # collective, plain Dict path
    empty!(ANGULAR_SAME_PLANS); empty!(ANGULAR_REMOTE_PLANS)
    n_same = 0; n_tab = 0; n_dict = 0; bad = 0; maxd = 0.0

    function compare(clo, out, irange, jrange, krange)
        nb = 0; md = 0.0
        for k in krange, j in jrange, i in irange
            x = clo(i,j,k); y = out[lidx(i,j,k,ngl)]
            (isnan(x) && isnan(y)) && continue
            reinterpret(UInt64, x) == reinterpret(UInt64, y) || (nb += 1)
            dd = abs(x - y); isnan(dd) || (md = max(md, dd))
        end
        return nb, md
    end

    cnt = 0
    for ((iel, d), deps) in setup.angular_deps
        cnt >= max_items && break
        for (dep_jel, dep_nodes) in deps
            length(dep_nodes) > 1 || continue
            jel = Int(dep_jel)
            for fid = 1:6
                cnt >= max_items && break
                Int(setup.neighbor_elem[iel][fid]) == jel || continue
                all(dn -> haskey(I_cg_all, (jel, dn)), dep_nodes) || continue
                em_jel = mesh.extra_mesh[jel]
                src_e = find_containing_angular_element(em_jel, Set(dep_nodes))
                src_e == 0 && continue
                clo = build_angular_ncf_upwind_unified(iel, fid, jel, d, src_e, mesh, em_jel, I_cg_all, setup.face_correspondence, ngl)
                pl = angular_same_plan!(mesh, setup.angular_deps, iel, fid, jel, d, dep_nodes, setup.face_correspondence, ngl)
                out = zeros(ngl^3); eval_angular_same_plan!(out, pl, I_cg_all)
                ir, jr, kr = face_ip_range(fid, ngl)
                nb, md = compare(clo, out, ir, jr, kr)
                n_same += 1; bad += nb; maxd = max(maxd, md); cnt += 1
            end
        end
    end

    cnt = 0
    for ((iel, fid, d), (remote_gel, remote_nodes)) in setup.angular_deps_cr
        cnt >= max_items && break
        length(remote_nodes) > 1 || continue
        haskey(setup.remote_angular_mesh, remote_gel) || continue
        em_remote = setup.remote_angular_mesh[remote_gel]
        src_e = find_containing_angular_element_ragged(em_remote, Set(remote_nodes))
        src_e == 0 && continue
        pl = angular_remote_plan!(mesh, setup.angular_deps_cr, iel, fid, d, remote_gel, remote_nodes, setup.remote_angular_mesh, ngl)
        pl === nothing && continue
        ir, jr, kr = face_ip_range(fid, ngl)
        for (which, gr) in ((1, gr_tab), (2, gr_dict))
            clo = build_angular_ncf_upwind_unified_remote(mesh, iel, remote_gel, d, src_e, em_remote, gr)
            out = zeros(ngl^3); eval_angular_remote_plan!(out, pl, gr)
            nb, md = compare(clo, out, ir, jr, kr)
            which == 1 ? (n_tab += 1) : (n_dict += 1)
            bad += nb; maxd = max(maxd, md)
        end
        cnt += 1
    end

    tot = MPI.Allreduce(Float64[n_same, n_tab, n_dict, bad], +, comm); mx = MPI.Allreduce(maxd, MPI.MAX, comm)
    MPI.Comm_rank(comm) == 0 || return String[]
    mxs = @sprintf("%.2e", mx)
    return String["angular plan verification: same-rank $(Int(tot[1])) faces, cross-rank $(Int(tot[2])) faces via table and " *
                  "$(Int(tot[3])) via Dict, $(Int(tot[4])) bitwise mismatches, max |Δ| $mxs (expected 0; a category with 0 faces verified nothing)"]
end
