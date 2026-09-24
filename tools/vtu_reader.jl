#---------------------------------------------------------------------------------
# Minimal reader for the VTK unstructured-grid files that Jexpresso writes.
#
# Jexpresso's only 2D/3D field output is VTK (`:outformat => "vtk"`), written by
# `write_vtk` in src/io/write_output.jl through WriteVTK's `pvtk_grid`. Nothing
# in the package reads those files back, so any post-processing that renders a
# figure from a finished run needs a reader. Adding ReadVTK.jl as a dependency
# for one plotting script is not worth the resolver churn, and the subset of the
# format that Jexpresso actually emits is small:
#
#   * `.pvtu` master file listing one `.vtu` piece per MPI rank,
#   * `format="appended"` `DataArray`s with `encoding="raw"`,
#   * `header_type="UInt64"` block sizes,
#   * no compressor (`write_vtk` passes `compress=false`).
#
# `format="ascii"` and inline base64 (`format="binary"`) are handled too, since
# they cost a few lines and make the reader usable on files written by other
# tools. A compressed appended blob is rejected with an explicit message rather
# than mis-decoded.
#
# Usage:
#
#     include("tools/vtu_reader.jl")
#     grid = read_vtk_grid("output/.../iter_21.pvtu")
#     grid.x, grid.y          # nodal coordinates, one entry per point
#     grid.data["rho"]        # a point-data field
#     grid.time               # simulation time, or `nothing` if not recorded
#---------------------------------------------------------------------------------

module VTUReader

using Base64

export read_vtk_grid, VTKGrid

"""
    VTKGrid

Nodal data of one output time: coordinates `x`, `y`, `z` (one entry per point),
the point-data fields in `data`, and the simulation `time` if the writer
recorded it as `TimeValue` field data (`nothing` for files written before that
was added).
"""
struct VTKGrid
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    data::Dict{String,Vector{Float64}}
    time::Union{Nothing,Float64}
end

#---------------------------------------------------------------------------------
# XML scanning
#
# The header of a VTU is plain ASCII up to `<AppendedData ...>_`; the raw binary
# blob follows the underscore. Read the whole file as bytes, cut it there, and
# scan the header with a cursor rather than a full XML parser — the tag set is
# fixed and shallow.
#---------------------------------------------------------------------------------

_attrs(tag::AbstractString) =
    Dict{String,String}(m.captures[1] => m.captures[2]
                        for m in eachmatch(r"([A-Za-z_][\w:.-]*)\s*=\s*\"([^\"]*)\"", tag))

function _split_header(raw::Vector{UInt8})
    marker = codeunits("<AppendedData")
    mi     = _findbytes(raw, marker)
    mi === nothing && return String(copy(raw)), nothing

    gt = findnext(==(UInt8('>')), raw, mi)
    gt === nothing && error("malformed <AppendedData> tag")
    us = findnext(==(UInt8('_')), raw, gt)
    us === nothing && error("<AppendedData> has no '_' data marker")

    return String(copy(raw[1:mi-1])), us + 1   # blob starts one byte past '_'
end

function _findbytes(hay::AbstractVector{UInt8}, needle::AbstractVector{UInt8})
    n, m = length(hay), length(needle)
    m == 0 && return 1
    @inbounds for i = 1:(n - m + 1)
        ok = true
        for j = 1:m
            if hay[i+j-1] != needle[j]
                ok = false
                break
            end
        end
        ok && return i
    end
    return nothing
end

"One `<DataArray>` and the section (`Points`, `PointData`, ...) it sits in."
struct _Array
    section::String
    attrs::Dict{String,String}
    content::String
end

const _SECTIONS = ("Points", "PointData", "CellData", "FieldData", "Cells",
                   "PPoints", "PPointData", "PCellData", "PFieldData", "PCells")

function _scan_arrays(xml::AbstractString)
    out     = _Array[]
    section = ""
    pos     = firstindex(xml)

    while true
        i = findnext(==('<'), xml, pos)
        i === nothing && break
        j = findnext(==('>'), xml, i)
        j === nothing && break
        tag = xml[i:j]

        m = match(r"^</?([A-Za-z]+)", tag)
        name = m === nothing ? "" : m.captures[1]

        if name in _SECTIONS
            section = startswith(tag, "</") ? "" : name
            pos = nextind(xml, j)
        elseif name == "DataArray" || name == "PDataArray"
            attrs = _attrs(tag)
            content = ""
            if endswith(tag, "/>")
                pos = nextind(xml, j)
            else
                k = findnext("</DataArray>", xml, nextind(xml, j))
                if k === nothing
                    pos = nextind(xml, j)
                else
                    content = xml[nextind(xml, j):prevind(xml, first(k))]
                    pos = nextind(xml, last(k))
                end
            end
            push!(out, _Array(section, attrs, content))
        else
            pos = nextind(xml, j)
        end
    end
    return out
end

#---------------------------------------------------------------------------------
# Payload decoding
#---------------------------------------------------------------------------------

function _jltype(t::AbstractString)
    t == "Float64" && return Float64
    t == "Float32" && return Float32
    t == "Int64"   && return Int64
    t == "Int32"   && return Int32
    t == "UInt64"  && return UInt64
    t == "UInt32"  && return UInt32
    t == "Int8"    && return Int8
    t == "UInt8"   && return UInt8
    error("VTUReader: unsupported DataArray type \"$t\"")
end

"Decode a length-prefixed binary block: `header_type` byte count, then payload."
function _decode_block(bytes::Vector{UInt8}, start::Int, T::Type, HT::Type)
    hsz    = sizeof(HT)
    nbytes = Int(only(reinterpret(HT, bytes[start:start+hsz-1])))
    first_ = start + hsz
    last_  = first_ + nbytes - 1
    last_ <= length(bytes) ||
        error("VTUReader: appended block runs past the end of the file")
    return Float64.(reinterpret(T, bytes[first_:last_]))
end

function _values(a::_Array, raw::Vector{UInt8}, blob::Union{Nothing,Int}, HT::Type)
    fmt = get(a.attrs, "format", "appended")
    T   = _jltype(get(a.attrs, "type", "Float64"))

    if fmt == "appended"
        blob === nothing && error("VTUReader: appended DataArray but no <AppendedData> block")
        off = parse(Int, get(a.attrs, "offset", "0"))
        return _decode_block(raw, blob + off, T, HT)
    elseif fmt == "binary"
        return _decode_block(base64decode(strip(a.content)), 1, T, HT)
    elseif fmt == "ascii"
        return Float64[parse(Float64, tok) for tok in split(a.content)]
    end
    error("VTUReader: unsupported DataArray format \"$fmt\"")
end

#---------------------------------------------------------------------------------
# Public entry points
#---------------------------------------------------------------------------------

"""
    read_vtu(path) -> VTKGrid

Read a single serial `.vtu` piece.
"""
function read_vtu(path::AbstractString)
    raw          = read(path)
    xml, blob    = _split_header(raw)

    vtkfile = match(r"<VTKFile\b[^>]*>", xml)
    vtkfile === nothing && error("VTUReader: $path is not a VTKFile")
    fattrs = _attrs(vtkfile.match)
    haskey(fattrs, "compressor") &&
        error("VTUReader: $path is compressed with $(fattrs["compressor"]); " *
              "Jexpresso writes with compress=false, so this file did not come " *
              "from write_vtk")
    HT = _jltype(get(fattrs, "header_type", "UInt64"))

    arrays = _scan_arrays(xml)

    x = Float64[]; y = Float64[]; z = Float64[]
    data = Dict{String,Vector{Float64}}()
    time = nothing

    for a in arrays
        if a.section == "Points"
            v = _values(a, raw, blob, HT)
            nc = parse(Int, get(a.attrs, "NumberOfComponents", "3"))
            np = length(v) ÷ nc
            x = [v[(ip-1)*nc + 1] for ip = 1:np]
            y = nc >= 2 ? [v[(ip-1)*nc + 2] for ip = 1:np] : zeros(np)
            z = nc >= 3 ? [v[(ip-1)*nc + 3] for ip = 1:np] : zeros(np)
        elseif a.section == "PointData"
            name = get(a.attrs, "Name", "")
            isempty(name) && continue
            data[name] = _values(a, raw, blob, HT)
        elseif a.section == "FieldData"
            if get(a.attrs, "Name", "") == "TimeValue"
                v = _values(a, raw, blob, HT)
                isempty(v) || (time = v[1])
            end
        end
    end

    return VTKGrid(x, y, z, data, time)
end

"""
    read_pvtu(path) -> VTKGrid

Read a `.pvtu` master file and concatenate its pieces. Nodes shared between
ranks appear once per piece; that is harmless for nearest-neighbour plotting,
which is all this reader is used for.
"""
function read_pvtu(path::AbstractString)
    xml  = String(read(path))
    dir  = dirname(path)
    srcs = [m.captures[1] for m in eachmatch(r"<Piece\b[^>]*Source\s*=\s*\"([^\"]*)\"", xml)]
    isempty(srcs) && error("VTUReader: $path lists no <Piece Source=...>")

    pieces = [read_vtu(isabspath(s) ? s : joinpath(dir, s)) for s in srcs]

    x = vcat((p.x for p in pieces)...)
    y = vcat((p.y for p in pieces)...)
    z = vcat((p.z for p in pieces)...)

    names = intersect((Set(keys(p.data)) for p in pieces)...)
    data  = Dict{String,Vector{Float64}}(n => vcat((p.data[n] for p in pieces)...)
                                         for n in names)

    time = nothing
    for p in pieces
        if p.time !== nothing
            time = p.time
            break
        end
    end

    return VTKGrid(x, y, z, data, time)
end

"""
    read_vtk_grid(path) -> VTKGrid

Read either a `.pvtu` master file or a serial `.vtu` piece.
"""
read_vtk_grid(path::AbstractString) =
    endswith(lowercase(path), ".pvtu") ? read_pvtu(path) : read_vtu(path)

end # module

using .VTUReader
