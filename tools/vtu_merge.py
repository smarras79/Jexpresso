#!/usr/bin/env python3
"""Merge the per-rank pieces of one Jexpresso snapshot into a single compact .vtu.

Usage: vtu_merge.py <output_dir> <iter> <out.vtu> [field ...]

Keeps only the point fields named (default: u v w θ), stores everything as
Float32 / Int32 and zlib-compresses each array, so a >1 GB 225-piece dump
comes out around 100 MB and opens in ParaView as one dataset. Points shared
between ranks are duplicated, which is harmless for visualisation.

Env: ZMAX=<m> drops every cell whose lowest point is above ZMAX.
"""
import sys, io, re, struct, glob, os, zlib
import numpy as np

def read_vtu(path):
    raw = io.open(path, 'rb').read()
    hdr_end = raw.find(b'<AppendedData')
    hdr = raw[:hdr_end].decode('utf-8', errors='replace')
    j = raw.find(b'_', hdr_end) + 1
    blob = raw[j:]
    out = {}
    for m in re.finditer(r'<DataArray type="(\w+)" Name="([^"]*)" NumberOfComponents="(\d+)" format="appended" offset="(\d+)"', hdr):
        typ, name, nc, off = m.group(1), m.group(2), int(m.group(3)), int(m.group(4))
        dt = {'Float64': '<f8', 'Float32': '<f4', 'Int64': '<i8', 'Int32': '<i4', 'UInt8': 'u1'}[typ]
        n = struct.unpack('<Q', blob[off:off+8])[0]
        a = np.frombuffer(blob[off+8:off+8+n], dtype=dt)
        out[name] = a.reshape(-1, nc) if nc > 1 else a
    return out

def main():
    outdir, it, outfile = sys.argv[1], sys.argv[2], sys.argv[3]
    fields = sys.argv[4:] or ['u', 'v', 'w', 'θ']
    zmax = float(os.environ.get('ZMAX', 'inf'))
    pieces = sorted(glob.glob(os.path.join(outdir, f'iter_{it}', f'iter_{it}_*.vtu')))
    if not pieces:
        raise SystemExit(f'no pieces under {outdir}/iter_{it}')

    P, C, O, T = [], [], [], []
    F = {f: [] for f in fields}
    npts = 0
    for k, p in enumerate(pieces):
        a = read_vtu(p)
        pts  = a['Points'].astype(np.float32)
        conn = a['connectivity'].astype(np.int64)
        offs = a['offsets'].astype(np.int64)
        typs = a['types'].astype(np.uint8)
        if np.isfinite(zmax):
            starts = np.concatenate([[0], offs[:-1]])
            keep = np.array([pts[conn[s:e], 2].min() <= zmax for s, e in zip(starts, offs)])
            if not keep.all():
                newconn, newoffs, acc = [], [], 0
                for s, e, kp in zip(starts, offs, keep):
                    if kp:
                        newconn.append(conn[s:e]); acc += e - s; newoffs.append(acc)
                conn = np.concatenate(newconn) if newconn else conn[:0]
                offs = np.array(newoffs, dtype=np.int64); typs = typs[keep]
        P.append(pts); C.append(conn + npts); O.append(offs + (O[-1][-1] if O else 0)); T.append(typs)
        for f in fields:
            if f not in a:
                raise SystemExit(f'{p}: no field {f!r}; have {sorted(a)}')
            F[f].append(a[f].astype(np.float32))
        npts += len(pts)
        if k % 25 == 0:
            print(f'  {k+1}/{len(pieces)} pieces, {npts} points', flush=True)

    pts  = np.concatenate(P); conn = np.concatenate(C).astype(np.int32)
    offs = np.concatenate(O).astype(np.int32); typs = np.concatenate(T)
    ncell = len(offs)

    # appended section, one zlib block per array, UInt64 headers
    app, offsets = io.BytesIO(), {}
    def put(name, arr):
        raw = arr.tobytes(); comp = zlib.compress(raw, 6)
        offsets[name] = app.tell()
        app.write(struct.pack('<QQQQ', 1, len(raw), len(raw), len(comp)))
        app.write(comp)
    put('Points', pts); put('connectivity', conn); put('offsets', offs); put('types', typs)
    for f in fields:
        put(f, np.concatenate(F[f]))

    def da(typ, name, nc, key):
        return f'        <DataArray type="{typ}" Name="{name}" NumberOfComponents="{nc}" format="appended" offset="{offsets[key]}"/>\n'
    with io.open(outfile, 'wb') as fh:
        fh.write(('<?xml version="1.0"?>\n'
                  '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" '
                  'header_type="UInt64" compressor="vtkZLibDataCompressor">\n'
                  '  <UnstructuredGrid>\n'
                  f'    <Piece NumberOfPoints="{len(pts)}" NumberOfCells="{ncell}">\n'
                  '      <Points>\n' + da('Float32', 'Points', 3, 'Points') + '      </Points>\n'
                  '      <Cells>\n' + da('Int32', 'connectivity', 1, 'connectivity')
                  + da('Int32', 'offsets', 1, 'offsets') + da('UInt8', 'types', 1, 'types')
                  + '      </Cells>\n      <PointData>\n'
                  + ''.join(da('Float32', f, 1, f) for f in fields)
                  + '      </PointData>\n    </Piece>\n  </UnstructuredGrid>\n'
                  '  <AppendedData encoding="raw">\n   _').encode('utf-8'))
        fh.write(app.getvalue())
        fh.write(b'\n  </AppendedData>\n</VTKFile>\n')
    print(f'wrote {outfile}: {len(pts)} points, {ncell} cells, {len(fields)} fields, '
          f'{os.path.getsize(outfile)/1e6:.1f} MB')

if __name__ == '__main__':
    main()
