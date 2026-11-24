#!/usr/bin/env bash
#
# Remove every Jexpresso mesh/SEM preprocess cache in the tree.
#
# Caches are per-case: each case keeps its own <case_dir>/.jexpresso_cache/
# directory (see _cache_dir in src/kernel/coupling/couplingStructs.jl), so
# there is no single directory to clear. They are also invalidated
# automatically — by the input fingerprint, by the St_mesh / St_metrics field
# signature, and by a shape check against the mesh actually in hand — so this
# script is a convenience, not part of the normal workflow.
#
# Run it from anywhere; it works on the repository this file lives in.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

found=0
while IFS= read -r -d '' dir; do
    echo " # removing $dir"
    rm -rf "$dir"
    found=$((found + 1))
done < <(find "$root" -type d -name '.jexpresso_cache' -print0)

# Legacy layout: caches written next to the .msh before they became per-case.
while IFS= read -r -d '' f; do
    echo " # removing $f"
    rm -f "$f"
    found=$((found + 1))
done < <(find "$root" \( -name 'MESH_*.jld2' -o -name 'PREPROCESS_*.jld2' \) -type f -print0)

if [ "$found" -eq 0 ]; then
    echo " # no Jexpresso caches found under $root"
fi
