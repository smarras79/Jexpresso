#!/bin/bash
# SessionStart hook: install Julia + Jexpresso dependencies so that tests
# and lint commands can run during Claude Code on the web sessions.
#
# Idempotent: safe to re-run. The container caches state after the first
# successful run, so subsequent sessions skip the heavy installs.

set -euo pipefail

# Only run inside Claude Code on the web; locally, do nothing.
if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

REPO="${CLAUDE_PROJECT_DIR:-$(pwd)}"
JULIA_CHANNEL="1.11"

log() { printf '[session-start] %s\n' "$*"; }

# ---------------------------------------------------------------------------
# 1. Install Julia via juliaup (if not present)
# ---------------------------------------------------------------------------
JULIAUP_BIN="$HOME/.juliaup/bin"
export PATH="$JULIAUP_BIN:$PATH"

if ! command -v julia >/dev/null 2>&1; then
  log "Installing juliaup + Julia $JULIA_CHANNEL ..."
  curl -fsSL https://install.julialang.org \
    | sh -s -- --yes --default-channel "$JULIA_CHANNEL" --add-to-path=no
else
  log "Julia already available: $(julia --version)"
fi

# Ensure juliaup tracks 1.11
if command -v juliaup >/dev/null 2>&1; then
  juliaup add "$JULIA_CHANNEL" >/dev/null 2>&1 || true
  juliaup default "$JULIA_CHANNEL" >/dev/null 2>&1 || true
fi

# Persist PATH for the rest of the session so subsequent Bash tool calls
# (and any tests/linters Claude runs) see julia on PATH.
if [ -n "${CLAUDE_ENV_FILE:-}" ]; then
  echo "export PATH=\"$JULIAUP_BIN:\$PATH\"" >> "$CLAUDE_ENV_FILE"
fi

log "julia: $(command -v julia)"
julia --version

# ---------------------------------------------------------------------------
# 2. Instantiate the Jexpresso project (resolves Manifest.toml)
# ---------------------------------------------------------------------------
cd "$REPO"

log "Running Pkg.instantiate() against pinned Manifest.toml ..."
julia --project=. --color=no -e '
using Pkg
Pkg.instantiate()
Pkg.precompile()
println("Pkg.status():")
Pkg.status()
'

log "Session-start hook complete."
