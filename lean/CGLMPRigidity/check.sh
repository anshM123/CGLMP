#!/usr/bin/env bash
# Check the CGLMPRigidity files. Run from the `lean/` directory of this repository with the Lean
# toolchain on PATH, after `lake exe cache get`, e.g.
#   export PATH="$HOME/.elan/bin:$PATH"; lake exe cache get; bash CGLMPRigidity/check.sh
# Each file is elaborated by `lake env lean`; its .olean goes to the (git-ignored) build directory
# so that the next file can import it.
set -euo pipefail
out=.lake/build/lib/lean/CGLMPRigidity
mkdir -p "$out"
for f in Rigidity Links Model Bridge; do
  echo "== CGLMPRigidity/$f.lean"
  lake env lean -o "$out/$f.olean" -i "$out/$f.ilean" "CGLMPRigidity/$f.lean"
done
echo "== CGLMPRigidity/Axioms.lean"
lake env lean CGLMPRigidity/Axioms.lean
echo "== forbidden keywords (sorry/admit/axiom/native_decide):"
grep -nE '\b(sorry|admit|native_decide)\b|^\s*axiom\b' CGLMPRigidity/*.lean || echo "none"
