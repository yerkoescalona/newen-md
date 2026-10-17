#!/usr/bin/env bash
# Regenerate spec/GromosSpec/Extracted/ from crates/gromos-logic:
#   Charon  (Rust -> LLBC, a serialized MIR)  then  Aeneas (LLBC -> pure Lean functions).
# Needs the prebuilt Aeneas release (it bundles the matching Charon) in $AENEAS_DIR and the
# Rust nightly named in $AENEAS_DIR/rust-toolchain. See spec/README.md.
set -euo pipefail

AENEAS_DIR=${AENEAS_DIR:-$HOME/.local/aeneas/aeneas}
export PATH="$AENEAS_DIR:$PATH"
spec=$(cd "$(dirname "$0")/.." && pwd)
root=$(cd "$spec/.." && pwd)
llbc_dir="$root/target/charon"
out="$spec/GromosSpec/Extracted"

# The Aeneas binary and the Aeneas Lean library (lakefile.lean) must be the same commit.
# `aeneas -version` prints e.g. `aeneas nightly-2026.10.17-e62dd83` (a short commit hash).
bin_rev=$(aeneas -version | grep -oE '[0-9a-f]{7,}$')
lib_rev=$(grep -oE '"[0-9a-f]{40}"' "$spec/lakefile.lean" | tr -d '"')
if [ -z "$bin_rev" ] || [ "${lib_rev#"$bin_rev"}" = "$lib_rev" ]; then
  echo "aeneas binary is $bin_rev but spec/lakefile.lean pins $lib_rev" >&2
  exit 1
fi

mkdir -p "$llbc_dir" "$out"
(cd "$root/crates/gromos-logic" &&
  CARGO_TARGET_DIR="$llbc_dir" charon cargo --preset=aeneas --dest-file "$llbc_dir/gromos_logic.llbc")
aeneas -backend lean -use-lean-modules false -dest "$out" "$llbc_dir/gromos_logic.llbc" \
  >"$llbc_dir/aeneas.log" 2>&1 || { cat "$llbc_dir/aeneas.log" >&2; exit 1; }

for f in "$out"/*.lean; do
  { echo "-- Do not edit: generated from crates/gromos-logic by spec/scripts/extract.sh."
    cat "$f"; } > "$f.tmp" && mv "$f.tmp" "$f"
done
echo "extracted: $(ls "$out")"
