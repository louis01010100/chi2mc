#!/usr/bin/env bash
# End-to-end manual test: runs chi2mc on fixture.tsv and prints results.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Build release binary
echo "Building chi2mc (release)..."
cargo build --release --manifest-path "$PROJECT_ROOT/Cargo.toml" 2>&1

BINARY="$PROJECT_ROOT/target/release/chi2mc"
INPUT="$SCRIPT_DIR/fixture.tsv"
OUTPUT="$SCRIPT_DIR/output.tsv"

echo ""
echo "Running chi2mc on fixture.tsv..."
echo ""

"$BINARY" "$INPUT" "$OUTPUT" \
    --min-expected-count 5 \
    --n-simulations 10000 \
    --random-seed 42

echo ""
echo "--- Output TSV (first 20 lines) ---"
head -20 "$OUTPUT"

echo ""
echo "--- Summary Report ---"
cat "${OUTPUT%.tsv}_summary.txt"
