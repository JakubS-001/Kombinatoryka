#!/usr/bin/env bash
set -euo pipefail

# Build
make -C Generator
make -C AlgorytmDokladny compute_eigs

EDGE_FILE=tests/example_edges.txt
JSON_FILE=tests/example.json
TARGET=tests/target.txt
COMPUTED=tests/computed.txt
RESULTS=tests/results.json
EPS=1e-6

# convert
Generator/graph_to_json "$EDGE_FILE" > "$JSON_FILE"
# compute eigenvalues
AlgorytmDokladny/compute_eigs "$JSON_FILE" > "$COMPUTED"
# compare
python3 tests/compare.py "$TARGET" "$COMPUTED" "$RESULTS" "$EPS"

echo "Result written to $RESULTS"
cat "$RESULTS"
