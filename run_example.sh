#!/usr/bin/env bash
set -euo pipefail

# Build local cJSON + generator + algdok
make -C Generator
make -C AlgorytmDokladny

# Prepare files
EDGE_FILE=tests/example_edges.txt
JSON_FILE=tests/example.json
TARGET=tests/target.txt

# Convert edge list to JSON
Generator/graph_to_json "$EDGE_FILE" > "$JSON_FILE"

echo "JSON written to $JSON_FILE"

# Run algorithm
AlgorytmDokladny/algdok "$JSON_FILE" "$TARGET"
