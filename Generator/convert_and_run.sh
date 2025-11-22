#!/usr/bin/env bash
# Skrypt: konwertuje prosty plik edge-list -> JSON i wywołuje program algorytmu
# Usage: ./convert_and_run.sh input_edges.txt target_spectrum.txt

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 input_edges.txt target_spectrum.txt" >&2
  exit 2
fi

INPUT="$1"
SPECTRUM="$2"

GENERATOR_CLI="$(dirname "$0")/graph_to_json"
ALGDOK_BIN="$(dirname "$0")/../AlgorytmDokladny/algdok"

if [ ! -x "$GENERATOR_CLI" ]; then
  echo "Generator not built. Compile it: gcc -o $GENERATOR_CLI $(basename $GENERATOR_CLI).c" >&2
  exit 3
fi
if [ ! -x "$ALGDOK_BIN" ]; then
  echo "AlgorytmDokladny not built. Run make in AlgorytmDokladny." >&2
  exit 4
fi

# Wygeneruj JSON do tymczasowego pliku
tmpjson=$(mktemp /tmp/graph_json.XXXXXX)
trap 'rm -f "$tmpjson"' EXIT

"$GENERATOR_CLI" "$INPUT" > "$tmpjson"

echo "Converted $INPUT -> $tmpjson. Calling AlgorytmDokladny..."

"$ALGDOK_BIN" "$tmpjson" "$SPECTRUM"
