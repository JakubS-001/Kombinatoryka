#!/usr/bin/env bash
# Wrapper: use nauty's genbg to generate bipartite graphs and convert each to JSON
# Requires: genbg and showg from nauty, and Generator/graph_to_json (compiled)

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 n1 n2 [genbg-args...]" >&2
  echo "Example: $0 2 2 -q" >&2
  exit 2
fi

n1="$1"; shift
n2="$1"; shift

GENBG=$(command -v genbg || true)
SHOWG=$(command -v showg || true)
GRAPH_TO_JSON="$(dirname "$0")/graph_to_json"

if [ -z "$GENBG" ] || [ -z "$SHOWG" ]; then
  echo "genbg or showg not found in PATH. Please install nauty and ensure genbg/showg are available." >&2
  exit 3
fi
if [ ! -x "$GRAPH_TO_JSON" ]; then
  echo "Generator binary $GRAPH_TO_JSON not found. Build it with 'make' in Generator." >&2
  exit 4
fi

# genbg will output bipartite graphs (graph6 or sparse6). We'll pipe each graph through showg -e to get edge list.
# Assumption: genbg produces graphs with vertex ordering such that first n1 vertices are left part.

tempfile=$(mktemp /tmp/genbg_out.XXXXXX)
trap 'rm -f "$tempfile"' EXIT

echo "Running: genbg -b $n1 $n2 $@"
genbg -b "$n1" "$n2" "$@" > "$tempfile"

count=0
while IFS= read -r line; do
  # send single-graph line to showg to get edge-list (-e) and then convert
  echo "$line" | showg -e - | awk 'BEGIN{m=0; n=0} /n=/{n=$0} {print}' >/dev/null 2>&1 || true
  # showg -e prints edges; we need to convert to our edge-list format: n1 n2 m ; edges
  # We'll create a small pipeline: convert showg -e output into our input format
  jsonout="$(mktemp /tmp/graph_json.XXXXXX)"
  echo "$line" | showg -e - | awk -v n1="$n1" -v n2="$n2" '
    BEGIN{m=0}
    /^\(/ { # edge lines like (u,v)
      gsub("[(),]"," ",$0);
      split($0,a," "); u=a[2]; v=a[3]; edges[++m]=u " " v
    }
    END{
      if (m==0) { print n1 " " n2 " " 0; exit }
      print n1 " " n2 " " m;
      for(i=1;i<=m;i++) print edges[i];
    }' > "$jsonout"
  # convert edge-list to json using graph_to_json
  "$GRAPH_TO_JSON" "$jsonout" || true
  rm -f "$jsonout"
  count=$((count+1))
done < "$tempfile"

echo "Processed $count graphs. Note: this script assumes genbg/showg produce canonical vertex ordering where first $n1 vertices are left partition." 
