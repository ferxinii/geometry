#!/usr/bin/env bash
# Compile robust_predicates with ROBUST_PREDICATES_PRINT_SIZE to emit stage_d
# buffer sizes as deprecation warnings, then parse them into a clean table.
set -euo pipefail

RPDIR="$(cd "$(dirname "$0")" && pwd)"
CXX="${CXX:-c++}"
INCLUDES="-I$RPDIR/includes/geometry -I$RPDIR/includes/mp11"
FLAGS="-std=c++20 -fsyntax-only -Wno-error -Wdeprecated-declarations -DROBUST_PREDICATES_PRINT_SIZE"
FILES=(
    "$RPDIR/predicates_basic_n2.cpp"
    "$RPDIR/predicates_basic_n3.cpp"
    "$RPDIR/predicates_basic_other.cpp"
    "$RPDIR/predicates_power_n1.cpp"
    "$RPDIR/predicates_power_n2.cpp"
    "$RPDIR/predicates_power_n3.cpp"
    "$RPDIR/predicates_lp_n2.cpp"
    "$RPDIR/predicates_lp_n3.cpp"
)

echo "Compiling to extract stage_d sizes (may take a minute)..."
echo ""

RAW=$(for f in "${FILES[@]}"; do
    "$CXX" $FLAGS $INCLUDES "$f" 2>&1 || true
done)

echo "$RAW" | awk '
/warning:.*show_stage_[db]_size<[0-9]+>/ {
    stage = ($0 ~ /show_stage_d_size/) ? "d" : "b"

    match($0, /show_stage_[db]_size<[0-9]+>/)
    chunk = substr($0, RSTART, RLENGTH)
    gsub(/[^0-9]/, "", chunk)
    size = chunk + 0

    getline srcline

    if (match(srcline, /_size_[A-Za-z0-9_]+/)) {
        name = substr(srcline, RSTART + 6, RLENGTH - 6)
        if (stage == "b") sub(/_b$/, "", name)
    } else {
        name = "?"
    }

    bytes = size * 8
    if (bytes < 1024)
        human = sprintf("%d B", bytes)
    else if (bytes < 1048576)
        human = sprintf("%.1f KB", bytes / 1024.0)
    else
        human = sprintf("%.1f MB", bytes / 1048576.0)
    printf "  %-32s  stage_%s  %9d doubles   %s\n", name, stage, size, human
}
'
