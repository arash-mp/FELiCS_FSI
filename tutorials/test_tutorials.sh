#!/usr/bin/env bash
# Runs all FELiCS tutorials referenced in documentation/Tutorials/index.md
# and prints a pass/fail summary.
#
# Logs for each run are written to test_tutorials_output/<timestamp>/ next to this script.
#
# Usage: ./test_tutorials [--timeout SECONDS]
#   Default timeout: 600 s (10 min) per tutorial.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
MAIN_PY="$REPO_ROOT/src/main.py"
TIMEOUT=600

while [[ $# -gt 0 ]]; do
    case "$1" in
        --timeout) TIMEOUT="$2"; shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 1 ;;
    esac
done

OUTPUT_DIR="$SCRIPT_DIR/test_tutorials_output/$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

pass=()
fail=()
fail_logs=()

run_tutorial() {
    local name="$1"
    local dir="$2"
    shift 2   # remaining args are the command + its arguments

    local slug="${name// /_}"
    local log="$OUTPUT_DIR/${slug}.log"

    printf "Running %-42s ... " "$name"

    (cd "$SCRIPT_DIR/$dir" && timeout "$TIMEOUT" "$@") 2>&1 \
        | sed 's/\x1b\[[0-9;]*m//g' > "$log"
    local exit_code="${PIPESTATUS[0]}"

    if [[ $exit_code -eq 0 ]]; then
        pass+=("$name")
        echo "PASSED"
    else
        fail+=("$name")
        fail_logs+=("$log")
        echo "FAILED (exit $exit_code)"
    fi
}

echo "============================================================"
echo "  FELiCS Tutorial Tests"
echo "  Logs: $OUTPUT_DIR"
echo "============================================================"
echo ""

run_tutorial "T1:_Cylinder_Wake_mesh"  "cylinder_wake_tutorial"  python cylinder_wake_mesh.py
run_tutorial "T1:_Cylinder_Wake"       "cylinder_wake_tutorial"  python solveBaseFlow.py
run_tutorial "T2:_Modal_Analysis"      "modal_analysis_tutorial" python "$MAIN_PY" -f modal.json -p -d
run_tutorial "T3:_Resolvent_geo_mesh"  "resolvent_tutorial"      python geoMesh.py
run_tutorial "T3:_Resolvent_mean_flow" "resolvent_tutorial"      python meanFlow.py
run_tutorial "T3:_Resolvent_Analysis"  "resolvent_tutorial"      python "$MAIN_PY" -f resolvent_settings_only_one_omega.json -p -d
run_tutorial "T4:_Input-Output"        "input_ouput_tutorial"    python "$MAIN_PY" -f input_output.json -p -d

echo ""
echo "============================================================"
echo "  Summary"
echo "============================================================"
echo ""

if [[ ${#pass[@]} -gt 0 ]]; then
    for t in "${pass[@]}"; do
        echo "  PASS  $t"
    done
fi

if [[ ${#fail[@]} -gt 0 ]]; then
    for i in "${!fail[@]}"; do
        echo "  FAIL  ${fail[$i]}"
        echo "        Last 15 lines of output (full log: ${fail_logs[$i]}):"
        tail -15 "${fail_logs[$i]}" | sed 's/^/          /'
        echo ""
    done
    exit 1
fi

echo ""
echo "All tutorials passed."
