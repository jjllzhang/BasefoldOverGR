#!/usr/bin/env bash
set -uo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNNER="$ROOT_DIR/scripts/run_release_c4_lambda128.sh"

if [[ ! -f "$RUNNER" ]]; then
  echo "Runner script not found: $RUNNER" >&2
  exit 1
fi

TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
RUN_ID="${RUN_ID:-${TIMESTAMP}_pid$$}"
BENCH_THREADS="${BENCH_THREADS:-8}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-1}"

declare -a CASES=(
  "results-new/8-thread-gr2p32-64-ring-switch compiler_eval_ring_switch ring-gr-2p32-64-ext 6 6 28"
  "results-new/8-thread-gr2p32-64-frobenius compiler_eval_frobenius ring-gr-2p32-64-ext 6 6 28"
  "results-new/8-thread-gr2p32-128-ring-switch compiler_eval_ring_switch ring-gr-2p32-128-ext 7 7 27"
  "results-new/8-thread-gr2p32-128-frobenius compiler_eval_frobenius ring-gr-2p32-128-ext 7 7 27"
)

SUMMARY_TSV="${SUMMARY_TSV:-$ROOT_DIR/results-new/gr2p32_compiler_eval_batch_${RUN_ID}.tsv}"
mkdir -p "$(dirname "$SUMMARY_TSV")"
printf 'out_dir\trun_suite\tcontext\tkappa\tell_min\tell_max\tstatus\n' > "$SUMMARY_TSV"

failures=0

run_case() {
  local out_dir="$1"
  local run_suite="$2"
  local context="$3"
  local kappa="$4"
  local ell_min="$5"
  local ell_max="$6"

  echo "==> OUT_DIR=$out_dir RUN_SUITE=$run_suite CONTEXTS=$context KAPPA=$kappa ELL=[$ell_min,$ell_max]"
  if (
    cd "$ROOT_DIR" && \
    OUT_DIR="$out_dir" \
    RUN_ID="$RUN_ID" \
    BENCH_THREADS="$BENCH_THREADS" \
    CONTINUE_ON_ERROR="$CONTINUE_ON_ERROR" \
    RUN_SUITE="$run_suite" \
    CONTEXTS="$context" \
    COMPILER_KAPPA="$kappa" \
    COMPILER_ELL_MIN="$ell_min" \
    COMPILER_ELL_MAX="$ell_max" \
    bash "$RUNNER"
  ); then
    printf '%s\t%s\t%s\t%s\t%s\t%s\tok\n' \
      "$out_dir" "$run_suite" "$context" "$kappa" "$ell_min" "$ell_max" >> "$SUMMARY_TSV"
    return 0
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\tfailed\n' \
    "$out_dir" "$run_suite" "$context" "$kappa" "$ell_min" "$ell_max" >> "$SUMMARY_TSV"
  return 1
}

for cfg in "${CASES[@]}"; do
  read -r out_dir run_suite context kappa ell_min ell_max <<< "$cfg"
  if ! run_case "$out_dir" "$run_suite" "$context" "$kappa" "$ell_min" "$ell_max"; then
    failures=$((failures + 1))
    echo "Failed: $out_dir" >&2
  fi
done

echo "Summary TSV: $SUMMARY_TSV"

if (( failures > 0 )); then
  echo "$failures case(s) failed." >&2
  exit 1
fi

echo "All 4 cases completed successfully."
