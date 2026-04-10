#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUNNER="$ROOT_DIR/scripts/run_release_c4_lambda128.sh"

if [[ ! -x "$RUNNER" ]]; then
  echo "Runner script not found or not executable: $RUNNER" >&2
  exit 1
fi

for cmd in bash lscpu python3 taskset; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Missing required command: $cmd" >&2
    exit 1
  fi
done

TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
RUN_ID="${RUN_ID:-${TIMESTAMP}_pid$$}"
OUT_ROOT="${OUT_ROOT:-$ROOT_DIR/results}"
OUT_DIR="${OUT_DIR:-$OUT_ROOT/backend_eval_single_thread_7ctx_${RUN_ID}}"

D_MIN="${D_MIN:-4}"
D_MAX="${D_MAX:-22}"
K0="${K0:-1}"
C="${C:-4}"
LAMBDA="${LAMBDA:-128}"
COMMIT_WARMUP="${COMMIT_WARMUP:-1}"
COMMIT_REPS="${COMMIT_REPS:-3}"
EVAL_WARMUP="${EVAL_WARMUP:-1}"
EVAL_REPS="${EVAL_REPS:-3}"
SEED="${SEED:-0}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-1}"
CMD_TIMEOUT_SEC="${CMD_TIMEOUT_SEC:-0}"
PIN_BUILD="${PIN_BUILD:-1}"
BUILD_PARALLEL="${BUILD_PARALLEL:-1}"

CONTEXTS_CSV="${CONTEXTS_CSV:-ring-gr-2p16-162,ring-gr-2p16-64-ext,ring-gr-2p16-128-ext,field-prime64-ext,field-f2p64-ext,field-prime128-ext,field-f2p128-ext}"
# Default CPUs are first hardware threads of seven distinct physical cores,
# spread across NUMA nodes on the current 256-logical-CPU machine.
CPU_IDS_CSV="${CPU_IDS_CSV:-0,16,32,48,64,80,96}"

IFS=',' read -r -a CONTEXT_IDS <<< "$CONTEXTS_CSV"
IFS=',' read -r -a CPU_IDS <<< "$CPU_IDS_CSV"

trim_array_in_place() {
  local -n arr_ref="$1"
  local i
  for i in "${!arr_ref[@]}"; do
    arr_ref[$i]="${arr_ref[$i]//[[:space:]]/}"
  done
}

trim_array_in_place CONTEXT_IDS
trim_array_in_place CPU_IDS

if (( ${#CONTEXT_IDS[@]} == 0 )); then
  echo "No contexts configured." >&2
  exit 1
fi

if (( ${#CONTEXT_IDS[@]} != ${#CPU_IDS[@]} )); then
  echo "CONTEXTS_CSV count (${#CONTEXT_IDS[@]}) must match CPU_IDS_CSV count (${#CPU_IDS[@]})." >&2
  exit 1
fi

declare -A CPU_TO_CORE_KEY=()
declare -A CPU_TO_NODE=()
while IFS=',' read -r cpu core socket node; do
  [[ "$cpu" =~ ^# ]] && continue
  CPU_TO_CORE_KEY["$cpu"]="${core}-${socket}"
  CPU_TO_NODE["$cpu"]="$node"
done < <(lscpu -p=CPU,CORE,SOCKET,NODE)

declare -A SEEN_CORE_KEY=()
for cpu in "${CPU_IDS[@]}"; do
  if ! [[ "$cpu" =~ ^[0-9]+$ ]]; then
    echo "Invalid CPU id: $cpu" >&2
    exit 1
  fi
  if [[ -z "${CPU_TO_CORE_KEY[$cpu]:-}" ]]; then
    echo "CPU $cpu is not present in lscpu topology output." >&2
    exit 1
  fi
  if ! taskset -c "$cpu" true >/dev/null 2>&1; then
    echo "CPU $cpu is not available to the current shell. Check /proc/self/status." >&2
    exit 1
  fi
  core_key="${CPU_TO_CORE_KEY[$cpu]}"
  if [[ -n "${SEEN_CORE_KEY[$core_key]:-}" ]]; then
    echo "CPU $cpu shares a physical core with CPU ${SEEN_CORE_KEY[$core_key]}." >&2
    echo "Choose one hardware thread per physical core." >&2
    exit 1
  fi
  SEEN_CORE_KEY["$core_key"]="$cpu"
done

mkdir -p "$OUT_DIR"

SUMMARY_TSV="$OUT_DIR/context_cpu_map.tsv"
{
  echo -e "context_id\tcpu_id\tnuma_node\tout_dir\tlog_path"
} > "$SUMMARY_TSV"

echo "Run root: $OUT_DIR"
echo "Cpus_allowed_list: $(grep Cpus_allowed_list /proc/self/status | awk '{print $2}')"
echo "Context to CPU mapping:"

declare -a PIDS=()
declare -a CHILD_CONTEXTS=()
declare -a CHILD_OUT_DIRS=()
declare -a CHILD_LOGS=()

terminate_children() {
  local signal_name="$1"
  trap - INT TERM
  echo "Received $signal_name, terminating child jobs..." >&2
  local pid
  for pid in "${PIDS[@]:-}"; do
    kill "$pid" 2>/dev/null || true
  done
  wait || true
  exit 130
}

trap 'terminate_children INT' INT
trap 'terminate_children TERM' TERM

for idx in "${!CONTEXT_IDS[@]}"; do
  context_id="${CONTEXT_IDS[$idx]}"
  cpu_id="${CPU_IDS[$idx]}"
  node_id="${CPU_TO_NODE[$cpu_id]}"
  child_name="$(printf "%02d_%s" "$idx" "$context_id")"
  child_out_dir="$OUT_DIR/$child_name"
  child_log="$child_out_dir/run.log"
  child_run_id="${RUN_ID}_${child_name}"

  mkdir -p "$child_out_dir"
  echo "  $context_id -> CPU $cpu_id (NUMA node $node_id)"
  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$context_id" "$cpu_id" "$node_id" "$child_out_dir" "$child_log" >> "$SUMMARY_TSV"

  (
    cd "$ROOT_DIR"
    exec env \
      RUN_SUITE=basefold_release \
      RUN_ID="$child_run_id" \
      OUT_DIR="$child_out_dir" \
      D_MIN="$D_MIN" \
      D_MAX="$D_MAX" \
      K0="$K0" \
      C="$C" \
      LAMBDA="$LAMBDA" \
      CONTEXTS="$context_id" \
      COMMIT_WARMUP="$COMMIT_WARMUP" \
      COMMIT_REPS="$COMMIT_REPS" \
      EVAL_WARMUP="$EVAL_WARMUP" \
      EVAL_REPS="$EVAL_REPS" \
      SEED="$SEED" \
      CONTINUE_ON_ERROR="$CONTINUE_ON_ERROR" \
      CMD_TIMEOUT_SEC="$CMD_TIMEOUT_SEC" \
      BENCH_THREADS=1 \
      CPU_PIN_MODE=manual \
      CPU_SET="$cpu_id" \
      PIN_BUILD="$PIN_BUILD" \
      ISOLATE_BUILD_DIR=1 \
      CMAKE_BUILD_PARALLEL_LEVEL="$BUILD_PARALLEL" \
      bash "$RUNNER"
  ) >"$child_log" 2>&1 &

  PIDS+=("$!")
  CHILD_CONTEXTS+=("$context_id")
  CHILD_OUT_DIRS+=("$child_out_dir")
  CHILD_LOGS+=("$child_log")
done

failures=0
for idx in "${!PIDS[@]}"; do
  pid="${PIDS[$idx]}"
  context_id="${CHILD_CONTEXTS[$idx]}"
  if wait "$pid"; then
    echo "Completed: $context_id"
  else
    echo "Failed: $context_id (see ${CHILD_LOGS[$idx]})" >&2
    failures=$((failures + 1))
  fi
done

combined_csv="$OUT_DIR/backend_eval_results.csv"
header_written=0
rm -f "$combined_csv"
for child_out_dir in "${CHILD_OUT_DIRS[@]}"; do
  child_csv="$child_out_dir/backend_eval_results.csv"
  if [[ ! -f "$child_csv" ]]; then
    continue
  fi
  if (( header_written == 0 )); then
    cat "$child_csv" > "$combined_csv"
    header_written=1
  else
    tail -n +2 "$child_csv" >> "$combined_csv"
  fi
done

if (( header_written == 0 )); then
  echo "No backend_eval_results.csv files were produced." >&2
  exit 1
fi

echo "Combined CSV: $combined_csv"
echo "Mapping TSV:  $SUMMARY_TSV"

if (( failures > 0 )); then
  echo "$failures child job(s) failed." >&2
  exit 1
fi

