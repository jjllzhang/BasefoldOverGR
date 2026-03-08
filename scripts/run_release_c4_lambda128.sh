#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_ROOT="${OUT_ROOT:-$ROOT_DIR/results}"
TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
RUN_ID="${RUN_ID:-${TIMESTAMP}_pid$$}"
ISOLATE_BUILD_DIR="${ISOLATE_BUILD_DIR:-0}"  # set 1 to use build-release-<RUN_ID> per run
DEFAULT_BUILD_DIR="${ROOT_DIR}/build-release"
if [[ -n "${BUILD_DIR:-}" ]]; then
  BUILD_DIR="$BUILD_DIR"
elif [[ "$ISOLATE_BUILD_DIR" == "1" ]]; then
  BUILD_DIR="${DEFAULT_BUILD_DIR}-${RUN_ID}"
else
  BUILD_DIR="$DEFAULT_BUILD_DIR"
fi
OUT_DIR="${OUT_DIR:-$OUT_ROOT/release_c4_lambda128_sweep_${RUN_ID}}"

poly_degree_from_coeff_list() {
  local coeffs="$1"
  awk -F',' '
    BEGIN { deg = -1 }
    {
      for (i = 1; i <= NF; i++) {
        gsub(/[[:space:]]/, "", $i)
        if ($i == "" || $i !~ /^-?[0-9]+$/) {
          print "ERR"
          exit
        }
        if (($i + 0) != 0) deg = i - 1
      }
    }
    END {
      if (deg < 0) print 0
      else print deg
    }' <<< "$coeffs"
}

is_power_of_two() {
  local value="$1"
  [[ "$value" =~ ^[1-9][0-9]*$ ]] || return 1
  (( (value & (value - 1)) == 0 ))
}

# Target profile: rate = 1/4 (c=4), security = 128 bits.
C="${C:-4}"
K0="${K0:-1}"
LAMBDA="${LAMBDA:-128}"
D_MIN="${D_MIN:-3}"
D_MAX="${D_MAX:-29}"

COMMIT_WARMUP="${COMMIT_WARMUP:-1}"
COMMIT_REPS="${COMMIT_REPS:-3}"
EVAL_WARMUP="${EVAL_WARMUP:-1}"
EVAL_REPS="${EVAL_REPS:-3}"
SEED="${SEED:-0}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-1}"
CMD_TIMEOUT_SEC="${CMD_TIMEOUT_SEC:-0}"
CONTEXTS="${CONTEXTS:-all}"  # all or comma list, see valid ids in parsing block
BENCH_THREADS="${BENCH_THREADS:-8}"  # set 0 to keep runtime defaults
CPU_PIN_MODE="${CPU_PIN_MODE:-none}"  # none|manual|slot
CPU_SET="${CPU_SET:-}"                # manual mode: e.g. 0-31 or 0,2,4-10
RUN_SLOT="${RUN_SLOT:-0}"             # slot mode: 0-based slot index
RUN_SLOTS_TOTAL="${RUN_SLOTS_TOTAL:-1}"  # slot mode: number of concurrent runs
USE_SMT_IN_SLOT="${USE_SMT_IN_SLOT:-0}"  # slot mode: 0=one hw thread/core, 1=all hw threads
PIN_BUILD="${PIN_BUILD:-0}"              # set 1 to also pin cmake configure/build
CHALLENGE_FIELD_EXT_DEG2="${CHALLENGE_FIELD_EXT_DEG2:-0,1;1;1}"  # E(U)=zeta+U+U^2
CHALLENGE_FIELD_EXT_DEG3="${CHALLENGE_FIELD_EXT_DEG3:-1;1;0;1}"  # E(U)=1+U+U^3
CHALLENGE_RING_EXT_DEG2="${CHALLENGE_RING_EXT_DEG2:-0,1;1;1}"    # E(U)=zeta+U+U^2
CHALLENGE_RING_EXT_DEG3="${CHALLENGE_RING_EXT_DEG3:-1;1;0;1}"    # E(U)=1+U+U^3
EXT128_BETA_TRACE1_DEFAULT='1,0,0,1,1,1,1,1,0,0,0,0,1,0,0,1,1,0,0,0,0,1,1,0,0,0,1,1,1,0,1,1,0,0,0,0,0,1,1,0,0,1,0,1,0,0,1,0,1,1,0,1,0,1,0,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,1,1,1,1,1,0,1,0,1,1,0,1,0,0,0,1,1,0,1,1,0,0,1,0,0,0,1,0,1,0,0,0,1,0,0,0,0,1,1,1,0,1,1,0,1,0,1,1,0,0,0,0,1,0,0,1,1,1'
CHALLENGE_F2_128_EXT_DEG2="${CHALLENGE_F2_128_EXT_DEG2:-$EXT128_BETA_TRACE1_DEFAULT;1;1}"  # E(U)=U^2+U+beta, Tr(beta)=1
CHALLENGE_GR_128_EXT_DEG2="${CHALLENGE_GR_128_EXT_DEG2:-$EXT128_BETA_TRACE1_DEFAULT;1;1}"  # lift of the same Artin-Schreier form
CHALLENGE_F2_64_EXT_DEG3="${CHALLENGE_F2_64_EXT_DEG3:-$CHALLENGE_FIELD_EXT_DEG3}"  # E(U)=1+U+U^3
CHALLENGE_F3_40_EXT_DEG3="${CHALLENGE_F3_40_EXT_DEG3:-0,1;1;0;1}"  # E(U)=zeta+U+U^3 over F_(3^40)
CHALLENGE_F3_81_EXT_DEG2="${CHALLENGE_F3_81_EXT_DEG2:-0,1;1;1}"    # E(U)=zeta+U+U^2 over F_(3^81)

FIELD255_MOD="${FIELD255_MOD:-57896044618658097711785492504343953926634992332820282019728792003956564819949}"  # 2^255 - 19
FIELD255_F="${FIELD255_F:-1,1}"       # x + 1
FIELD255_ZETA="${FIELD255_ZETA:-0,1}"  # x

F2_256_MOD="${F2_256_MOD:-2}"
F2_256_F_DEFAULT='1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,1,1,1,0,1,1,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1,0,1,1,1,0,0,1,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,1,1,0,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,1,1,0,1,1,0,1,0,1,1,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,0,0,1,1,1,0,1,0,0,0,1,1,1,0,1,1,0,0,1,0,1,1,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,1,0,0,1,0,1'
F2_256_F="${F2_256_F:-$F2_256_F_DEFAULT}"
F2_256_ZETA="${F2_256_ZETA:-0,1}"  # x

FIELD128_MOD="${FIELD128_MOD:-326594724262804054738278293730872375507}"  # 128-bit prime
FIELD128_F="${FIELD128_F:-1,1}"  # x + 1
FIELD128_ZETA="${FIELD128_ZETA:-0,1}"  # x

FIELD64_MOD="${FIELD64_MOD:-18446744073709551557}"  # 64-bit prime: 2^64 - 59
FIELD64_F="${FIELD64_F:-1,1}"  # x + 1
FIELD64_ZETA="${FIELD64_ZETA:-0,1}"  # x

# F_(3^40) context (64-bit class): degree-3 extension challenges by default.
F3_40_MOD="${F3_40_MOD:-3}"
F3_40_R="${F3_40_R:-40}"
F3_40_F_DEFAULT='1,2,0,1,2,0,0,2,1,0,0,1,2,1,2,0,0,2,1,1,0,0,2,0,0,2,2,2,2,1,1,0,1,0,2,0,2,0,1,0,1'  # GF(3) irreducible, deg=40
F3_40_F="${F3_40_F:-$F3_40_F_DEFAULT}"
F3_40_ZETA="${F3_40_ZETA:-0,1}"  # x

# F_(3^81) context (128-bit class): degree-2 extension challenges by default.
F3_81_MOD="${F3_81_MOD:-3}"
F3_81_R="${F3_81_R:-81}"
F3_81_F_DEFAULT='1,1,2,2,2,1,2,2,0,2,2,0,2,1,2,1,1,1,0,0,0,0,1,0,0,1,1,1,1,2,0,1,1,2,1,0,1,1,1,0,0,0,2,2,2,2,2,0,0,2,2,0,2,1,0,2,0,1,1,0,1,1,1,1,1,1,2,1,0,0,2,0,1,1,2,0,0,1,1,0,1,1'  # GF(3) irreducible, deg=81
F3_81_F="${F3_81_F:-$F3_81_F_DEFAULT}"
F3_81_ZETA="${F3_81_ZETA:-0,1}"  # x

RING_F_64_DEFAULT='1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,0,1,1,0,0,0,0,0,1,0,1'
# x^128 + x^7 + x^2 + x + 1 (AES-GCM polynomial), also used as its 2-adic lift.
RING_F_128_DEFAULT='1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1'
RING_F_64="${RING_F_64:-$RING_F_64_DEFAULT}"
RING_F_128="${RING_F_128:-$RING_F_128_DEFAULT}"

F2_128_MOD="${F2_128_MOD:-2}"
F2_128_F="${F2_128_F:-$RING_F_128_DEFAULT}"
F2_128_ZETA="${F2_128_ZETA:-0,1}"  # x

F2_64_MOD="${F2_64_MOD:-2}"
F2_64_F="${F2_64_F:-$RING_F_64_DEFAULT}"
F2_64_ZETA="${F2_64_ZETA:-0,1}"  # x

RING_F_162_DEFAULT='1,1,0,0,1,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,0,0,0,0,1,1,1,1,0,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,0,0,1,0,1,1,0,0,1,0,0,1,0,1,0,0,1,1,1,0,1,1,0,1,0,1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,1,1,1,1,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,1,1,1'

RING2P16_MOD="${RING2P16_MOD:-65536}"  # 2^16
RING2P16_P="${RING2P16_P:-2}"
RING2P16_F="${RING2P16_F:-$RING_F_162_DEFAULT}"
RING2P16_ZETA="${RING2P16_ZETA:-0,1}"  # x

RING2P2_MOD="${RING2P2_MOD:-4}"  # 2^2
RING2P2_P="${RING2P2_P:-2}"
RING2P2_F="${RING2P2_F:-$RING_F_162_DEFAULT}"
RING2P2_ZETA="${RING2P2_ZETA:-0,1}"  # x

ENABLE_FIELD255=0
ENABLE_RING2P16_162=0
ENABLE_F2_256=0
ENABLE_RING2P2_162=0
ENABLE_FIELD64_EXT=0
ENABLE_F2_64_EXT=0
ENABLE_FIELD128_EXT=0
ENABLE_F2_128_EXT=0
ENABLE_F3_40_EXT=0
ENABLE_F3_81_EXT=0
ENABLE_RING2P16_64_EXT=0
ENABLE_RING2P16_128_EXT=0
ENABLE_RING2P2_64_EXT=0
ENABLE_RING2P2_128_EXT=0
SELECTED_CONTEXT_COUNT=0
EFFECTIVE_CPU_SET=""

if (( D_MIN < 0 || D_MAX < D_MIN )); then
  echo "Invalid dimension range: D_MIN=$D_MIN D_MAX=$D_MAX" >&2
  exit 2
fi
if [[ "$ISOLATE_BUILD_DIR" != "0" && "$ISOLATE_BUILD_DIR" != "1" ]]; then
  echo "ISOLATE_BUILD_DIR must be 0 or 1" >&2
  exit 2
fi
if [[ "$CONTINUE_ON_ERROR" != "0" && "$CONTINUE_ON_ERROR" != "1" ]]; then
  echo "CONTINUE_ON_ERROR must be 0 or 1" >&2
  exit 2
fi
if ! [[ "$BENCH_THREADS" =~ ^[0-9]+$ ]]; then
  echo "BENCH_THREADS must be a non-negative integer" >&2
  exit 2
fi
if ! [[ "$K0" =~ ^[1-9][0-9]*$ ]]; then
  echo "K0 must be a positive integer" >&2
  exit 2
fi
if ! is_power_of_two "$K0"; then
  echo "K0 must be a power of two" >&2
  exit 2
fi
if [[ "$CPU_PIN_MODE" != "none" && "$CPU_PIN_MODE" != "manual" && "$CPU_PIN_MODE" != "slot" ]]; then
  echo "CPU_PIN_MODE must be one of: none, manual, slot" >&2
  exit 2
fi
if ! [[ "$RUN_SLOT" =~ ^[0-9]+$ && "$RUN_SLOTS_TOTAL" =~ ^[0-9]+$ ]]; then
  echo "RUN_SLOT and RUN_SLOTS_TOTAL must be non-negative integers" >&2
  exit 2
fi
if (( RUN_SLOTS_TOTAL == 0 )); then
  echo "RUN_SLOTS_TOTAL must be >= 1" >&2
  exit 2
fi
if (( RUN_SLOT >= RUN_SLOTS_TOTAL )); then
  echo "RUN_SLOT must satisfy RUN_SLOT < RUN_SLOTS_TOTAL" >&2
  exit 2
fi
if [[ "$USE_SMT_IN_SLOT" != "0" && "$USE_SMT_IN_SLOT" != "1" ]]; then
  echo "USE_SMT_IN_SLOT must be 0 or 1" >&2
  exit 2
fi
if [[ "$PIN_BUILD" != "0" && "$PIN_BUILD" != "1" ]]; then
  echo "PIN_BUILD must be 0 or 1" >&2
  exit 2
fi
if ! [[ "$F3_40_R" =~ ^[1-9][0-9]*$ && "$F3_81_R" =~ ^[1-9][0-9]*$ ]]; then
  echo "F3_40_R and F3_81_R must be positive integers" >&2
  exit 2
fi
if ! [[ "$F3_40_MOD" =~ ^([2-9]|[1-9][0-9]+)$ && "$F3_81_MOD" =~ ^([2-9]|[1-9][0-9]+)$ ]]; then
  echo "F3_40_MOD and F3_81_MOD must be integers >= 2" >&2
  exit 2
fi

if [[ "$CPU_PIN_MODE" != "none" ]] && ! command -v taskset >/dev/null 2>&1; then
  echo "CPU pinning requested but 'taskset' is not available" >&2
  exit 2
fi

build_slot_cpuset() {
  local slot="$1"
  local total_slots="$2"
  local use_smt="$3"
  local -a cpu_ids
  if (( use_smt == 1 )); then
    if ! command -v nproc >/dev/null 2>&1; then
      echo "slot pinning (USE_SMT_IN_SLOT=1) requires nproc" >&2
      return 1
    fi
    mapfile -t cpu_ids < <(seq 0 $(( $(nproc --all) - 1 )))
  else
    if ! command -v lscpu >/dev/null 2>&1; then
      echo "slot pinning (USE_SMT_IN_SLOT=0) requires lscpu" >&2
      return 1
    fi
    mapfile -t cpu_ids < <(lscpu -p=CPU,CORE,SOCKET | awk -F',' '!/^#/ {key=$2 "-" $3; if (!(key in seen)) {seen[key]=1; print $1}}')
  fi

  local total_cpus="${#cpu_ids[@]}"
  if (( total_cpus == 0 )); then
    echo "Unable to derive CPU topology for slot pinning" >&2
    return 1
  fi

  local base=$(( total_cpus / total_slots ))
  local rem=$(( total_cpus % total_slots ))
  local start=$(( slot * base + (slot < rem ? slot : rem) ))
  local len="$base"
  if (( slot < rem )); then
    len=$((len + 1))
  fi
  if (( len == 0 )); then
    echo "Slot $slot has no CPUs assigned. Reduce RUN_SLOTS_TOTAL." >&2
    return 1
  fi

  local end=$(( start + len ))
  local cpuset=""
  local i
  for (( i = start; i < end; i++ )); do
    if [[ -n "$cpuset" ]]; then
      cpuset+=","
    fi
    cpuset+="${cpu_ids[$i]}"
  done
  echo "$cpuset"
}

case "$CPU_PIN_MODE" in
  none)
    ;;
  manual)
    if [[ -z "$CPU_SET" ]]; then
      echo "CPU_PIN_MODE=manual requires CPU_SET" >&2
      exit 2
    fi
    EFFECTIVE_CPU_SET="$CPU_SET"
    ;;
  slot)
    EFFECTIVE_CPU_SET="$(build_slot_cpuset "$RUN_SLOT" "$RUN_SLOTS_TOTAL" "$USE_SMT_IN_SLOT")"
    ;;
esac

if [[ -n "$EFFECTIVE_CPU_SET" ]]; then
  if ! taskset -c "$EFFECTIVE_CPU_SET" true >/dev/null 2>&1; then
    echo "Invalid CPU set: $EFFECTIVE_CPU_SET" >&2
    exit 2
  fi
fi

if (( BENCH_THREADS > 0 )); then
  export OMP_NUM_THREADS="$BENCH_THREADS"
  export BASEFOLD_MERKLE_MAX_THREADS="${BASEFOLD_MERKLE_MAX_THREADS:-$BENCH_THREADS}"
  export BASEFOLD_VERIFY_QUERY_MAX_THREADS="${BASEFOLD_VERIFY_QUERY_MAX_THREADS:-$BENCH_THREADS}"
fi
if [[ -n "$EFFECTIVE_CPU_SET" ]]; then
  export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
  export OMP_PLACES="${OMP_PLACES:-cores}"
fi

if [[ "$CONTEXTS" == "all" ]]; then
  ENABLE_FIELD255=1
  ENABLE_RING2P16_162=1
  ENABLE_F2_256=1
  ENABLE_RING2P2_162=1
  ENABLE_FIELD64_EXT=1
  ENABLE_F2_64_EXT=1
  ENABLE_FIELD128_EXT=1
  ENABLE_F2_128_EXT=1
  ENABLE_F3_40_EXT=1
  ENABLE_F3_81_EXT=1
  ENABLE_RING2P16_64_EXT=1
  ENABLE_RING2P16_128_EXT=1
  ENABLE_RING2P2_64_EXT=1
  ENABLE_RING2P2_128_EXT=1
else
  IFS=',' read -r -a context_tokens <<< "$CONTEXTS"
  for raw_token in "${context_tokens[@]}"; do
    token="${raw_token//[[:space:]]/}"
    case "$token" in
      field-255)
        ENABLE_FIELD255=1
        ;;
      ring-gr-2p16-162)
        ENABLE_RING2P16_162=1
        ;;
      field-f2p256)
        ENABLE_F2_256=1
        ;;
      ring-gr-2p2-162)
        ENABLE_RING2P2_162=1
        ;;
      field-prime64|field-prime64-ext)
        ENABLE_FIELD64_EXT=1
        ;;
      field-f2p64|field-f2p64-ext)
        ENABLE_F2_64_EXT=1
        ;;
      field-prime128|field-prime128-ext)
        ENABLE_FIELD128_EXT=1
        ;;
      field-f2p128|field-f2p128-ext)
        ENABLE_F2_128_EXT=1
        ;;
      field-f3p40|field-f3p40-ext)
        ENABLE_F3_40_EXT=1
        ;;
      field-f3p81|field-f3p81-ext)
        ENABLE_F3_81_EXT=1
        ;;
      ring-gr-2p16-64-ext)
        ENABLE_RING2P16_64_EXT=1
        ;;
      ring-gr-2p16-128-ext)
        ENABLE_RING2P16_128_EXT=1
        ;;
      ring-gr-2p2-64-ext)
        ENABLE_RING2P2_64_EXT=1
        ;;
      ring-gr-2p2-128-ext)
        ENABLE_RING2P2_128_EXT=1
        ;;
      "")
        ;;
      *)
        echo "Unknown context in CONTEXTS: $token" >&2
        echo "Valid: field-255,ring-gr-2p16-162,field-f2p256,ring-gr-2p2-162,field-prime64-ext,field-f2p64-ext,field-prime128-ext,field-f2p128-ext,field-f3p40-ext,field-f3p81-ext,ring-gr-2p16-64-ext,ring-gr-2p16-128-ext,ring-gr-2p2-64-ext,ring-gr-2p2-128-ext,all" >&2
        exit 2
        ;;
    esac
  done
fi

if (( ENABLE_F3_40_EXT )); then
  deg_f3_40="$(poly_degree_from_coeff_list "$F3_40_F")"
  if [[ "$deg_f3_40" == "ERR" ]]; then
    echo "Invalid F3_40_F: expect comma-separated integer coefficients" >&2
    exit 2
  fi
  if (( deg_f3_40 != F3_40_R )); then
    echo "F3_40 mismatch: deg(F)=$deg_f3_40 but F3_40_R=$F3_40_R" >&2
    echo "Set F3_40_F to an irreducible polynomial of degree F3_40_R." >&2
    exit 2
  fi
fi
if (( ENABLE_F3_81_EXT )); then
  deg_f3_81="$(poly_degree_from_coeff_list "$F3_81_F")"
  if [[ "$deg_f3_81" == "ERR" ]]; then
    echo "Invalid F3_81_F: expect comma-separated integer coefficients" >&2
    exit 2
  fi
  if (( deg_f3_81 != F3_81_R )); then
    echo "F3_81 mismatch: deg(F)=$deg_f3_81 but F3_81_R=$F3_81_R" >&2
    echo "Set F3_81_F to an irreducible polynomial of degree F3_81_R." >&2
    exit 2
  fi
fi

SELECTED_CONTEXT_COUNT=$((ENABLE_FIELD255 + ENABLE_RING2P16_162 + ENABLE_F2_256 + ENABLE_RING2P2_162 + ENABLE_FIELD64_EXT + ENABLE_F2_64_EXT + ENABLE_FIELD128_EXT + ENABLE_F2_128_EXT + ENABLE_F3_40_EXT + ENABLE_F3_81_EXT + ENABLE_RING2P16_64_EXT + ENABLE_RING2P16_128_EXT + ENABLE_RING2P2_64_EXT + ENABLE_RING2P2_128_EXT))
if (( SELECTED_CONTEXT_COUNT == 0 )); then
  echo "No context selected. Set CONTEXTS=all or provide at least one valid context id." >&2
  exit 2
fi

mkdir -p "$OUT_DIR/logs"

RESULT_CSV="$OUT_DIR/results.csv"
RESULT_MD="$OUT_DIR/RESULTS.md"
RUN_AT_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
PROOF_SIZE_SOURCE="bench_pcs_eval fixed-width payload count"

cat > "$RESULT_CSV" <<CSV
context_id,context_label,mode,d,poly_dim,c,k0,lambda,gamma,queries,commit_mean_ms,prove_phase_mean_ms,verifier_mean_ms,proof_size_kb,proof_size_bytes,status,error
CSV

run_and_log() {
  local log_file="$1"
  shift
  {
    printf '+'
    if (( CMD_TIMEOUT_SEC > 0 )); then
      printf ' %q %q' timeout "$CMD_TIMEOUT_SEC"
    fi
    if [[ -n "$EFFECTIVE_CPU_SET" ]]; then
      printf ' %q %q %q' taskset -c "$EFFECTIVE_CPU_SET"
    fi
    printf ' %q' "$@"
    printf '\n'
    if (( CMD_TIMEOUT_SEC > 0 )); then
      if [[ -n "$EFFECTIVE_CPU_SET" ]]; then
        timeout "$CMD_TIMEOUT_SEC" taskset -c "$EFFECTIVE_CPU_SET" "$@"
      else
        timeout "$CMD_TIMEOUT_SEC" "$@"
      fi
    else
      if [[ -n "$EFFECTIVE_CPU_SET" ]]; then
        taskset -c "$EFFECTIVE_CPU_SET" "$@"
      else
        "$@"
      fi
    fi
  } >"$log_file" 2>&1
}

parse_first() {
  local pattern="$1"
  local field_idx="$2"
  local file="$3"
  awk -v pat="$pattern" -v idx="$field_idx" '$0 ~ pat {print $idx; exit}' "$file"
}

first_error_line() {
  local file="$1"
  local line
  line="$(grep -E "Invalid|Error|Unhandled|failed|Killed|timed out|cannot|bad_alloc" "$file" | head -n 1 || true)"
  if [[ -z "$line" ]]; then
    line="$(tail -n 1 "$file" | tr '\r' ' ' || true)"
  fi
  line="${line//$','/;}"
  echo "$line"
}

maybe_abort() {
  local status="$1"
  local msg="$2"
  if [[ "$status" != "ok" && "$CONTINUE_ON_ERROR" == "0" ]]; then
    echo "Abort on error: $msg" >&2
    exit 1
  fi
}

write_row() {
  local context_id="$1"
  local context_label="$2"
  local mode="$3"
  local d="$4"
  local poly_dim="$5"
  local gamma="$6"
  local queries="$7"
  local commit_ms="$8"
  local prove_phase_ms="$9"
  local verifier_ms="${10}"
  local proof_kb="${11}"
  local proof_bytes="${12}"
  local status="${13}"
  local error="${14}"
  local context_label_csv="${context_label//,/;}"
  local error_csv="${error//,/;}"

  echo "${context_id},${context_label_csv},${mode},${d},${poly_dim},${C},${K0},${LAMBDA},${gamma},${queries},${commit_ms},${prove_phase_ms},${verifier_ms},${proof_kb},${proof_bytes},${status},${error_csv}" >> "$RESULT_CSV"
}

run_one_context_d() {
  local context_id="$1"
  local context_label="$2"
  local mode="$3"
  local d="$4"
  shift 4
  local -a bench_args=("$@")

  local poly_dim=$((K0 * (1 << d)))
  local status="ok"
  local error=""
  local gamma="NA"
  local queries="NA"
  local commit_ms="NA"
  local prove_phase_ms="NA"
  local verifier_ms="NA"
  local proof_kb="NA"
  local proof_bytes="NA"
  local calc_p=""
  local calc_r=""
  local calc_m="1"
  local -a eval_extra_args=()

  case "$context_id" in
    field-255)
      calc_p="$FIELD255_MOD"
      calc_r="1"
      ;;
    ring-gr-2p16-162)
      calc_p="2"
      calc_r="162"
      ;;
    field-f2p256)
      calc_p="2"
      calc_r="256"
      ;;
    ring-gr-2p2-162)
      calc_p="2"
      calc_r="162"
      ;;
    field-prime64-ext)
      calc_p="$FIELD64_MOD"
      calc_r="1"
      calc_m="3"
      eval_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_FIELD_EXT_DEG3")
      ;;
    field-f2p64-ext)
      calc_p="2"
      calc_r="64"
      calc_m="3"
      eval_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_F2_64_EXT_DEG3")
      ;;
    field-prime128-ext)
      calc_p="$FIELD128_MOD"
      calc_r="1"
      calc_m="2"
      eval_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_FIELD_EXT_DEG2")
      ;;
    field-f2p128-ext)
      calc_p="2"
      calc_r="128"
      calc_m="2"
      eval_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_F2_128_EXT_DEG2")
      ;;
    field-f3p40-ext)
      calc_p="$F3_40_MOD"
      calc_r="$F3_40_R"
      calc_m="3"
      eval_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_F3_40_EXT_DEG3")
      ;;
    field-f3p81-ext)
      calc_p="$F3_81_MOD"
      calc_r="$F3_81_R"
      calc_m="2"
      eval_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_F3_81_EXT_DEG2")
      ;;
    ring-gr-2p16-64-ext)
      calc_p="2"
      calc_r="64"
      calc_m="3"
      eval_extra_args+=(--use-extension-challenges --ring-challenge-ext "$CHALLENGE_RING_EXT_DEG3")
      ;;
    ring-gr-2p16-128-ext)
      calc_p="2"
      calc_r="128"
      calc_m="2"
      eval_extra_args+=(--use-extension-challenges --ring-challenge-ext "$CHALLENGE_GR_128_EXT_DEG2")
      ;;
    ring-gr-2p2-64-ext)
      calc_p="2"
      calc_r="64"
      calc_m="3"
      eval_extra_args+=(--use-extension-challenges --ring-challenge-ext "$CHALLENGE_RING_EXT_DEG3")
      ;;
    ring-gr-2p2-128-ext)
      calc_p="2"
      calc_r="128"
      calc_m="2"
      eval_extra_args+=(--use-extension-challenges --ring-challenge-ext "$CHALLENGE_GR_128_EXT_DEG2")
      ;;
    *)
      status="context_config_error"
      error="unknown context"
      write_row "$context_id" "$context_label" "$mode" "$d" "$poly_dim" "$gamma" "$queries" "$commit_ms" "$prove_phase_ms" "$verifier_ms" "$proof_kb" "$proof_bytes" "$status" "$error"
      maybe_abort "$status" "unknown context_id=$context_id"
      return
      ;;
  esac

  local prefix="${context_id}_d${d}"
  local calc_log="$OUT_DIR/logs/${prefix}_calc_iopp.log"
  local commit_log="$OUT_DIR/logs/${prefix}_commit.log"
  local eval_log="$OUT_DIR/logs/${prefix}_eval.log"

  echo "  - ${context_label}: d=${d} (poly_dim=${poly_dim} = ${K0}*2^${d})"

  if ! run_and_log "$calc_log" \
      "$BUILD_DIR/calc_iopp_params" \
      --d "$d" --c "$C" --k0 "$K0" --lambda "$LAMBDA" \
      --p "$calc_p" --r "$calc_r" --m "$calc_m" \
      --auto-gamma; then
    status="calc_failed"
    error="$(first_error_line "$calc_log")"
  else
    gamma="$(parse_first "^  gamma    = " 3 "$calc_log")"
    queries="$(parse_first "^  l_min_for_PCS" 3 "$calc_log")"
    if [[ -z "$queries" || "$queries" == "N/A" || ! "$queries" =~ ^[0-9]+$ ]]; then
      status="no_feasible_queries"
      error="l_min_for_PCS unavailable"
      queries="NA"
    fi
  fi

  if ! run_and_log "$commit_log" \
      "$BUILD_DIR/bench_pcs_commit" \
      --mode "$mode" \
      "${bench_args[@]}" \
      --c "$C" --k0 "$K0" --d "$d" \
      --warmup "$COMMIT_WARMUP" --reps "$COMMIT_REPS" \
      --seed "$SEED"; then
    if [[ "$status" == "ok" ]]; then
      status="commit_failed"
      error="$(first_error_line "$commit_log")"
    fi
  else
    commit_ms="$(parse_first "^  commit[[:space:]]+mean" 3 "$commit_log")"
    if [[ -z "$commit_ms" ]]; then
      commit_ms="NA"
    fi
  fi

  if [[ "$status" == "ok" ]]; then
    if ! run_and_log "$eval_log" \
        "$BUILD_DIR/bench_pcs_eval" \
        --mode "$mode" \
        "${eval_extra_args[@]}" \
        "${bench_args[@]}" \
        --c "$C" --k0 "$K0" --d "$d" \
        --queries "$queries" \
        --warmup "$EVAL_WARMUP" --reps "$EVAL_REPS" \
        --seed "$SEED"; then
      status="eval_failed"
      error="$(first_error_line "$eval_log")"
    else
      prove_phase_ms="$(parse_first "prove-phase mean" 3 "$eval_log")"
      verifier_ms="$(parse_first "verifier mean" 3 "$eval_log")"
      if [[ -z "$prove_phase_ms" ]]; then
        prove_phase_ms="NA"
      fi
      if [[ -z "$verifier_ms" ]]; then
        verifier_ms="NA"
      fi
      proof_kb="$(parse_first "proof size" 3 "$eval_log")"
      proof_bytes="$(awk '/proof size/{gsub(/[^0-9]/, "", $5); print $5; exit}' "$eval_log")"
      if [[ -z "$proof_kb" || -z "$proof_bytes" ]]; then
        status="proof_parse_failed"
        error="missing proof size line in bench_pcs_eval output"
      fi
    fi
  fi

  if [[ -z "$error" ]]; then
    error="-"
  fi
  write_row "$context_id" "$context_label" "$mode" "$d" "$poly_dim" "$gamma" "$queries" "$commit_ms" "$prove_phase_ms" "$verifier_ms" "$proof_kb" "$proof_bytes" "$status" "$error"
  maybe_abort "$status" "${context_label} d=${d} status=${status}"
}

echo "[1/4] Configure Release build in: $BUILD_DIR"
if [[ -n "$EFFECTIVE_CPU_SET" && "$PIN_BUILD" == "1" ]]; then
  taskset -c "$EFFECTIVE_CPU_SET" cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
else
  cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
fi

echo "[2/4] Build required benchmarks/tools"
if [[ -n "$EFFECTIVE_CPU_SET" && "$PIN_BUILD" == "1" ]]; then
  taskset -c "$EFFECTIVE_CPU_SET" cmake --build "$BUILD_DIR" \
    --target bench_pcs_commit bench_pcs_eval calc_iopp_params \
    --parallel
else
  cmake --build "$BUILD_DIR" \
    --target bench_pcs_commit bench_pcs_eval calc_iopp_params \
    --parallel
fi

echo "[3/4] Run sweep: d in [$D_MIN, $D_MAX], selected contexts = $SELECTED_CONTEXT_COUNT (contexts serial, each bench may use threads)"
for d in $(seq "$D_MIN" "$D_MAX"); do
  echo "[d=$d]"
  if (( ENABLE_FIELD255 )); then
    run_one_context_d \
      "field-255" "Field-255" "field" "$d" \
      --field-mod "$FIELD255_MOD" \
      --field-F "$FIELD255_F" \
      --field-zeta "$FIELD255_ZETA"
  fi
  if (( ENABLE_RING2P16_162 )); then
    run_one_context_d \
      "ring-gr-2p16-162" "GR(2^16,162)" "ring" "$d" \
      --ring-mod "$RING2P16_MOD" \
      --ring-p "$RING2P16_P" \
      --ring-F "$RING2P16_F" \
      --ring-zeta "$RING2P16_ZETA"
  fi
  if (( ENABLE_F2_256 )); then
    run_one_context_d \
      "field-f2p256" "F_2^256" "field" "$d" \
      --field-mod "$F2_256_MOD" \
      --field-F "$F2_256_F" \
      --field-zeta "$F2_256_ZETA"
  fi
  if (( ENABLE_RING2P2_162 )); then
    run_one_context_d \
      "ring-gr-2p2-162" "GR(2^2,162)" "ring" "$d" \
      --ring-mod "$RING2P2_MOD" \
      --ring-p "$RING2P2_P" \
      --ring-F "$RING2P2_F" \
      --ring-zeta "$RING2P2_ZETA"
  fi
  if (( ENABLE_FIELD64_EXT )); then
    run_one_context_d \
      "field-prime64-ext" "Field-Prime64 (ext-challenge)" "field" "$d" \
      --field-mod "$FIELD64_MOD" \
      --field-F "$FIELD64_F" \
      --field-zeta "$FIELD64_ZETA"
  fi
  if (( ENABLE_F2_64_EXT )); then
    run_one_context_d \
      "field-f2p64-ext" "F_2^64 (ext-challenge)" "field" "$d" \
      --field-mod "$F2_64_MOD" \
      --field-F "$F2_64_F" \
      --field-zeta "$F2_64_ZETA"
  fi
  if (( ENABLE_FIELD128_EXT )); then
    run_one_context_d \
      "field-prime128-ext" "Field-Prime128 (ext-challenge)" "field" "$d" \
      --field-mod "$FIELD128_MOD" \
      --field-F "$FIELD128_F" \
      --field-zeta "$FIELD128_ZETA"
  fi
  if (( ENABLE_F2_128_EXT )); then
    run_one_context_d \
      "field-f2p128-ext" "F_2^128 (ext-challenge)" "field" "$d" \
      --field-mod "$F2_128_MOD" \
      --field-F "$F2_128_F" \
      --field-zeta "$F2_128_ZETA"
  fi
  if (( ENABLE_F3_40_EXT )); then
    run_one_context_d \
      "field-f3p40-ext" "F_3^40 (ext-challenge)" "field" "$d" \
      --field-mod "$F3_40_MOD" \
      --field-F "$F3_40_F" \
      --field-zeta "$F3_40_ZETA"
  fi
  if (( ENABLE_F3_81_EXT )); then
    run_one_context_d \
      "field-f3p81-ext" "F_3^81 (ext-challenge)" "field" "$d" \
      --field-mod "$F3_81_MOD" \
      --field-F "$F3_81_F" \
      --field-zeta "$F3_81_ZETA"
  fi
  if (( ENABLE_RING2P16_64_EXT )); then
    run_one_context_d \
      "ring-gr-2p16-64-ext" "GR(2^16,64) (ext-challenge)" "ring" "$d" \
      --ring-mod "$RING2P16_MOD" \
      --ring-p "$RING2P16_P" \
      --ring-F "$RING_F_64" \
      --ring-zeta "$RING2P16_ZETA"
  fi
  if (( ENABLE_RING2P16_128_EXT )); then
    run_one_context_d \
      "ring-gr-2p16-128-ext" "GR(2^16,128) (ext-challenge)" "ring" "$d" \
      --ring-mod "$RING2P16_MOD" \
      --ring-p "$RING2P16_P" \
      --ring-F "$RING_F_128" \
      --ring-zeta "$RING2P16_ZETA"
  fi
  if (( ENABLE_RING2P2_64_EXT )); then
    run_one_context_d \
      "ring-gr-2p2-64-ext" "GR(2^2,64) (ext-challenge)" "ring" "$d" \
      --ring-mod "$RING2P2_MOD" \
      --ring-p "$RING2P2_P" \
      --ring-F "$RING_F_64" \
      --ring-zeta "$RING2P2_ZETA"
  fi
  if (( ENABLE_RING2P2_128_EXT )); then
    run_one_context_d \
      "ring-gr-2p2-128-ext" "GR(2^2,128) (ext-challenge)" "ring" "$d" \
      --ring-mod "$RING2P2_MOD" \
      --ring-p "$RING2P2_P" \
      --ring-F "$RING_F_128" \
      --ring-zeta "$RING2P2_ZETA"
  fi
done

echo "[4/4] Build markdown summary"
{
  echo "# Release Sweep Results: c=$C, lambda=$LAMBDA"
  echo ""
  echo "- run_at_utc: $RUN_AT_UTC"
  echo "- run_id: $RUN_ID"
  echo "- build_dir: $BUILD_DIR"
  echo "- isolate_build_dir: $ISOLATE_BUILD_DIR"
  echo "- output_dir: $OUT_DIR"
  echo "- d_range: [$D_MIN, $D_MAX] (poly_dim = k_d = ${K0}*2^d)"
  echo "- contexts: $CONTEXTS"
  echo "- bench_threads: $BENCH_THREADS (set 0 to use runtime default)"
  echo "- proof_size_source: $PROOF_SIZE_SOURCE"
  echo "- cpu_pin_mode: $CPU_PIN_MODE"
  if [[ -n "$EFFECTIVE_CPU_SET" ]]; then
    echo "- cpu_set: $EFFECTIVE_CPU_SET"
    if [[ "$CPU_PIN_MODE" == "slot" ]]; then
      echo "- run_slot: $RUN_SLOT / $RUN_SLOTS_TOTAL (use_smt_in_slot=$USE_SMT_IN_SLOT)"
    fi
  fi
  echo "- pin_build: $PIN_BUILD"
  echo "- continue_on_error: $CONTINUE_ON_ERROR"
  echo "- cmd_timeout_sec: $CMD_TIMEOUT_SEC"
  echo ""
  echo "## Results"
  echo ""
  echo "| context | d | poly_dim | gamma | queries | commit mean ms | prove-phase mean ms | verifier mean ms | proof size KB | status |"
  echo "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"
  tail -n +2 "$RESULT_CSV" | awk -F',' '{printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", $2, $4, $5, $9, $10, $11, $12, $13, $14, $16}'
  echo ""
  echo "## Status Counts"
  echo ""
  tail -n +2 "$RESULT_CSV" | awk -F',' '{cnt[$16]++} END {for (k in cnt) printf "- %s: %d\n", k, cnt[k]}'
  echo ""
  echo "Raw csv: $RESULT_CSV"
  echo "Raw logs dir: $OUT_DIR/logs"
} > "$RESULT_MD"

echo "Done."
echo "Summary markdown: $RESULT_MD"
echo "Summary csv:      $RESULT_CSV"
