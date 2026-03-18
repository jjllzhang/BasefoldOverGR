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

log2_power_of_two() {
  local value="$1"
  local log2=0
  while (( value > 1 )); do
    value=$((value >> 1))
    ((log2 += 1))
  done
  echo "$log2"
}

# Target profile: rate = 1/4 (c=4), security = 128 bits.
C="${C:-4}"
K0="${K0:-1}"
LAMBDA="${LAMBDA:-128}"
D_MIN="${D_MIN:-3}"
D_MAX="${D_MAX:-29}"
COMPILER_KAPPA="${COMPILER_KAPPA:-}"
COMPILER_ELL_MIN="${COMPILER_ELL_MIN:-}"
COMPILER_ELL_MAX="${COMPILER_ELL_MAX:-}"
COMPILER_FAMILY="${COMPILER_FAMILY:-}"

COMMIT_WARMUP="${COMMIT_WARMUP:-1}"
COMMIT_REPS="${COMMIT_REPS:-3}"
PROVE_WARMUP="${PROVE_WARMUP:-1}"
PROVE_REPS="${PROVE_REPS:-3}"
VERIFY_WARMUP="${VERIFY_WARMUP:-1}"
VERIFY_REPS="${VERIFY_REPS:-3}"
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

mkdir -p "$OUT_DIR/logs"
RUN_SUITE="${RUN_SUITE:-basefold_release}"

RUN_BASEFOLD_RELEASE=0
RUN_COMPILER_EVAL=0
COMPILER_EVAL_FAMILY=""
case "$RUN_SUITE" in
  basefold_release)
    RUN_BASEFOLD_RELEASE=1
    ;;
  compiler_eval_ring_switch)
    RUN_COMPILER_EVAL=1
    COMPILER_EVAL_FAMILY="ring_switch"
    ;;
  compiler_eval_frobenius)
    RUN_COMPILER_EVAL=1
    COMPILER_EVAL_FAMILY="frobenius"
    ;;
  compiler_eval)
    RUN_COMPILER_EVAL=1
    ;;
  *)
    echo "RUN_SUITE must be one of: basefold_release, compiler_eval_ring_switch, compiler_eval_frobenius, compiler_eval" >&2
    exit 2
    ;;
esac

if (( RUN_COMPILER_EVAL )) && [[ -z "$COMPILER_EVAL_FAMILY" ]]; then
  if [[ "$COMPILER_FAMILY" == "ring_switch" || "$COMPILER_FAMILY" == "frobenius" ]]; then
    COMPILER_EVAL_FAMILY="$COMPILER_FAMILY"
  else
    echo "Set RUN_SUITE=compiler_eval_ring_switch or compiler_eval_frobenius (or keep RUN_SUITE=compiler_eval and set COMPILER_FAMILY=ring_switch|frobenius)." >&2
    exit 2
  fi
fi

if (( RUN_COMPILER_EVAL )) && [[ -n "$COMPILER_FAMILY" ]] && [[ "$COMPILER_FAMILY" != "$COMPILER_EVAL_FAMILY" ]]; then
  echo "COMPILER_FAMILY=$COMPILER_FAMILY conflicts with RUN_SUITE-selected family $COMPILER_EVAL_FAMILY" >&2
  exit 2
fi

if (( RUN_BASEFOLD_RELEASE )); then
  if ! [[ "$D_MIN" =~ ^[0-9]+$ && "$D_MAX" =~ ^[0-9]+$ ]]; then
    echo "D_MIN and D_MAX must be non-negative integers" >&2
    exit 2
  fi
  if (( D_MAX < D_MIN )); then
    echo "Invalid dimension range: D_MIN=$D_MIN D_MAX=$D_MAX" >&2
    exit 2
  fi
fi

if (( RUN_COMPILER_EVAL )); then
  if ! [[ "$COMPILER_KAPPA" =~ ^[1-9][0-9]*$ ]]; then
    echo "COMPILER_KAPPA must be a positive integer when compiler_eval suite is enabled" >&2
    exit 2
  fi
  if ! [[ "$COMPILER_ELL_MIN" =~ ^[0-9]+$ && "$COMPILER_ELL_MAX" =~ ^[0-9]+$ ]]; then
    echo "COMPILER_ELL_MIN and COMPILER_ELL_MAX must be non-negative integers when compiler_eval suite is enabled" >&2
    exit 2
  fi
  if (( COMPILER_ELL_MIN < COMPILER_KAPPA )); then
    echo "COMPILER_ELL_MIN must satisfy COMPILER_ELL_MIN >= COMPILER_KAPPA" >&2
    exit 2
  fi
  if (( COMPILER_ELL_MAX < COMPILER_ELL_MIN )); then
    echo "COMPILER_ELL_MAX must satisfy COMPILER_ELL_MAX >= COMPILER_ELL_MIN" >&2
    exit 2
  fi
fi

if (( SELECTED_CONTEXT_COUNT == 0 )); then
  echo "No context selected. Set CONTEXTS=all or provide at least one valid context id." >&2
  exit 2
fi

COMPILER_EVAL_RESULT_CSV="$OUT_DIR/compiler_eval_results.csv"
BACKEND_EVAL_RESULT_CSV="$OUT_DIR/backend_eval_results.csv"
RESULT_MD="$OUT_DIR/RESULTS.md"
RUN_AT_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
BASEFOLD_RELEASE_SOURCE="family=basefold rows use bench_basefold_pcs_commit + bench_basefold_pcs_eval"
RELEASE_QUERY_SOURCE="basefold_release/compiler_eval* use calc_iopp_params and parse default l_min_for_IOPP"
if (( RUN_COMPILER_EVAL )); then
  COMPILER_EVAL_SOURCE="family=${COMPILER_EVAL_FAMILY} rows use bench_z2k_${COMPILER_EVAL_FAMILY}_eval with commit split"
else
  COMPILER_EVAL_SOURCE="n/a"
fi
BACKEND_EVAL_SOURCE="family=basefold rows use bench_basefold_pcs_commit + bench_basefold_pcs_eval only"
COMPILER_EVAL_LAYOUT_SOURCE="compiler_eval_* rows use fixed COMPILER_KAPPA with ell range, derive d=ell-kappa, and record poly_dim=2^ell"
COMPILER_BACKEND_SOURCE="current z2k compiler benches build BaseFold backend params with k0=1 only"

cat > "$COMPILER_EVAL_RESULT_CSV" <<CSV
family,context_id,context_label,mode,d,poly_dim,c,k0,lambda,gamma,queries,compiler_ell,ell_prime,compiler_kappa,outer_commit_mean_ms,backend_commit_mean_ms,commit_total_mean_ms,prove_total_mean_ms,prove_outer_mean_ms,prove_backend_mean_ms,verify_total_mean_ms,verify_outer_mean_ms,verify_backend_mean_ms,outer_proof_size_kb,outer_proof_size_bytes,proof_size_kb,proof_size_bytes,status,error
CSV

cat > "$BACKEND_EVAL_RESULT_CSV" <<CSV
family,context_id,context_label,mode,d,poly_dim,c,k0,lambda,gamma,queries,commit_mean_ms,prove_phase_mean_ms,verifier_mean_ms,proof_size_kb,proof_size_bytes,status,error
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

parse_size_kb() {
  local pattern="$1"
  local file="$2"
  awk -v pat="$pattern" '
    $0 ~ pat {
      for (i = 1; i <= NF; ++i) {
        if ($i == "KB") {
          print $(i - 1)
          exit
        }
      }
    }' "$file"
}

parse_size_bytes() {
  local pattern="$1"
  local file="$2"
  awk -v pat="$pattern" '
    $0 ~ pat {
      if (match($0, /\(([0-9]+) B\)/)) {
        print substr($0, RSTART + 1, RLENGTH - 4)
        exit
      }
    }' "$file"
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

compiler_precheck() {
  local family="$1"
  local context_extension_degree="$2"
  local compiler_kappa="$3"
  local context_scalar_modulus="$4"
  if [[ "$context_extension_degree" == "ERR" || ! "$context_extension_degree" =~ ^[0-9]+$ ]]; then
    echo "context_config_error,invalid ring polynomial coefficients"
  elif [[ "$compiler_kappa" == "NA" ]] || (( compiler_kappa < 1 )); then
    echo "unsupported_kappa,compiler family requires kappa >= 1"
  elif ! is_power_of_two "$context_extension_degree"; then
    echo "unsupported_context_degree,compiler family requires ring extension degree deg(F)=2^kappa"
  elif (( context_extension_degree != (1 << compiler_kappa) )); then
    echo "unsupported_context_degree,selected ring context degree deg(F)=${context_extension_degree} does not match compiler_kappa=${compiler_kappa}"
  elif [[ "$family" == "frobenius" && "$context_scalar_modulus" =~ ^[0-9]+$ ]] && (( context_scalar_modulus == 4 )); then
    echo "disabled_gr2p2_context,frobenius compiler bench disabled for GR(2^2;r) contexts"
  elif (( K0 != 1 )); then
    echo "unsupported_k0,current compiler benches build BaseFold backend params with k0=1 only"
  else
    echo "ok,-"
  fi
}

write_compiler_eval_row() {
  local family="$1"
  local context_id="$2"
  local context_label="$3"
  local mode="$4"
  local d="$5"
  local poly_dim="$6"
  local gamma="$7"
  local queries="$8"
  local compiler_ell="$9"
  local ell_prime="${10}"
  local compiler_kappa="${11}"
  local outer_commit_ms="${12}"
  local backend_commit_ms="${13}"
  local commit_total_ms="${14}"
  local prove_total_ms="${15}"
  local prove_outer_ms="${16}"
  local prove_backend_ms="${17}"
  local verify_total_ms="${18}"
  local verify_outer_ms="${19}"
  local verify_backend_ms="${20}"
  local outer_proof_kb="${21}"
  local outer_proof_bytes="${22}"
  local proof_kb="${23}"
  local proof_bytes="${24}"
  local status="${25}"
  local error="${26}"
  local context_label_csv="${context_label//,/;}"
  local error_csv="${error//,/;}"

  echo "${family},${context_id},${context_label_csv},${mode},${d},${poly_dim},${C},${K0},${LAMBDA},${gamma},${queries},${compiler_ell},${ell_prime},${compiler_kappa},${outer_commit_ms},${backend_commit_ms},${commit_total_ms},${prove_total_ms},${prove_outer_ms},${prove_backend_ms},${verify_total_ms},${verify_outer_ms},${verify_backend_ms},${outer_proof_kb},${outer_proof_bytes},${proof_kb},${proof_bytes},${status},${error_csv}" >> "$COMPILER_EVAL_RESULT_CSV"
}

write_backend_eval_row() {
  local family="$1"
  local context_id="$2"
  local context_label="$3"
  local mode="$4"
  local d="$5"
  local poly_dim="$6"
  local gamma="$7"
  local queries="$8"
  local commit_ms="$9"
  local prove_ms="${10}"
  local verify_ms="${11}"
  local proof_kb="${12}"
  local proof_bytes="${13}"
  local status="${14}"
  local error="${15}"
  local context_label_csv="${context_label//,/;}"
  local error_csv="${error//,/;}"

  echo "${family},${context_id},${context_label_csv},${mode},${d},${poly_dim},${C},${K0},${LAMBDA},${gamma},${queries},${commit_ms},${prove_ms},${verify_ms},${proof_kb},${proof_bytes},${status},${error_csv}" >> "$BACKEND_EVAL_RESULT_CSV"
}

run_compiler_eval_row() {
  local family="$1"
  local context_id="$2"
  local context_label="$3"
  local mode="$4"
  local d="$5"
  local poly_dim="$6"
  local gamma="$7"
  local queries="$8"
  local query_status="$9"
  local query_error="${10}"
  local compiler_ell="${11}"
  local ell_prime="${12}"
  local compiler_kappa="${13}"
  local context_extension_degree="${14}"
  local context_scalar_modulus="${15}"
  shift 15
  local -a bench_args=("$@")

  local eval_bin=""
  case "$family" in
    ring_switch)
      eval_bin="bench_z2k_ring_switch_eval"
      ;;
    frobenius)
      eval_bin="bench_z2k_frobenius_eval"
      ;;
    *)
      echo "Unknown compiler family: $family" >&2
      exit 1
      ;;
  esac

  local status=""
  local error=""
  IFS=',' read -r status error < <(compiler_precheck "$family" "$context_extension_degree" "$compiler_kappa" "$context_scalar_modulus")
  if [[ "$status" != "ok" ]]; then
    write_compiler_eval_row "$family" "$context_id" "$context_label" "$mode" "$d" "$poly_dim" "$gamma" "$queries" "$compiler_ell" "$ell_prime" "$compiler_kappa" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "$status" "$error"
    maybe_abort "$status" "${context_label} d=${d} ${family}:eval status=${status}"
    return
  fi
  if [[ "$query_status" != "ok" ]]; then
    write_compiler_eval_row "$family" "$context_id" "$context_label" "$mode" "$d" "$poly_dim" "$gamma" "$queries" "$compiler_ell" "$ell_prime" "$compiler_kappa" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA" "$query_status" "$query_error"
    maybe_abort "$query_status" "${context_label} d=${d} ${family}:eval status=${query_status}"
    return
  fi

  local log_file="$OUT_DIR/logs/${context_id}_d${d}_${family}_eval.log"
  local outer_commit_ms="NA"
  local backend_commit_ms="NA"
  local commit_total_ms="NA"
  local prove_total_ms="NA"
  local prove_outer_ms="NA"
  local prove_backend_ms="NA"
  local verify_total_ms="NA"
  local verify_outer_ms="NA"
  local verify_backend_ms="NA"
  local outer_proof_kb="NA"
  local outer_proof_bytes="NA"
  local proof_kb="NA"
  local proof_bytes="NA"
  status="ok"
  error="-"

  if ! run_and_log "$log_file" \
      "$BUILD_DIR/$eval_bin" \
      "${bench_args[@]}" \
      --c "$C" --ell "$compiler_ell" --kappa "$compiler_kappa" \
      --queries "$queries" \
      --warmup "$EVAL_WARMUP" --reps "$EVAL_REPS" \
      --seed "$SEED"; then
    status="${family}_eval_failed"
    error="$(first_error_line "$log_file")"
  else
    outer_commit_ms="$(parse_first "outer commit mean" 4 "$log_file")"
    backend_commit_ms="$(parse_first "backend commit mean" 4 "$log_file")"
    commit_total_ms="$(parse_first "commit total mean" 4 "$log_file")"
    prove_total_ms="$(parse_first "prove-phase mean" 3 "$log_file")"
    prove_outer_ms="$(parse_first "outer prover mean" 4 "$log_file")"
    prove_backend_ms="$(parse_first "backend prover mean" 4 "$log_file")"
    verify_total_ms="$(parse_first "verifier mean" 3 "$log_file")"
    verify_outer_ms="$(parse_first "outer verifier mean" 4 "$log_file")"
    verify_backend_ms="$(parse_first "backend verifier mean" 4 "$log_file")"
    outer_proof_kb="$(parse_size_kb "outer proof size" "$log_file")"
    outer_proof_bytes="$(parse_size_bytes "outer proof size" "$log_file")"
    proof_kb="$(parse_size_kb "^  proof size" "$log_file")"
    proof_bytes="$(parse_size_bytes "^  proof size" "$log_file")"
    [[ -n "$outer_commit_ms" ]] || outer_commit_ms="NA"
    [[ -n "$backend_commit_ms" ]] || backend_commit_ms="NA"
    [[ -n "$commit_total_ms" ]] || commit_total_ms="NA"
    [[ -n "$prove_total_ms" ]] || prove_total_ms="NA"
    [[ -n "$prove_outer_ms" ]] || prove_outer_ms="NA"
    [[ -n "$prove_backend_ms" ]] || prove_backend_ms="NA"
    [[ -n "$verify_total_ms" ]] || verify_total_ms="NA"
    [[ -n "$verify_outer_ms" ]] || verify_outer_ms="NA"
    [[ -n "$verify_backend_ms" ]] || verify_backend_ms="NA"
    [[ -n "$outer_proof_kb" ]] || outer_proof_kb="NA"
    [[ -n "$outer_proof_bytes" ]] || outer_proof_bytes="NA"
    [[ -n "$proof_kb" ]] || proof_kb="NA"
    [[ -n "$proof_bytes" ]] || proof_bytes="NA"
  fi

  write_compiler_eval_row "$family" "$context_id" "$context_label" "$mode" "$d" "$poly_dim" "$gamma" "$queries" "$compiler_ell" "$ell_prime" "$compiler_kappa" "$outer_commit_ms" "$backend_commit_ms" "$commit_total_ms" "$prove_total_ms" "$prove_outer_ms" "$prove_backend_ms" "$verify_total_ms" "$verify_outer_ms" "$verify_backend_ms" "$outer_proof_kb" "$outer_proof_bytes" "$proof_kb" "$proof_bytes" "$status" "$error"
  maybe_abort "$status" "${context_label} d=${d} ${family}:eval status=${status}"
}

run_one_context_d() {
  local context_id="$1"
  local context_label="$2"
  local mode="$3"
  local d="$4"
  shift 4
  local -a bench_args=("$@")

  local poly_dim=$((K0 * (1 << d)))
  local calc_p=""
  local calc_r=""
  local calc_m="1"
  local context_extension_degree="NA"
  local -a phase_extra_args=()

  case "$context_id" in
    field-255)
      calc_p="$FIELD255_MOD"
      calc_r="1"
      ;;
    ring-gr-2p16-162)
      calc_p="2"
      calc_r="162"
      context_extension_degree="$(poly_degree_from_coeff_list "$RING2P16_F")"
      ;;
    field-f2p256)
      calc_p="2"
      calc_r="256"
      ;;
    ring-gr-2p2-162)
      calc_p="2"
      calc_r="162"
      context_extension_degree="$(poly_degree_from_coeff_list "$RING2P2_F")"
      ;;
    field-prime64-ext)
      calc_p="$FIELD64_MOD"
      calc_r="1"
      calc_m="3"
      phase_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_FIELD_EXT_DEG3")
      ;;
    field-f2p64-ext)
      calc_p="2"
      calc_r="64"
      calc_m="3"
      phase_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_F2_64_EXT_DEG3")
      ;;
    field-prime128-ext)
      calc_p="$FIELD128_MOD"
      calc_r="1"
      calc_m="2"
      phase_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_FIELD_EXT_DEG2")
      ;;
    field-f2p128-ext)
      calc_p="2"
      calc_r="128"
      calc_m="2"
      phase_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_F2_128_EXT_DEG2")
      ;;
    field-f3p40-ext)
      calc_p="$F3_40_MOD"
      calc_r="$F3_40_R"
      calc_m="3"
      phase_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_F3_40_EXT_DEG3")
      ;;
    field-f3p81-ext)
      calc_p="$F3_81_MOD"
      calc_r="$F3_81_R"
      calc_m="2"
      phase_extra_args+=(--use-extension-challenges --field-challenge-ext "$CHALLENGE_F3_81_EXT_DEG2")
      ;;
    ring-gr-2p16-64-ext)
      calc_p="2"
      calc_r="64"
      calc_m="3"
      context_extension_degree="$(poly_degree_from_coeff_list "$RING_F_64")"
      phase_extra_args+=(--use-extension-challenges --ring-challenge-ext "$CHALLENGE_RING_EXT_DEG3")
      ;;
    ring-gr-2p16-128-ext)
      calc_p="2"
      calc_r="128"
      calc_m="2"
      context_extension_degree="$(poly_degree_from_coeff_list "$RING_F_128")"
      phase_extra_args+=(--use-extension-challenges --ring-challenge-ext "$CHALLENGE_GR_128_EXT_DEG2")
      ;;
    ring-gr-2p2-64-ext)
      calc_p="2"
      calc_r="64"
      calc_m="3"
      context_extension_degree="$(poly_degree_from_coeff_list "$RING_F_64")"
      phase_extra_args+=(--use-extension-challenges --ring-challenge-ext "$CHALLENGE_RING_EXT_DEG3")
      ;;
    ring-gr-2p2-128-ext)
      calc_p="2"
      calc_r="128"
      calc_m="2"
      context_extension_degree="$(poly_degree_from_coeff_list "$RING_F_128")"
      phase_extra_args+=(--use-extension-challenges --ring-challenge-ext "$CHALLENGE_GR_128_EXT_DEG2")
      ;;
    *)
      echo "unknown context_id=$context_id" >&2
      exit 1
      ;;
  esac

  local gamma="NA"
  local release_queries="NA"
  local release_query_status="ok"
  local release_query_error="-"
  local prefix="${context_id}_d${d}"

  echo "  - ${context_label}: d=${d} (poly_dim=${poly_dim} = ${K0}*2^${d})"
  local calc_log="$OUT_DIR/logs/${prefix}_calc_iopp.log"
  if ! run_and_log "$calc_log" \
      "$BUILD_DIR/calc_iopp_params" \
      --d "$d" --c "$C" --k0 "$K0" --lambda "$LAMBDA" \
      --p "$calc_p" --r "$calc_r" --m "$calc_m" \
      --auto-gamma; then
    release_query_status="calc_failed"
    release_query_error="$(first_error_line "$calc_log")"
  else
    gamma="$(parse_first "^  gamma    = " 3 "$calc_log")"
    [[ -n "$gamma" ]] || gamma="NA"
    release_queries="$(parse_first "^  l_min_for_IOPP" 3 "$calc_log")"
    if [[ -z "$release_queries" || "$release_queries" == "N/A" || ! "$release_queries" =~ ^[0-9]+$ ]]; then
      release_query_status="no_feasible_queries"
      release_query_error="l_min_for_IOPP unavailable"
      release_queries="NA"
    fi
  fi

  local basefold_commit_ms="NA"
  local basefold_prove_ms="NA"
  local basefold_verify_ms="NA"
  local basefold_proof_kb="NA"
  local basefold_proof_bytes="NA"
  local basefold_status="ok"
  local basefold_error="-"
  local commit_log="$OUT_DIR/logs/${prefix}_basefold_commit.log"
  local eval_log="$OUT_DIR/logs/${prefix}_basefold_eval.log"

  if [[ "$release_query_status" != "ok" ]]; then
    basefold_status="$release_query_status"
    basefold_error="$release_query_error"
  fi

  if ! run_and_log "$commit_log" \
      "$BUILD_DIR/bench_basefold_pcs_commit" \
      --mode "$mode" \
      "${bench_args[@]}" \
      --c "$C" --k0 "$K0" --d "$d" \
      --warmup "$COMMIT_WARMUP" --reps "$COMMIT_REPS" \
      --seed "$SEED"; then
    if [[ "$basefold_status" == "ok" ]]; then
      basefold_status="basefold_commit_failed"
      basefold_error="$(first_error_line "$commit_log")"
    fi
  else
    basefold_commit_ms="$(parse_first "^  commit[[:space:]]+mean" 3 "$commit_log")"
    [[ -n "$basefold_commit_ms" ]] || basefold_commit_ms="NA"
  fi

  if [[ "$release_query_status" == "ok" ]]; then
    if ! run_and_log "$eval_log" \
        "$BUILD_DIR/bench_basefold_pcs_eval" \
        --mode "$mode" \
        "${phase_extra_args[@]}" \
        "${bench_args[@]}" \
        --c "$C" --k0 "$K0" --d "$d" \
        --queries "$release_queries" \
        --warmup "$EVAL_WARMUP" --reps "$EVAL_REPS" \
        --seed "$SEED"; then
      if [[ "$basefold_status" == "ok" ]]; then
        basefold_status="basefold_eval_failed"
        basefold_error="$(first_error_line "$eval_log")"
      fi
    else
      basefold_prove_ms="$(parse_first "prove-phase mean" 3 "$eval_log")"
      basefold_verify_ms="$(parse_first "verifier mean" 3 "$eval_log")"
      basefold_proof_kb="$(parse_size_kb "^  proof size" "$eval_log")"
      basefold_proof_bytes="$(parse_size_bytes "^  proof size" "$eval_log")"
      [[ -n "$basefold_prove_ms" ]] || basefold_prove_ms="NA"
      [[ -n "$basefold_verify_ms" ]] || basefold_verify_ms="NA"
      [[ -n "$basefold_proof_kb" ]] || basefold_proof_kb="NA"
      [[ -n "$basefold_proof_bytes" ]] || basefold_proof_bytes="NA"
    fi
  fi

  write_backend_eval_row "basefold" "$context_id" "$context_label" "$mode" "$d" "$poly_dim" "$gamma" "$release_queries" "$basefold_commit_ms" "$basefold_prove_ms" "$basefold_verify_ms" "$basefold_proof_kb" "$basefold_proof_bytes" "$basefold_status" "$basefold_error"
  maybe_abort "$basefold_status" "${context_label} d=${d} basefold status=${basefold_status}"
}

resolve_compiler_context_metadata() {
  local context_id="$1"
  COMPILER_CONTEXT_CALC_P=""
  COMPILER_CONTEXT_CALC_R=""
  COMPILER_CONTEXT_CALC_M="1"
  COMPILER_CONTEXT_EXTENSION_DEGREE="NA"
  COMPILER_CONTEXT_SCALAR_MODULUS="NA"

  case "$context_id" in
    ring-gr-2p16-162)
      COMPILER_CONTEXT_CALC_P="2"
      COMPILER_CONTEXT_CALC_R="162"
      COMPILER_CONTEXT_EXTENSION_DEGREE="$(poly_degree_from_coeff_list "$RING2P16_F")"
      COMPILER_CONTEXT_SCALAR_MODULUS="$RING2P16_MOD"
      ;;
    ring-gr-2p2-162)
      COMPILER_CONTEXT_CALC_P="2"
      COMPILER_CONTEXT_CALC_R="162"
      COMPILER_CONTEXT_EXTENSION_DEGREE="$(poly_degree_from_coeff_list "$RING2P2_F")"
      COMPILER_CONTEXT_SCALAR_MODULUS="$RING2P2_MOD"
      ;;
    ring-gr-2p16-64-ext)
      COMPILER_CONTEXT_CALC_P="2"
      COMPILER_CONTEXT_CALC_R="64"
      COMPILER_CONTEXT_CALC_M="3"
      COMPILER_CONTEXT_EXTENSION_DEGREE="$(poly_degree_from_coeff_list "$RING_F_64")"
      COMPILER_CONTEXT_SCALAR_MODULUS="$RING2P16_MOD"
      ;;
    ring-gr-2p16-128-ext)
      COMPILER_CONTEXT_CALC_P="2"
      COMPILER_CONTEXT_CALC_R="128"
      COMPILER_CONTEXT_CALC_M="2"
      COMPILER_CONTEXT_EXTENSION_DEGREE="$(poly_degree_from_coeff_list "$RING_F_128")"
      COMPILER_CONTEXT_SCALAR_MODULUS="$RING2P16_MOD"
      ;;
    ring-gr-2p2-64-ext)
      COMPILER_CONTEXT_CALC_P="2"
      COMPILER_CONTEXT_CALC_R="64"
      COMPILER_CONTEXT_CALC_M="3"
      COMPILER_CONTEXT_EXTENSION_DEGREE="$(poly_degree_from_coeff_list "$RING_F_64")"
      COMPILER_CONTEXT_SCALAR_MODULUS="$RING2P2_MOD"
      ;;
    ring-gr-2p2-128-ext)
      COMPILER_CONTEXT_CALC_P="2"
      COMPILER_CONTEXT_CALC_R="128"
      COMPILER_CONTEXT_CALC_M="2"
      COMPILER_CONTEXT_EXTENSION_DEGREE="$(poly_degree_from_coeff_list "$RING_F_128")"
      COMPILER_CONTEXT_SCALAR_MODULUS="$RING2P2_MOD"
      ;;
    *)
      echo "Compiler suites only support ring contexts, got context_id=$context_id" >&2
      exit 1
      ;;
  esac
}

run_one_context_compiler_eval_ell() {
  local context_id="$1"
  local context_label="$2"
  local compiler_family="$3"
  local compiler_ell="$4"
  shift 4
  local -a bench_args=("$@")

  local compiler_kappa="$COMPILER_KAPPA"
  local d=$((compiler_ell - compiler_kappa))
  local poly_dim=$((1 << compiler_ell))
  local gamma="NA"
  local release_queries="NA"
  local release_query_status="ok"
  local release_query_error="-"
  local precheck_status="ok"
  local precheck_error="-"
  local prefix="${context_id}_ell${compiler_ell}"

  resolve_compiler_context_metadata "$context_id"

  echo "  - ${context_label}: ell=${compiler_ell} (d=${d}, poly_dim=2^${compiler_ell}=${poly_dim})"

  IFS=',' read -r precheck_status precheck_error < <(compiler_precheck "$compiler_family" "$COMPILER_CONTEXT_EXTENSION_DEGREE" "$compiler_kappa" "$COMPILER_CONTEXT_SCALAR_MODULUS")
  if [[ "$precheck_status" == "ok" ]]; then
    local calc_log="$OUT_DIR/logs/${prefix}_calc_iopp.log"
    if ! run_and_log "$calc_log" \
        "$BUILD_DIR/calc_iopp_params" \
        --d "$d" --c "$C" --k0 "$K0" --lambda "$LAMBDA" \
        --p "$COMPILER_CONTEXT_CALC_P" --r "$COMPILER_CONTEXT_CALC_R" --m "$COMPILER_CONTEXT_CALC_M" \
        --auto-gamma; then
      release_query_status="calc_failed"
      release_query_error="$(first_error_line "$calc_log")"
    else
      gamma="$(parse_first "^  gamma    = " 3 "$calc_log")"
      [[ -n "$gamma" ]] || gamma="NA"
      release_queries="$(parse_first "^  l_min_for_IOPP" 3 "$calc_log")"
      if [[ -z "$release_queries" || "$release_queries" == "N/A" || ! "$release_queries" =~ ^[0-9]+$ ]]; then
        release_query_status="no_feasible_queries"
        release_query_error="l_min_for_IOPP unavailable"
        release_queries="NA"
      fi
    fi
  fi

  run_compiler_eval_row "$compiler_family" "$context_id" "$context_label" "ring" "$d" "$poly_dim" "$gamma" "$release_queries" "$release_query_status" "$release_query_error" "$compiler_ell" "$d" "$compiler_kappa" "$COMPILER_CONTEXT_EXTENSION_DEGREE" "$COMPILER_CONTEXT_SCALAR_MODULUS" "${bench_args[@]}"
}

echo "[1/4] Configure Release build in: $BUILD_DIR"
if [[ -n "$EFFECTIVE_CPU_SET" && "$PIN_BUILD" == "1" ]]; then
  taskset -c "$EFFECTIVE_CPU_SET" cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
else
  cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
fi

BUILD_TARGETS=()
if (( RUN_BASEFOLD_RELEASE || RUN_COMPILER_EVAL )); then
  BUILD_TARGETS+=(calc_iopp_params)
fi
if (( RUN_BASEFOLD_RELEASE )); then
  BUILD_TARGETS+=(bench_basefold_pcs_commit bench_basefold_pcs_eval)
fi
if (( RUN_COMPILER_EVAL )); then
  case "$COMPILER_EVAL_FAMILY" in
    ring_switch)
      BUILD_TARGETS+=(bench_z2k_ring_switch_eval)
      ;;
    frobenius)
      BUILD_TARGETS+=(bench_z2k_frobenius_eval)
      ;;
    *)
      echo "Unknown compiler eval family: $COMPILER_EVAL_FAMILY" >&2
      exit 1
      ;;
  esac
fi

echo "[2/4] Build required benchmarks/tools"
if [[ -n "$EFFECTIVE_CPU_SET" && "$PIN_BUILD" == "1" ]]; then
  taskset -c "$EFFECTIVE_CPU_SET" cmake --build "$BUILD_DIR" --target "${BUILD_TARGETS[@]}" --parallel
else
  cmake --build "$BUILD_DIR" --target "${BUILD_TARGETS[@]}" --parallel
fi

echo "[3/4] Run sweep: suite=$RUN_SUITE, selected contexts = $SELECTED_CONTEXT_COUNT"
if (( RUN_BASEFOLD_RELEASE )); then
  echo "  basefold_release d in [$D_MIN, $D_MAX]"
fi
if (( RUN_COMPILER_EVAL )); then
  echo "  compiler_eval family=$COMPILER_EVAL_FAMILY ell in [$COMPILER_ELL_MIN, $COMPILER_ELL_MAX] with kappa=$COMPILER_KAPPA (d=ell-kappa)"
fi

if (( RUN_BASEFOLD_RELEASE )); then
  for d in $(seq "$D_MIN" "$D_MAX"); do
    echo "[basefold d=$d]"
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
fi

if (( RUN_COMPILER_EVAL )); then
  for compiler_ell in $(seq "$COMPILER_ELL_MIN" "$COMPILER_ELL_MAX"); do
    echo "[compiler_eval/$COMPILER_EVAL_FAMILY ell=$compiler_ell d=$((compiler_ell - COMPILER_KAPPA))]"
    if (( ENABLE_RING2P16_162 )); then
      run_one_context_compiler_eval_ell \
        "ring-gr-2p16-162" "GR(2^16,162)" "$COMPILER_EVAL_FAMILY" "$compiler_ell" \
        --ring-mod "$RING2P16_MOD" \
        --ring-p "$RING2P16_P" \
        --ring-F "$RING2P16_F" \
        --ring-zeta "$RING2P16_ZETA"
    fi
    if (( ENABLE_RING2P2_162 )); then
      run_one_context_compiler_eval_ell \
        "ring-gr-2p2-162" "GR(2^2,162)" "$COMPILER_EVAL_FAMILY" "$compiler_ell" \
        --ring-mod "$RING2P2_MOD" \
        --ring-p "$RING2P2_P" \
        --ring-F "$RING2P2_F" \
        --ring-zeta "$RING2P2_ZETA"
    fi
    if (( ENABLE_RING2P16_64_EXT )); then
      run_one_context_compiler_eval_ell \
        "ring-gr-2p16-64-ext" "GR(2^16,64) (ext-challenge)" "$COMPILER_EVAL_FAMILY" "$compiler_ell" \
        --ring-mod "$RING2P16_MOD" \
        --ring-p "$RING2P16_P" \
        --ring-F "$RING_F_64" \
        --ring-zeta "$RING2P16_ZETA"
    fi
    if (( ENABLE_RING2P16_128_EXT )); then
      run_one_context_compiler_eval_ell \
        "ring-gr-2p16-128-ext" "GR(2^16,128) (ext-challenge)" "$COMPILER_EVAL_FAMILY" "$compiler_ell" \
        --ring-mod "$RING2P16_MOD" \
        --ring-p "$RING2P16_P" \
        --ring-F "$RING_F_128" \
        --ring-zeta "$RING2P16_ZETA"
    fi
    if (( ENABLE_RING2P2_64_EXT )); then
      run_one_context_compiler_eval_ell \
        "ring-gr-2p2-64-ext" "GR(2^2,64) (ext-challenge)" "$COMPILER_EVAL_FAMILY" "$compiler_ell" \
        --ring-mod "$RING2P2_MOD" \
        --ring-p "$RING2P2_P" \
        --ring-F "$RING_F_64" \
        --ring-zeta "$RING2P2_ZETA"
    fi
    if (( ENABLE_RING2P2_128_EXT )); then
      run_one_context_compiler_eval_ell \
        "ring-gr-2p2-128-ext" "GR(2^2,128) (ext-challenge)" "$COMPILER_EVAL_FAMILY" "$compiler_ell" \
        --ring-mod "$RING2P2_MOD" \
        --ring-p "$RING2P2_P" \
        --ring-F "$RING_F_128" \
        --ring-zeta "$RING2P2_ZETA"
    fi
  done
fi

echo "[4/4] Build markdown summary"
{
  echo "# Release Sweep Results: c=$C, lambda=$LAMBDA"
  echo ""
  echo "- run_at_utc: $RUN_AT_UTC"
  echo "- run_id: $RUN_ID"
  echo "- suite: $RUN_SUITE"
  echo "- build_dir: $BUILD_DIR"
  echo "- isolate_build_dir: $ISOLATE_BUILD_DIR"
  echo "- output_dir: $OUT_DIR"
  if (( RUN_BASEFOLD_RELEASE )); then
    echo "- basefold_d_range: [$D_MIN, $D_MAX] (backend_eval poly_dim = ${K0}*2^d)"
  else
    echo "- basefold_d_range: n/a"
  fi
  if (( RUN_COMPILER_EVAL )); then
    echo "- compiler_eval_family: $COMPILER_EVAL_FAMILY"
    echo "- compiler_eval_kappa: $COMPILER_KAPPA"
    echo "- compiler_eval_ell_range: [$COMPILER_ELL_MIN, $COMPILER_ELL_MAX] (compiler_eval/${COMPILER_EVAL_FAMILY} rows use d=ell-kappa and poly_dim=2^ell)"
  else
    echo "- compiler_eval_family: n/a"
    echo "- compiler_eval_kappa: n/a"
    echo "- compiler_eval_ell_range: n/a"
  fi
  echo "- contexts: $CONTEXTS"
  echo "- bench_threads: $BENCH_THREADS (set 0 to use runtime default)"
  echo "- commit_warmup/reps: $COMMIT_WARMUP / $COMMIT_REPS"
  echo "- eval_warmup/reps: $EVAL_WARMUP / $EVAL_REPS"
  echo "- basefold_release_source: $BASEFOLD_RELEASE_SOURCE"
  echo "- release_query_source: $RELEASE_QUERY_SOURCE"
  echo "- compiler_eval_source: $COMPILER_EVAL_SOURCE"
  echo "- backend_eval_source: $BACKEND_EVAL_SOURCE"
  echo "- compiler_eval_layout_source: $COMPILER_EVAL_LAYOUT_SOURCE"
  echo "- compiler_backend_source: $COMPILER_BACKEND_SOURCE"
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
  echo "## Backend Eval"
  echo ""
  echo "| family | context | d | poly_dim | gamma | queries | commit ms | prove ms | verify ms | proof KB | status |"
  echo "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"
  tail -n +2 "$BACKEND_EVAL_RESULT_CSV" | awk -F',' '{printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", $1, $3, $5, $6, $10, $11, $12, $13, $14, $15, $17}'
  echo ""
  echo "## Compiler Eval"
  echo ""
  echo "| family | context | d | poly_dim | queries | ell | ell' | kappa | outer commit ms | backend commit ms | commit total ms | prove total ms | verify total ms | outer proof KB | total proof KB | status |"
  echo "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |"
  tail -n +2 "$COMPILER_EVAL_RESULT_CSV" | awk -F',' '{printf "| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", $1, $3, $5, $6, $11, $12, $13, $14, $15, $16, $17, $18, $21, $24, $26, $28}'
  echo ""
  echo "## Status Counts"
  echo ""
  echo "### Backend Eval"
  tail -n +2 "$BACKEND_EVAL_RESULT_CSV" | awk -F',' '{cnt[$17]++} END {for (k in cnt) printf "- %s: %d\n", k, cnt[k]}'
  echo ""
  echo "### Compiler Eval"
  tail -n +2 "$COMPILER_EVAL_RESULT_CSV" | awk -F',' '{cnt[$28]++} END {for (k in cnt) printf "- %s: %d\n", k, cnt[k]}'
  echo ""
  echo "Backend eval csv: $BACKEND_EVAL_RESULT_CSV"
  echo "Compiler eval csv: $COMPILER_EVAL_RESULT_CSV"
  echo "Raw logs dir: $OUT_DIR/logs"
} > "$RESULT_MD"

echo "Done."
echo "Summary markdown:   $RESULT_MD"
echo "Backend eval csv:   $BACKEND_EVAL_RESULT_CSV"
echo "Compiler eval csv:  $COMPILER_EVAL_RESULT_CSV"
