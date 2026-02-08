#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build-release}"
OUT_ROOT="${OUT_ROOT:-$ROOT_DIR/results}"
TIMESTAMP="${TIMESTAMP:-$(date +%Y%m%d_%H%M%S)}"
OUT_DIR="${OUT_DIR:-$OUT_ROOT/release_c4_lambda128_sweep_${TIMESTAMP}}"

# Target profile: rate = 1/4 (c=4), security = 128 bits.
C="${C:-4}"
K0="${K0:-1}"
LAMBDA="${LAMBDA:-128}"
D_MIN="${D_MIN:-3}"
D_MAX="${D_MAX:-29}"

COMMIT_WARMUP="${COMMIT_WARMUP:-1}"
COMMIT_REPS="${COMMIT_REPS:-3}"
EVAL_WARMUP="${EVAL_WARMUP:-0}"
EVAL_REPS="${EVAL_REPS:-1}"
SEED="${SEED:-0}"
RUN_PROOF_SIZE="${RUN_PROOF_SIZE:-1}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-1}"
CMD_TIMEOUT_SEC="${CMD_TIMEOUT_SEC:-0}"

FIELD255_MOD="${FIELD255_MOD:-57896044618658097711785492504343953926634992332820282019728792003956564819949}"  # 2^255 - 19
FIELD255_F="${FIELD255_F:-1,1}"       # x + 1
FIELD255_ZETA="${FIELD255_ZETA:-0,1}"  # x

F2_256_MOD="${F2_256_MOD:-2}"
F2_256_F_DEFAULT='1,0,1,0,0,1,0,0,1,0,1,0,0,1,1,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,1,0,1,1,0,0,0,0,0,1,1,1,0,1,1,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,1,1,1,0,1,0,1,1,1,0,0,1,1,0,1,0,1,0,1,0,1,0,0,0,0,0,0,0,1,1,0,1,1,0,1,0,1,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,0,0,0,1,0,1,1,0,1,1,0,1,1,0,1,0,1,1,0,0,1,0,0,0,0,1,1,0,0,0,0,1,0,0,1,0,1,0,0,0,0,0,1,1,1,0,1,1,0,1,1,0,1,1,0,0,0,0,1,0,0,0,0,1,1,0,0,1,1,0,0,1,1,1,0,1,0,0,0,1,1,1,0,1,1,0,0,1,0,1,1,0,1,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,1,1,0,1,1,1,0,1,1,0,0,1,0,1'
F2_256_F="${F2_256_F:-$F2_256_F_DEFAULT}"
F2_256_ZETA="${F2_256_ZETA:-0,1}"  # x

RING_F_162_DEFAULT='1,1,0,0,1,1,1,1,0,0,0,1,0,0,1,0,1,1,1,0,0,0,0,0,0,1,1,1,1,0,0,0,1,1,1,0,0,0,1,1,0,0,1,1,1,0,0,0,0,1,1,1,1,0,1,1,0,1,0,0,0,0,1,1,1,0,1,1,1,0,0,1,0,1,1,0,0,1,0,0,1,0,1,0,0,1,1,1,0,1,1,0,1,0,1,1,1,0,1,1,0,1,1,0,0,1,1,1,0,1,0,1,1,0,1,0,0,0,1,1,1,1,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,1,1,1,1,0,0,1,0,1,1,0,0,1,1,1,1,0,0,1,1,1,1,0,1,1,1'

RING2P16_MOD="${RING2P16_MOD:-65536}"  # 2^16
RING2P16_P="${RING2P16_P:-2}"
RING2P16_F="${RING2P16_F:-$RING_F_162_DEFAULT}"
RING2P16_ZETA="${RING2P16_ZETA:-0,1}"  # x

RING2P2_MOD="${RING2P2_MOD:-4}"  # 2^2
RING2P2_P="${RING2P2_P:-2}"
RING2P2_F="${RING2P2_F:-$RING_F_162_DEFAULT}"
RING2P2_ZETA="${RING2P2_ZETA:-0,1}"  # x

if (( D_MIN < 0 || D_MAX < D_MIN )); then
  echo "Invalid dimension range: D_MIN=$D_MIN D_MAX=$D_MAX" >&2
  exit 2
fi
if [[ "$RUN_PROOF_SIZE" != "0" && "$RUN_PROOF_SIZE" != "1" ]]; then
  echo "RUN_PROOF_SIZE must be 0 or 1" >&2
  exit 2
fi
if [[ "$CONTINUE_ON_ERROR" != "0" && "$CONTINUE_ON_ERROR" != "1" ]]; then
  echo "CONTINUE_ON_ERROR must be 0 or 1" >&2
  exit 2
fi

mkdir -p "$OUT_DIR/logs"

RESULT_CSV="$OUT_DIR/results.csv"
RESULT_MD="$OUT_DIR/RESULTS.md"
RUN_AT_UTC="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

cat > "$RESULT_CSV" <<CSV
context_id,context_label,mode,d,poly_dim,c,k0,lambda,gamma,queries,commit_mean_ms,prover_mean_ms,verifier_mean_ms,proof_size_kb,proof_size_bytes,status,error
CSV

run_and_log() {
  local log_file="$1"
  shift
  {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
    if (( CMD_TIMEOUT_SEC > 0 )); then
      timeout "$CMD_TIMEOUT_SEC" "$@"
    else
      "$@"
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
  local prover_ms="$9"
  local verifier_ms="${10}"
  local proof_kb="${11}"
  local proof_bytes="${12}"
  local status="${13}"
  local error="${14}"
  local context_label_csv="${context_label//,/;}"
  local error_csv="${error//,/;}"

  echo "${context_id},${context_label_csv},${mode},${d},${poly_dim},${C},${K0},${LAMBDA},${gamma},${queries},${commit_ms},${prover_ms},${verifier_ms},${proof_kb},${proof_bytes},${status},${error_csv}" >> "$RESULT_CSV"
}

run_one_context_d() {
  local context_id="$1"
  local context_label="$2"
  local mode="$3"
  local d="$4"
  shift 4
  local -a bench_args=("$@")

  local poly_dim=$((1 << d))
  local status="ok"
  local error=""
  local gamma="NA"
  local queries="NA"
  local commit_ms="NA"
  local prover_ms="NA"
  local verifier_ms="NA"
  local proof_kb="NA"
  local proof_bytes="NA"
  local calc_p=""
  local calc_r=""
  local calc_m="1"

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
    *)
      status="context_config_error"
      error="unknown context"
      write_row "$context_id" "$context_label" "$mode" "$d" "$poly_dim" "$gamma" "$queries" "$commit_ms" "$prover_ms" "$verifier_ms" "$proof_kb" "$proof_bytes" "$status" "$error"
      maybe_abort "$status" "unknown context_id=$context_id"
      return
      ;;
  esac

  local prefix="${context_id}_d${d}"
  local calc_log="$OUT_DIR/logs/${prefix}_calc_iopp.log"
  local commit_log="$OUT_DIR/logs/${prefix}_commit.log"
  local eval_log="$OUT_DIR/logs/${prefix}_eval.log"
  local proof_log="$OUT_DIR/logs/${prefix}_proof_size.log"

  echo "  - ${context_label}: d=${d} (poly_dim=2^${d})"

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
    commit_ms="$(parse_first "encode-only mean" 3 "$commit_log")"
    if [[ -z "$commit_ms" ]]; then
      commit_ms="NA"
    fi
  fi

  if [[ "$status" == "ok" ]]; then
    if ! run_and_log "$eval_log" \
        "$BUILD_DIR/bench_pcs_eval" \
        --mode "$mode" \
        "${bench_args[@]}" \
        --c "$C" --k0 "$K0" --d "$d" \
        --queries "$queries" \
        --warmup "$EVAL_WARMUP" --reps "$EVAL_REPS" \
        --seed "$SEED"; then
      status="eval_failed"
      error="$(first_error_line "$eval_log")"
    else
      prover_ms="$(parse_first "prover   mean" 3 "$eval_log")"
      verifier_ms="$(parse_first "verifier mean" 3 "$eval_log")"
      if [[ -z "$prover_ms" ]]; then
        prover_ms="NA"
      fi
      if [[ -z "$verifier_ms" ]]; then
        verifier_ms="NA"
      fi
    fi
  fi

  if [[ "$status" == "ok" && "$RUN_PROOF_SIZE" == "1" ]]; then
    if ! run_and_log "$proof_log" \
        "$BUILD_DIR/bench_pcs_proof_size" \
        --mode "$mode" \
        "${bench_args[@]}" \
        --c "$C" --k0 "$K0" --d "$d" \
        --queries "$queries" \
        --formula \
        --seed "$SEED"; then
      status="proof_failed"
      error="$(first_error_line "$proof_log")"
    else
      proof_kb="$(parse_first "proof size" 3 "$proof_log")"
      proof_bytes="$(awk '/proof size/{gsub(/[^0-9]/, "", $5); print $5; exit}' "$proof_log")"
      if [[ -z "$proof_kb" ]]; then
        proof_kb="NA"
      fi
      if [[ -z "$proof_bytes" ]]; then
        proof_bytes="NA"
      fi
    fi
  fi

  if [[ -z "$error" ]]; then
    error="-"
  fi
  write_row "$context_id" "$context_label" "$mode" "$d" "$poly_dim" "$gamma" "$queries" "$commit_ms" "$prover_ms" "$verifier_ms" "$proof_kb" "$proof_bytes" "$status" "$error"
  maybe_abort "$status" "${context_label} d=${d} status=${status}"
}

echo "[1/4] Configure Release build in: $BUILD_DIR"
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release

echo "[2/4] Build required benchmarks/tools"
cmake --build "$BUILD_DIR" \
  --target bench_pcs_commit bench_pcs_eval bench_pcs_proof_size calc_iopp_params \
  --parallel

echo "[3/4] Run sweep: d in [$D_MIN, $D_MAX], contexts = 4"
for d in $(seq "$D_MIN" "$D_MAX"); do
  echo "[d=$d]"
  run_one_context_d \
    "field-255" "Field-255" "field" "$d" \
    --field-mod "$FIELD255_MOD" \
    --field-F "$FIELD255_F" \
    --field-zeta "$FIELD255_ZETA"
  run_one_context_d \
    "ring-gr-2p16-162" "GR(2^16,162)" "ring" "$d" \
    --ring-mod "$RING2P16_MOD" \
    --ring-p "$RING2P16_P" \
    --ring-F "$RING2P16_F" \
    --ring-zeta "$RING2P16_ZETA"
  run_one_context_d \
    "field-f2p256" "F_2^256" "field" "$d" \
    --field-mod "$F2_256_MOD" \
    --field-F "$F2_256_F" \
    --field-zeta "$F2_256_ZETA"
  run_one_context_d \
    "ring-gr-2p2-162" "GR(2^2,162)" "ring" "$d" \
    --ring-mod "$RING2P2_MOD" \
    --ring-p "$RING2P2_P" \
    --ring-F "$RING2P2_F" \
    --ring-zeta "$RING2P2_ZETA"
done

echo "[4/4] Build markdown summary"
{
  echo "# Release Sweep Results: c=$C, lambda=$LAMBDA"
  echo ""
  echo "- run_at_utc: $RUN_AT_UTC"
  echo "- build_dir: $BUILD_DIR"
  echo "- output_dir: $OUT_DIR"
  echo "- d_range: [$D_MIN, $D_MAX] (poly_dim = 2^d)"
  echo "- run_proof_size: $RUN_PROOF_SIZE"
  echo "- continue_on_error: $CONTINUE_ON_ERROR"
  echo "- cmd_timeout_sec: $CMD_TIMEOUT_SEC"
  echo ""
  echo "## Results"
  echo ""
  echo "| context | d | poly_dim | gamma | queries | commit mean ms | prover mean ms | verifier mean ms | proof size KB | status |"
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
