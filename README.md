# BasefoldOverGR

C++ implementations over finite fields and Galois rings for BaseFold-style encoding, IOPP/PCS, and `Z_{2^k}` compiler benchmarks.

The benchmark workflow in this repo is centered on [scripts/run_release_c4_lambda128.sh](scripts/run_release_c4_lambda128.sh).

## Repository Structure

```text
.
├── CMakeLists.txt
├── README.md
├── README_zh.md
├── LICENSE
├── bench/
│   ├── bench_basefold_pcs_*.cpp
│   ├── bench_z2k_ring_switch_*.cpp
│   ├── bench_z2k_frobenius_*.cpp
│   ├── calc_iopp_params.cpp
│   └── exp_params_release_c4_lambda128.md
├── include/
├── src/
├── tests/
├── scripts/
│   ├── run_release_c4_lambda128.sh
│   └── plot_benchmark_results.py
├── results-legacy/
├── build/             # generated
├── build-release/     # generated
└── results/           # generated benchmark outputs
```

## Running Benchmarks

`scripts/run_release_c4_lambda128.sh` is configured entirely through environment variables. It configures/builds `build-release` automatically and writes outputs to `results/release_c4_lambda128_sweep_<RUN_ID>/` unless `OUT_DIR` is set explicitly.

### Default BaseFold Release Sweep

```bash
scripts/run_release_c4_lambda128.sh
```

Current defaults:

- `RUN_SUITE=basefold_release`
- `CONTEXTS=all`
- `D_MIN=3`
- `D_MAX=29`
- `C=4`
- `K0=1`
- `LAMBDA=128`
- `COMMIT_WARMUP=1`, `COMMIT_REPS=3`
- `EVAL_WARMUP=1`, `EVAL_REPS=3` (shared eval-loop warmup/reps for prove + verify)
- `BENCH_THREADS=8`

### Run One BaseFold Context with a Smaller `d` Range

```bash
CONTEXTS=field-prime128-ext \
D_MIN=10 D_MAX=20 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh
```

### Run Ring-Switch Compiler Eval

```bash
RUN_SUITE=compiler_eval_ring_switch \
CONTEXTS=ring-gr-2p16-64-ext \
COMPILER_KAPPA=6 \
COMPILER_ELL_MIN=9 \
COMPILER_ELL_MAX=12 \
EVAL_WARMUP=0 EVAL_REPS=1 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh
```

Notes:

- Compiler suites currently support ring contexts only.
- The selected ring context must satisfy `deg(F) = 2^COMPILER_KAPPA`.
- Compiler benches currently build only `K0=1` backends; `K0 != 1` rows are emitted as `status=unsupported_k0`.
- In `bench_z2k_ring_switch_eval`, default verifier split is measured as standalone outer replay plus standalone backend-only verify; `--profiled-backend-verify` restores the old subcall-timed split.

### Run Frobenius Compiler Eval

```bash
RUN_SUITE=compiler_eval_frobenius \
CONTEXTS=ring-gr-2p16-64-ext \
COMPILER_KAPPA=6 \
COMPILER_ELL_MIN=9 \
COMPILER_ELL_MAX=12 \
EVAL_WARMUP=0 EVAL_REPS=1 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh
```

Notes:

- Compiler suites currently support ring contexts only.
- The selected ring context must satisfy `deg(F) = 2^COMPILER_KAPPA`.
- Compiler benches currently build only `K0=1` backends; `K0 != 1` rows are emitted as `status=unsupported_k0`.
- `RUN_SUITE=compiler_eval_frobenius` still skips `GR(2^2,r)` contexts and emits `status=disabled_gr2p2_context`.
- In `bench_z2k_frobenius_eval`, default verifier split is measured as standalone outer replay plus standalone backend-only verify; `--profiled-backend-verify` restores the old subcall-timed split.

### Most Useful Environment Variables

- `RUN_SUITE`: `basefold_release`, `compiler_eval_ring_switch`, `compiler_eval_frobenius`
- `CONTEXTS`: `all` or a comma-separated subset of supported contexts
- `D_MIN`, `D_MAX`: BaseFold sweep range, used by `basefold_release`
- `COMPILER_KAPPA`, `COMPILER_ELL_MIN`, `COMPILER_ELL_MAX`: required by compiler suites
- Compiler suites accept ring contexts only and require `deg(F) = 2^COMPILER_KAPPA`
- `K0`: BaseFold message base dimension, default `1`
- Compiler suites currently require `K0=1`
- `BENCH_THREADS`: threads used inside each bench process
- `OUT_DIR`: explicit output directory
- `BUILD_DIR`: explicit build directory
- `ISOLATE_BUILD_DIR=1`: use `build-release-<RUN_ID>` instead of the shared `build-release`
- `CMD_TIMEOUT_SEC`: per-benchmark timeout, default `0` means no timeout
- `CONTINUE_ON_ERROR`: keep sweeping after per-point failures, default `1`

The full environment-variable reference and supported context list are in [bench/exp_params_release_c4_lambda128.md](bench/exp_params_release_c4_lambda128.md).

## Output Files

Each run writes:

- `backend_eval_results.csv`: BaseFold release rows from `bench_basefold_pcs_commit` and `bench_basefold_pcs_eval`
- `compiler_eval_results.csv`: compiler-eval rows for the selected family (`ring_switch` or `frobenius`), including `outer_proof_size_*` and `total_proof_size_*`
- `RESULTS.md`: markdown summary table
- `logs/*.log`: raw logs for `calc_iopp_params` and benchmark binaries

Note:

- `RUN_SUITE=compiler_eval_frobenius` does not execute Frobenius full eval on `GR(2^2,r)` contexts; those rows are emitted as `status=disabled_gr2p2_context`.
