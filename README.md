# BasefoldOverGR

For the Chinese version, see [README_zh.md](README_zh.md).

BasefoldOverGR contains C++ implementations of BaseFold-style encoding, IOPP,
and PCS over finite fields and Galois rings, together with `Z_{2^k}` compiler
benchmarks for ring-switch and Frobenius variants.

The main benchmark entrypoint is
[scripts/run_release_c4_lambda128.sh](scripts/run_release_c4_lambda128.sh).
It auto-configures a Release build in `build-release/` and writes fresh sweep
outputs under `results/` unless `OUT_DIR` is set explicitly.

## Build and Test

Requirements:

- CMake 3.16 or newer
- A C++17 compiler
- NTL and GMP
- OpenMP if you want multithreaded benches
- `python3` and `matplotlib` for
  [scripts/plot_benchmark_results.py](scripts/plot_benchmark_results.py)

Manual build and test:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

The sweep script manages its own `build-release/` tree. Local directories such
as `build/`, `build-release/`, `results/`, and `.venv/` are git-ignored and may
exist in a working tree without being part of the tracked repository layout.

## Repository Layout

```text
.
├── CMakeLists.txt
├── README.md
├── README_zh.md
├── LICENSE
├── bench/                      # benchmark binaries and calc_iopp_params
├── include/                    # public headers
├── src/                        # library implementation
├── tests/                      # unit tests and CLI/path hygiene checks
├── scripts/
│   ├── run_release_c4_lambda128.sh
│   ├── run_backend_eval_single_thread_7contexts.sh
│   └── plot_benchmark_results.py
├── FRI_Ligero-based_results.md # tracked FRI/Ligero baseline tables
├── results-new/                # tracked current-format CSVs and figures
└── results-legacy/             # tracked older CSVs and plots
```

Tracked data currently lives in:

- `results-new/1-thread/backend_eval_results.csv`
- `results-new/8-thread/backend_eval_results.csv`
- `results-new/8-thread/compiler_eval_results.csv`
- `results-new/fri_ligero_based_eval_results.csv`
- `results-new/figures/*.png`
- `results-legacy/` for older CSV names and `*_vs_d.png` plots

## Benchmark Workflows

`scripts/run_release_c4_lambda128.sh` is fully environment-variable driven.
Current suites are:

- `basefold_release`
- `compiler_eval_ring_switch`
- `compiler_eval_frobenius`
- `compiler_outer_commit_ring_switch`
- `compiler_outer_commit_frobenius`

Current defaults:

- `C=4`
- `K0=1`
- `LAMBDA=128`
- `D_MIN=3`, `D_MAX=29`
- `COMMIT_WARMUP=1`, `COMMIT_REPS=3`
- `EVAL_WARMUP=1`, `EVAL_REPS=3`
- `BENCH_THREADS=8`
- `CONTEXTS=all` (16 contexts)

Current context ids:

- `field-255`
- `ring-gr-2p16-162`
- `field-f2p256`
- `ring-gr-2p2-162`
- `field-prime64-ext`
- `field-f2p64-ext`
- `field-prime128-ext`
- `field-f2p128-ext`
- `field-f3p40-ext`
- `field-f3p81-ext`
- `ring-gr-2p16-64-ext`
- `ring-gr-2p16-128-ext`
- `ring-gr-2p32-64-ext`
- `ring-gr-2p32-128-ext`
- `ring-gr-2p2-64-ext`
- `ring-gr-2p2-128-ext`

When `BENCH_THREADS > 0`, the runner also exports:

- `OMP_NUM_THREADS`
- `BASEFOLD_MERKLE_MAX_THREADS`
- `BASEFOLD_VERIFY_QUERY_MAX_THREADS`

to the same value. Set `BENCH_THREADS=0` if you want to keep the runtime's
default threading policy.

Representative commands:

```bash
# Default BaseFold release sweep
scripts/run_release_c4_lambda128.sh

# One BaseFold context with a smaller d range
CONTEXTS=field-prime128-ext \
D_MIN=10 D_MAX=20 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh

# Ring-switch compiler full eval on GR(2^32,64)
RUN_SUITE=compiler_eval_ring_switch \
CONTEXTS=ring-gr-2p32-64-ext \
COMPILER_KAPPA=6 \
COMPILER_ELL_MIN=9 \
COMPILER_ELL_MAX=12 \
EVAL_WARMUP=0 EVAL_REPS=1 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh

# Frobenius compiler full eval on GR(2^32,128)
RUN_SUITE=compiler_eval_frobenius \
CONTEXTS=ring-gr-2p32-128-ext \
COMPILER_KAPPA=7 \
COMPILER_ELL_MIN=9 \
COMPILER_ELL_MAX=12 \
EVAL_WARMUP=0 EVAL_REPS=1 \
BENCH_THREADS=1 \
scripts/run_release_c4_lambda128.sh
```

Important current constraints:

- Compiler suites accept ring contexts only.
- The selected ring context must satisfy `deg(F) = 2^COMPILER_KAPPA`.
- Current compiler benches support `K0=1` only.
- Frobenius compiler suites skip `GR(2^2,r)` contexts and emit
  `status=disabled_gr2p2_context`.
- Sweep CSV `queries` are derived from `calc_iopp_params` by parsing
  `l_min_for_IOPP`.

Useful environment variables:

- `CONTEXTS` for selecting a subset of contexts
- `COMPILER_KAPPA`, `COMPILER_ELL_MIN`, `COMPILER_ELL_MAX` for compiler suites
- `OUT_DIR`, `OUT_ROOT`, `RUN_ID`, `TIMESTAMP` for output placement
- `BUILD_DIR`, `ISOLATE_BUILD_DIR`, `PIN_BUILD` for build-tree control
- `CPU_PIN_MODE`, `CPU_SET`, `RUN_SLOT`, `RUN_SLOTS_TOTAL`,
  `USE_SMT_IN_SLOT` for CPU pinning
- `CMD_TIMEOUT_SEC`, `CONTINUE_ON_ERROR` for failure handling

The full per-context parameter table and environment-variable reference are in
[bench/exp_params_release_c4_lambda128.md](bench/exp_params_release_c4_lambda128.md).

### Concurrent Single-Thread Backend Sweep

[scripts/run_backend_eval_single_thread_7contexts.sh](scripts/run_backend_eval_single_thread_7contexts.sh)
launches seven single-thread `basefold_release` runs in parallel, pins each
child to one CPU, and merges the child CSVs into one combined
`backend_eval_results.csv`.

Default example:

```bash
scripts/run_backend_eval_single_thread_7contexts.sh
```

This helper expects `bash`, `lscpu`, `python3`, and `taskset`.

## Output Files

The run directory defaults to `results/release_c4_lambda128_sweep_<RUN_ID>/`.
Depending on `RUN_SUITE`, it can contain:

- `backend_eval_results.csv` for `basefold_release`
- `compiler_eval_results.csv` for `compiler_eval_ring_switch` and
  `compiler_eval_frobenius`
- `compiler_outer_commit_results.csv` for
  `compiler_outer_commit_ring_switch` and
  `compiler_outer_commit_frobenius`
- `RESULTS.md` for markdown summaries
- `logs/*.log` for raw `calc_iopp_params` and benchmark logs

`backend_eval_results.csv` uses the current backend schema:

- `commit_mean_ms`
- `open_mean_ms`
- `prove_mean_ms`
- `verifier_mean_ms`
- `proof_size_kb`
- `proof_size_bytes`

`compiler_eval_results.csv` uses the current compiler schema, including:

- `outer_commit_mean_ms`
- `backend_commit_mean_ms`
- `commit_total_mean_ms`
- `open_total_mean_ms`
- `prove_total_mean_ms`
- `verify_total_mean_ms`
- `outer_proof_size_kb`
- `total_proof_size_kb`

## Plotting CSV Results

`scripts/plot_benchmark_results.py` reads one or more benchmark CSV files and
plots `commit`, `open`, `prover`, `verifier`, and `proof_size` against the
number of constraints.

If you omit `-o`, the script writes PNG files to `result/plots`. Use `-o` when
you want to regenerate the tracked figures under `results-new/figures`.

Examples:

```bash
# Regenerate the tracked single-thread backend figures
python3 scripts/plot_benchmark_results.py \
  results-new/1-thread/backend_eval_results.csv \
  -o results-new/figures \
  --prefix PCSoverGaloisRing_1_thread

# Regenerate the tracked 8-thread compiler figures
python3 scripts/plot_benchmark_results.py \
  results-new/8-thread/compiler_eval_results.csv \
  -o results-new/figures \
  --prefix PCSoverZ2K_8_threads
```

Useful plotting options:

- `--metrics commit prover verifier proof_size`
- `--proof-size-column <column>`
- `--legend-mode auto|inline|split`
- Multi-input mode for overlaying several CSVs in one figure
