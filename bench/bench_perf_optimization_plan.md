# Bench Performance Optimization Plan

## Goal

This plan targets only the benchmark surface under `bench/` and keeps the core protocol implementation aligned with the paper semantics.

The optimization objective is:

- Make the six benchmark binaries measure the intended hot path rather than setup, validation, or correctness cross-check work.
- Move all clearly precomputable work into setup or one-time pre-benchmark preparation.
- Keep checked/general APIs intact, and add or use unchecked / precomputed fast paths only where they serve benchmark measurement.
- Stop at moderate, high-leverage optimization; do not turn this prototype repo into a deployment-grade performance engineering project.

Target benchmarks:

- `bench/bench_basefold_pcs_eval.cpp`
- `bench/bench_basefold_pcs_commit.cpp`
- `bench/bench_z2k_ring_switch_eval.cpp`
- `bench/bench_z2k_ring_switch_commit.cpp`
- `bench/bench_z2k_frobenius_eval.cpp`
- `bench/bench_z2k_frobenius_commit.cpp`

## Non-Negotiable Constraints

1. Benchmark files should not do protocol correctness checking, commitment cross-checking, or repeated parameter validation in the measured workflow. The benchmark assumes honest inputs and correctly prepared parameters.
2. Setup-time work must not be charged to commit/eval benchmarks. Any significant precomputable work should be moved into setup or pre-loop preparation.
3. Protocol semantics and checked public APIs must remain available. Optimization should come from fast paths, unchecked entry points, and cached/precomputed data, not from weakening the generic implementation contract.
4. The optimization target is benchmark fidelity and moderate speedup, not maximal micro-optimization.
5. If a remaining hotspot is protocol-inherent or requires an invasive semantics-distorting refactor, stop.

Practical interpretation of item 1:

- Remove protocol-level validation, commitment/proof consistency checks, and honest-witness cross-checks from the six benchmark drivers.
- Keep only minimal CLI parsing / fatal input-shape handling that is needed for the binary to run at all. These are not performance targets and are outside timed regions.

## Current Confirmed State

Already landed:

- `bench_basefold_pcs_eval` already uses unchecked BaseFold commit artifacts and unchecked prove paths.
- `bench_basefold_pcs_commit` is already documented as an unchecked top-commit benchmark.
- `bench_z2k_ring_switch_eval` already uses unchecked prove and unchecked verify hot paths.
- `bench_z2k_frobenius_eval` now uses unchecked prove and unchecked verify hot paths.
- The compiler backend used by ring-switch now defaults to the unchecked BaseFold commit/prove hot path.
- Ring-switch hot paths already removed repeated internal params/basis validation and no longer build `r_table -> r_monomial_coeffs` just to evaluate `g*`.

Still misaligned with the desired benchmark contract:

- Ring-switch still has two major implementation hotspots after the previous cleanup:
  - `RecoverPartialEvaluationsFromSByU`
  - `BuildRingSwitchComponentTensorInternal`
- Ring-switch packing still uses the generic checked compose path, which is costly for commit-oriented benchmarks.

## Benchmark Timing Contract

The benchmark contract after this optimization pass should be:

Untimed:

- CLI parsing
- ring / field context initialization
- parameter setup
- basis discovery / normalization
- all reusable precomputed tables or linear transforms
- deterministic message / point generation
- claimed value generation
- one-time benchmark fixture construction

Timed in commit benchmarks:

- The headline commit work only
- If the benchmark reports `packing` separately, that packing path must also be hot-path-only and free of repeated validation
- Backend commit timing should exclude setup and debug cross-checks

Timed in eval benchmarks:

- Commit split only if the benchmark explicitly reports commit split
- Prove hot path only
- Verify hot path only
- Outer/backend splits only if the benchmark explicitly reports them

Not present in timed or untimed benchmark loops:

- repeated `Validate*` calls
- repeated checked wrapper calls when an unchecked equivalent exists
- commitment equality cross-checks against a second implementation
- honest-proof sanity checks inside the benchmark driver

## Work Packages

### Phase 0: Bench Contract Cleanup

Status: `[x]`

Tasks:

- `[x]` Remove benchmark-driver correctness / consistency checks from:
  - `bench/bench_basefold_pcs_eval.cpp`
  - `bench/bench_basefold_pcs_commit.cpp`
  - `bench/bench_z2k_ring_switch_eval.cpp`
  - `bench/bench_z2k_ring_switch_commit.cpp`
  - `bench/bench_z2k_frobenius_eval.cpp`
  - `bench/bench_z2k_frobenius_commit.cpp`
- `[x]` Delete compiler-commitment cross-checks in ring-switch and Frobenius commit benches.
- `[x]` Remove benchmark-driver post-run proof/verify sanity checks from the six target benches and rely on `tests/` for correctness coverage.
- `[x]` Update help text / notes so each benchmark explicitly states what is excluded from timing.

Acceptance criteria:

- No benchmark loop contains commitment mismatch checks or benchmark-local proof correctness checks.
- Benchmark files no longer call checked compiler wrappers by default when an unchecked hot path exists.

### Phase 1: Setup / Precompute Boundary Cleanup

Status: `[x]`

Tasks:

- `[x]` Audit the six target benchmarks and move all non-message-dependent preparation outside timed regions.
- `[x]` Ensure setup-heavy work is done once per process or once per benchmark run, not once per repetition.
- `[x]` If a benchmark still needs per-repetition artifact construction only for timing decomposition, construct only the minimum hot-path artifacts being measured.
- `[x]` Keep proof size computation and reporting outside timed regions.

Acceptance criteria:

- No setup-like work remains inside per-repetition timers.
- Repetition loops contain only the operation being measured plus anti-optimization bookkeeping.

### Phase 2: Ring-Switch Core Precompute and Hot-Path Simplification

Status: `[ ]`

Primary targets:

- `include/Compiler/Z2k/RingSwitchPCS.hpp`
- `src/Compiler/Z2k/RingSwitchPCS.cpp`
- `bench/bench_z2k_ring_switch_eval.cpp`
- `bench/bench_z2k_ring_switch_commit.cpp`

Tasks:

- `[ ]` Add a `RingSwitchPCS` setup-owned precomputed block, analogous in spirit to Frobenius precomputed tables.
- `[ ]` Precompute basis metadata that allows fast-path selection at setup time.
- `[ ]` Precompute the reusable linear data needed for:
  - recovering `s_by_u` into partial evaluations
  - packing base-ring coefficients into packed GR evaluations
  - decomposing equality-table values into alpha-basis coordinates
- `[ ]` Replace the current generic `RecoverPartialEvaluationsFromSByU` implementation with:
  - a polynomial-basis fast path for the default basis case
  - a generic precomputed transform path for provided-basis cases
- `[ ]` Replace generic checked packing in `PackZ2kCoeffsToGREvals` with a hot-path version that uses precomputed data and does not revalidate params per call.
- `[ ]` Rewrite `BuildRingSwitchComponentTensorInternal` so it computes `EqualityTableFromPoint(r_suffix)` once and then reuses precomputed decomposition data instead of:
  - rebuilding a Boolean point per `w`
  - calling `EqPolynomial` per `w`
  - doing a fresh generic basis recovery per `w`
- `[ ]` Keep checked APIs intact by layering the fast path under setup-owned trusted params and unchecked entry points.

Acceptance criteria:

- Bench hot paths no longer pay generic basis recovery / composition costs when the default polynomial basis is used.
- `RecoverPartialEvaluationsFromSByU` is no longer the obvious first implementation hotspot for ring-switch at small and medium benchmark sizes.
- `BuildRingSwitchComponentTensorInternal` no longer computes `EqPolynomial` point-by-point.

### Phase 3: Frobenius Bench Parity

Status: `[ ]`

Primary targets:

- `include/Compiler/Z2k/FrobeniusPCS.hpp`
- `src/Compiler/Z2k/FrobeniusPCS.cpp`
- `bench/bench_z2k_frobenius_eval.cpp`
- `bench/bench_z2k_frobenius_commit.cpp`

Tasks:

- `[x]` Add unchecked Frobenius verify entry points and use them by default in `bench_z2k_frobenius_eval.cpp`.
- `[ ]` Ensure Frobenius commit benchmarks do not call checked top-level commit wrappers in the timed path.
- `[ ]` Keep all basis discovery / provided-basis normalization / precomputed-table construction strictly in setup.
- `[ ]` Audit Frobenius eval/commit benchmarks for any remaining benchmark-local correctness cross-checks and remove them.

Acceptance criteria:

- `bench_z2k_frobenius_eval` and `bench_z2k_ring_switch_eval` follow the same benchmark contract shape: unchecked prove + unchecked verify + timed hot path only.
- `bench_z2k_frobenius_commit` measures only packing and backend commit work, not checked wrapper overhead.

### Phase 4: BaseFold Bench Hygiene

Status: `[ ]`

Primary targets:

- `bench/bench_basefold_pcs_eval.cpp`
- `bench/bench_basefold_pcs_commit.cpp`

Tasks:

- `[ ]` Audit the two BaseFold benchmark drivers and strip benchmark-local correctness / validation code that is not needed for the measurement contract.
- `[ ]` Keep current semantics that `bench_basefold_pcs_eval` excludes top encode+commit from prover timing.
- `[ ]` Do not modify BaseFold PCS core unless the audit reveals benchmark-visible checked-wrapper overhead that materially affects the published timing.

Acceptance criteria:

- BaseFold benchmarks remain the measurement baseline.
- BaseFold benchmark code is no stricter than the compiler benchmark code about correctness checking.

### Phase 5: Instrumentation and Stop Rules

Status: `[ ]`

Tasks:

- `[ ]` Add commit-path profiling only if needed to decide between two plausible optimization directions.
- `[ ]` Keep outer-only diagnostic benches as diagnostic tools, not optimization targets.
- `[ ]` After each phase, record a short before/after timing note in this file.
- `[ ]` Stop once the remaining dominant costs are either:
  - backend work
  - protocol-inherent outer work
  - or low-payoff constant-factor cleanup

Acceptance criteria:

- The repo ends with a cleaner benchmark contract and a shorter list of dominant hotspots.
- No optimization step widens scope into unnecessary deployment-grade engineering.

## Per-File Task Board

### Bench Files

- `[x]` `bench/bench_basefold_pcs_eval.cpp`
  - remove benchmark-local correctness checks not required for timing
  - keep current unchecked prove/commit-artifact flow
- `[x]` `bench/bench_basefold_pcs_commit.cpp`
  - audit and strip nonessential benchmark-local validation
- `[x]` `bench/bench_z2k_ring_switch_eval.cpp`
  - keep unchecked prove/verify
  - remove remaining benchmark-local correctness checks
  - ensure commit split construction uses hot-path-only wrappers
- `[x]` `bench/bench_z2k_ring_switch_commit.cpp`
  - stop calling checked `RingSwitchPCSCommit` in timed path
  - remove commitment cross-check
- `[x]` `bench/bench_z2k_frobenius_eval.cpp`
  - switch default benchmark verify to unchecked verify
  - remove benchmark-local correctness checks
- `[x]` `bench/bench_z2k_frobenius_commit.cpp`
  - stop calling checked `FrobeniusPCSCommit` in timed path
  - remove commitment cross-check

### Core Files

- `[ ]` `include/Compiler/Z2k/RingSwitchPCS.hpp`
  - add setup-owned precomputed fast-path data
- `[ ]` `src/Compiler/Z2k/RingSwitchPCS.cpp`
  - implement ring-switch setup precompute
  - implement fast packing / partial-recovery / tensor-build paths
- `[x]` `include/Compiler/Z2k/FrobeniusPCS.hpp`
  - declare unchecked verify APIs if needed
- `[x]` `src/Compiler/Z2k/FrobeniusPCS.cpp`
  - implement unchecked verify hot path and keep checked wrappers intact
- `[ ]` `include/PCS/Common/Profile.hpp`
  - touch only if commit-path profiling becomes necessary

## Execution Order

Recommended implementation order:

1. Bench contract cleanup in the six target benchmark files.
2. Add unchecked Frobenius verify and switch Frobenius benches to the same measurement contract as ring-switch.
3. Land ring-switch setup-owned precomputed data and fast-path packing / partial recovery.
4. Land ring-switch tensor-build cleanup based on one-shot equality table generation.
5. Re-audit BaseFold benches and stop unless a clear benchmark-visible overhead remains.

Rationale:

- Steps 1 and 2 are low-risk and immediately align benchmark semantics.
- Step 3 is likely the largest remaining real speed win.
- Step 4 is the second ring-switch hotspot after the previous `g*` cleanup.
- Step 5 prevents unnecessary churn in the BaseFold baseline.

## Validation Plan

Build targets:

```bash
cmake --build build-release --target \
  bench_basefold_pcs_commit \
  bench_basefold_pcs_eval \
  bench_z2k_ring_switch_commit \
  bench_z2k_ring_switch_eval \
  bench_z2k_frobenius_commit \
  bench_z2k_frobenius_eval \
  test_pcs \
  test_z2k_ring_switch_pcs \
  test_z2k_frobenius_pcs \
  test_z2k_ring_switch_bench_cli \
  test_z2k_frobenius_bench_cli -j 4
```

Minimum regression tests after core changes:

```bash
./build-release/test_pcs
./build-release/test_z2k_ring_switch_pcs
./build-release/test_z2k_frobenius_pcs
./build-release/test_z2k_ring_switch_bench_cli
./build-release/test_z2k_frobenius_bench_cli
```

Representative benchmark reruns:

```bash
OMP_NUM_THREADS=1 ./build-release/bench_basefold_pcs_commit ...
OMP_NUM_THREADS=1 ./build-release/bench_basefold_pcs_eval ...
OMP_NUM_THREADS=1 ./build-release/bench_z2k_ring_switch_commit ...
OMP_NUM_THREADS=1 ./build-release/bench_z2k_ring_switch_eval ...
OMP_NUM_THREADS=1 ./build-release/bench_z2k_frobenius_commit ...
OMP_NUM_THREADS=1 ./build-release/bench_z2k_frobenius_eval ...
```

## Before / After Log

Fill this section as work lands.

### Phase 0

- Before:
  - Commit benches still performed benchmark-driver commitment cross-checks.
  - Eval benches still performed benchmark-driver verify-success assertions.
  - Frobenius eval still defaulted to the checked verifier path.
- After:
  - All six target benches now drop benchmark-driver correctness/cross-check logic from the measured workflow.
  - Ring-switch and Frobenius commit benches now time backend commit on prepacked monomial data rather than checked top-level commit wrappers.
  - Frobenius eval now defaults to unchecked verify, matching the release-sweep hot-path contract.
- Notes:
  - Validated with `bench_*` target rebuild plus `test_z2k_frobenius_pcs`, `test_z2k_ring_switch_bench_cli`, and `test_z2k_frobenius_bench_cli`.

### Phase 1

- Before:
- `bench_basefold_pcs_eval` still rebuilt `BaseFoldPCSCommitArtifactsUnchecked(...)` once per repetition even though it depended only on the fixed benchmark fixture.
- Eval benches still recomputed serializer-backed proof sizes on every repetition, even though proof size is reported only once and is not part of any timed metric.
- After:
- `bench_basefold_pcs_eval` now constructs the committed top-oracle artifacts once per benchmark run and reuses them across all measured repetitions.
- `bench_basefold_pcs_eval`, `bench_z2k_ring_switch_eval`, and `bench_z2k_frobenius_eval` now compute proof size once from the last measured proof, after the timed repetitions have finished.
- Compiler eval benches keep per-repetition outer/backend commit construction only where it is the metric being reported, and otherwise reuse loop-local bookkeeping objects.
- Notes:
- No commit-benchmark timing contract changed in this phase; ring-switch and Frobenius commit/eval loops were left with only the explicit split metrics plus anti-optimization bookkeeping.

### Phase 2

- Before:
- After:
- Notes:

### Phase 3

- Before:
- After:
- Notes:

### Phase 4

- Before:
- After:
- Notes:

## Out of Scope

The following are intentionally out of scope for this optimization pass:

- changing protocol semantics to get better benchmark numbers
- removing checked public APIs from the library
- persistent cross-process caches
- deep thread-level / NUMA / allocator optimization
- deployment-grade observability or benchmarking infrastructure
- optimizing outer-only diagnostic benches as primary targets
