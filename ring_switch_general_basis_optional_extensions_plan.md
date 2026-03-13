# Ring-Switch General Basis Optional Extensions Plan

Status: WP0-WP4 completed; WP5-WP6 pending
Last updated: 2026-03-14
Scope: finish the two follow-on optional extensions after `ring_switch_general_basis_shared_layer_plan.md`

## Progress

- 2026-03-14: `WP0` landed with this plan file and froze the follow-up scope.
- 2026-03-14: `WP1` landed in the working tree.
- Ring-switch bench common helpers now own:
  - the shared basis CLI data model,
  - basis-mode parsing helpers,
  - semicolon-separated basis-element parsing,
  - basis-shape validation against `deg(F)`,
  - a shared setup-input builder for default vs provided basis mode,
  - a shared basis-mode help-text fragment.
- `bench_z2k_ring_switch_commit.cpp` and `bench_z2k_ring_switch_eval.cpp` no longer carry their own private copies of the ring-switch common helper layer.
- All five ring-switch bench binaries now consume the shared basis-mode note fragment from the common helper.
- 2026-03-14: `WP2` landed in the working tree.
- All five ring-switch bench binaries now accept:
  - `--basis-mode default|provided`
  - `--alpha-basis` / `--beta-basis`
  - optional `--alpha-dual-basis` / `--beta-dual-basis`
- Basis CLI decoding now happens only after the active `ZZ_p` / `ZZ_pE`
  contexts are initialized, so provided-basis parsing is aligned with the live
  ring context instead of the pre-context argument loop.
- Benchmark output now prints `basis mode ...`; in provided mode it also prints
  basis dimensions plus whether each dual basis was supplied or derived.
- README now documents the provided-basis bench surface and a minimal example.
- 2026-03-14: `WP3` landed in the working tree.
- `tests/test_z2k_ring_switch_pcs.cpp` now has a reusable ring-switch context
  spec layer for:
  - default `GR(4,2)`,
  - the existing `GR(4,2)` backend variant,
  - `GR(4,4)`.
- Provided-basis test helpers now synthesize bases from deterministic
  unitriangular or dense base-ring change-of-basis matrices instead of
  hand-coded one-off `GR(4,2)` vectors.
- Deterministic helper builders now cover:
  - extension polynomials,
  - backend params from context specs,
  - witness tables,
  - query points,
  - valid transformed bases,
  - singular bases,
  - malformed dual bases.
- Existing provided-basis setup / packing / semantic / end-to-end / proof-size
  tests now consume the shared helper layer instead of the previous ad hoc
  `BuildNonPolynomial*GR42` helpers.
- 2026-03-14: `WP4` landed in the working tree.
- Positive end-to-end provided-basis coverage now includes:
  - Context A `GR(4,2)`:
    - non-polynomial `alpha` only,
    - non-polynomial `beta` only,
    - non-polynomial `alpha != beta`,
    - explicit-dual `alpha != beta`,
    - auto-derived-dual `alpha != beta`.
  - Context B `GR(4,2)` variant backend:
    - one non-polynomial `alpha != beta` end-to-end case.
  - Context C `GR(4,4)`:
    - non-polynomial `alpha` only,
    - non-polynomial `beta` only,
    - non-polynomial `alpha != beta`.
- For Context A and Context C, the positive matrix now also locks:
  - packing against the provided `beta` basis,
  - `A_{u||w}` reconstruction against the provided `alpha` basis,
  - partial-evaluation recovery against direct evaluation,
  - composed proof size against serialized byte length.

## Goal

Land two bounded follow-up improvements for the ring-switch general-basis work:

- expose caller-provided `alpha/beta` basis data on all ring-switch bench CLIs
- raise provided-basis regression coverage from the current `GR(4,2)` witness pair to a high-confidence, deterministic coverage matrix

The plan assumes the current library-level protocol semantics are already correct and does not reopen the ring-switch proof logic itself.

## Non-Goals

- no changes to the ring-switch proof message shape
- no changes to the backend PCS abstraction boundary
- no multipoint opening work
- no Blaze-Orion composition work
- no WHIR/STIR work
- no generic benchmark artifact pipeline redesign
- no random or flaky fuzzing in the required test path

## Why This Follow-Up Exists

The current branch already completes the required shared-layer plan. This follow-up plan was created to finish two optional extensions:

1. Bench usability gap:
   expose caller-provided basis data directly on the `bench_z2k_ring_switch_*`
   CLIs so performance comparisons do not require ad hoc code edits.

2. Confidence gap:
   provided-basis correctness is no longer a one-off demo, but most positive end-to-end
   coverage still clusters around one non-polynomial `GR(4,2)` pair in
   `tests/test_z2k_ring_switch_pcs.cpp`.

These are good follow-up tasks because they improve performance comparability and test confidence without changing the main protocol implementation.

## Current State

Today the repo has:

- a completed library-level general-basis implementation for ring-switch
- fixed-width serializer/proof-size support under provided bases
- bench binaries that now support both default and caller-provided basis modes
- provided-basis tests concentrated in `GR(4,2)`

The current optional-extension surface mainly lives in:

- `bench/bench_z2k_ring_switch_common.hpp`
- `bench/bench_z2k_ring_switch_commit.cpp`
- `bench/bench_z2k_ring_switch_eval.cpp`
- `bench/bench_z2k_ring_switch_outer_commit.cpp`
- `bench/bench_z2k_ring_switch_outer_prove.cpp`
- `bench/bench_z2k_ring_switch_outer_verify.cpp`
- `tests/test_z2k_ring_switch_pcs.cpp`
- `README.md`

## Design Decisions To Freeze Up Front

### 1. Bench default remains polynomial mode

All existing bench invocations must remain valid with no new flags.

- default behavior stays `use_provided_basis=false`
- provided-basis mode is opt-in
- no existing release or smoke command should need to change

### 2. Provided-basis CLI uses inline coefficient syntax

For this pass, the clean and bounded interface is:

- `--basis-mode default|provided`
- `--alpha-basis <elem0;elem1;...>`
- `--beta-basis <elem0;elem1;...>`
- optional:
  - `--alpha-dual-basis <elem0;elem1;...>`
  - `--beta-dual-basis <elem0;elem1;...>`

Element encoding:

- each basis element is encoded in the active polynomial basis of `GR(2^k, 2^kappa)`
- element coefficients use the same low-to-high comma-separated convention already used by `--ring-F`
- basis elements are separated by semicolons

Example:

```bash
./build/bench_z2k_ring_switch_eval \
  --basis-mode provided \
  --alpha-basis '1,1;0,1' \
  --beta-basis '1,0;1,1' \
  --queries 1 --warmup 0 --reps 3
```

Rationale:

- this keeps the first CLI surface self-contained and scriptable
- it matches the repo's existing coefficient-list style
- it avoids inventing a new file format

Explicit basis-file loading can be a later extension if large-dimension use becomes common.

### 3. Dual bases stay optional at the CLI boundary

If the user passes `--basis-mode provided`:

- `--alpha-basis` and `--beta-basis` are required
- `--alpha-dual-basis` and `--beta-dual-basis` are optional
- if a dual basis is omitted, setup should derive it as today

This keeps the CLI practical while still allowing exact reproductions when the caller wants to pin a specific dual basis.

### 4. High-confidence tests must be deterministic

Required test coverage should come from:

- deterministic context specs
- deterministic basis-generation families
- deterministic witness tables and query points

The required path should not rely on probabilistic fuzzing or large random sweeps.

## Deliverables

By the end of this follow-up plan, the repo should have:

- all five ring-switch bench binaries supporting caller-provided basis mode
- one shared bench helper layer that owns basis CLI parsing and setup construction
- README documentation for the new bench usage
- a deterministic, multi-context provided-basis test matrix with both positive and negative cases
- explicit stop-rule validation commands for both default mode and provided-basis mode

## Recommended Work Breakdown

### WP0. Freeze CLI contract and coverage matrix

Files:

- `ring_switch_general_basis_optional_extensions_plan.md`
- `README.md`

Tasks:

- freeze the provided-basis bench CLI syntax listed above
- freeze the required context matrix for high-confidence tests
- freeze the reporting contract:
  - benches should print `basis mode default|provided`
  - in provided mode, print basis dimensions and whether dual bases were supplied or derived

Required context matrix:

- Context A: `GR(4,2)` with the current default backend
- Context B: `GR(4,2)` with the current backend variant already used in tests
- Context C: `GR(4,4)` reusing the existing feasible degree-4 extension pattern already present in Frobenius tests

Required basis-family matrix:

- polynomial baseline
- provided non-polynomial `alpha` only
- provided non-polynomial `beta` only
- provided non-polynomial `alpha != beta`
- provided basis with explicit duals
- provided basis with omitted duals so setup derives them

Required negative matrix:

- singular `alpha`
- singular `beta`
- wrong-size `alpha`
- wrong-size `beta`
- malformed `alpha` dual
- malformed `beta` dual
- missing `alpha` or missing `beta` in provided mode

Acceptance criteria:

- this plan is explicit enough that implementation can proceed without re-opening interface questions

Estimated effort:

- small

### WP1. Unify ring-switch bench setup plumbing in the common helper

Files:

- `bench/bench_z2k_ring_switch_common.hpp`
- `bench/bench_z2k_ring_switch_commit.cpp`
- `bench/bench_z2k_ring_switch_eval.cpp`
- `bench/bench_z2k_ring_switch_outer_commit.cpp`
- `bench/bench_z2k_ring_switch_outer_prove.cpp`
- `bench/bench_z2k_ring_switch_outer_verify.cpp`

Tasks:

- move the remaining duplicated ring-switch bench parsing/build helpers from
  `bench_z2k_ring_switch_commit.cpp` and `bench_z2k_ring_switch_eval.cpp`
  into `bench_z2k_ring_switch_common.hpp`
- define a shared CLI-facing basis config, for example:

```cpp
enum class BasisMode { kDefault, kProvided };

struct BasisCliData {
  bool has_basis = false;
  std::vector<NTL::ZZ_pE> basis;
  bool has_dual_basis = false;
  std::vector<NTL::ZZ_pE> dual_basis;
};

struct RingSwitchBenchCliConfig {
  ContextSpec context;
  BasisMode basis_mode = BasisMode::kDefault;
  BasisCliData alpha;
  BasisCliData beta;
};
```

- add shared parsers:
  - parse semicolon-separated basis lists
  - parse optional dual-basis lists
  - validate basis element counts against `deg(F)`
- add one shared setup builder that maps the CLI config into `RingSwitchPCSSetupInput`
- add one shared help-text fragment for basis flags so five binaries do not diverge

Acceptance criteria:

- all basis parsing and setup construction live in the common helper
- commit/eval benches no longer maintain a separate duplicated setup path

Estimated effort:

- medium

### WP2. Wire provided-basis mode into all five bench binaries

Files:

- `bench/bench_z2k_ring_switch_commit.cpp`
- `bench/bench_z2k_ring_switch_eval.cpp`
- `bench/bench_z2k_ring_switch_outer_commit.cpp`
- `bench/bench_z2k_ring_switch_outer_prove.cpp`
- `bench/bench_z2k_ring_switch_outer_verify.cpp`
- `README.md`

Tasks:

- add the new flags to each binary
- keep default polynomial mode unchanged
- in provided mode:
  - require both `alpha` and `beta`
  - accept optional duals
  - route all setup through the shared helper
- update printed output to include basis mode and provided-basis metadata
- document at least one provided-basis invocation per bench family in `README.md`

Acceptance criteria:

- each bench binary still runs with old default commands
- each bench binary also runs with `--basis-mode provided`
- help text across all five binaries uses the same basis-flag wording

Estimated effort:

- medium

### WP3. Add reusable ring-switch test context and basis-generation helpers

Files:

- `tests/test_z2k_ring_switch_pcs.cpp`
- optionally `tests/test_common.hpp` if a tiny shared helper is clearly useful

Tasks:

- factor the current ad hoc `GR(4,2)` provided-basis helpers into reusable test helpers
- introduce a context-spec layer similar in spirit to the existing Frobenius feasibility tests
- add deterministic builders for:
  - extension polynomials
  - backend params
  - witness tables
  - query points
- add deterministic basis synthesis from base-ring change-of-basis matrices

Recommended basis synthesis strategy:

- start from the active polynomial basis
- left-multiply it by deterministic invertible matrices over the base ring
- use at least two families:
  - unitriangular transforms for sparse, easy-to-read cases
  - dense invertible transforms for less structured cases

Recommended helper surface:

- `BuildPolynomialBasisCase(...)`
- `BuildTransformedBasisCase(..., seed, family)`
- `BuildSingularBasisCase(...)`
- `BuildBrokenDualBasisCase(...)`
- `BuildProvidedRingSwitchParamsFromSpec(...)`

Acceptance criteria:

- adding a new context or basis family no longer requires hand-copying one-off helpers

Estimated effort:

- medium

### WP4. Land the positive high-confidence provided-basis matrix

Files:

- `tests/test_z2k_ring_switch_pcs.cpp`

Tasks:

- keep the existing `GR(4,2)` witness pair as a named baseline
- add end-to-end `Commit/ProveEval/VerifyEval` positive tests for:
  - Context A with at least:
    - non-polynomial `alpha` only
    - non-polynomial `beta` only
    - non-polynomial `alpha != beta`
    - explicit-dual variant
    - auto-derived-dual variant
  - Context B with at least:
    - one `alpha != beta` end-to-end case
  - Context C with at least:
    - non-polynomial `alpha` only
    - non-polynomial `beta` only
    - non-polynomial `alpha != beta`
- for at least one provided-basis pair in Context A and one in Context C, also lock:
  - packing semantics against manual `beta`-basis composition
  - `A_{u||w}` recovery against manual `alpha`-basis decomposition
  - partial-evaluation recovery against direct evaluation
  - proof-size helper equals serialized byte length

Implementation rule:

- keep the total positive matrix bounded and named
- do not turn this into a giant nested sweep

Acceptance criteria:

- provided-basis correctness is no longer evidenced by only one `GR(4,2)` pair
- at least one larger-degree context (`GR(4,4)`) passes end-to-end

Estimated effort:

- medium

### WP5. Land the negative and boundary high-confidence matrix

Files:

- `tests/test_z2k_ring_switch_pcs.cpp`

Tasks:

- ensure setup rejection exists for every malformed-basis class in the matrix
- add at least one malformed-dual check in Context C, not only in `GR(4,2)`
- add one boundary case where provided mode uses explicit duals that are correct
- preserve the existing `ell_prime=0` and tamper-rejection regressions
- make sure serializer/proof-size regressions still run under a provided-basis case outside the original witness pair

Acceptance criteria:

- setup validation is covered by more than one context shape
- the test suite would catch both semantic mix-ups and basis-validation regressions

Estimated effort:

- small to medium

### WP6. Validation, docs, and stop rule

Files:

- `README.md`
- optionally `README_zh.md`

Tasks:

- update the ring-switch section with:
  - bench provided-basis syntax
  - a short example command
  - a note that inline basis encoding is in polynomial-basis coordinates
- run the required validation matrix
- record any deliberately deferred items explicitly

Required validation:

```bash
cmake -S . -B build
cmake --build build --target \
  test_z2k_ring_switch_pcs \
  bench_z2k_ring_switch_commit \
  bench_z2k_ring_switch_eval \
  bench_z2k_ring_switch_outer_commit \
  bench_z2k_ring_switch_outer_prove \
  bench_z2k_ring_switch_outer_verify -j

ctest --test-dir build --output-on-failure -R test_z2k_ring_switch_pcs

./build/bench_z2k_ring_switch_commit --warmup 0 --reps 1
./build/bench_z2k_ring_switch_eval --warmup 0 --reps 1 --queries 1
./build/bench_z2k_ring_switch_outer_commit --warmup 0 --reps 1
./build/bench_z2k_ring_switch_outer_prove --warmup 0 --reps 1 --queries 1
./build/bench_z2k_ring_switch_outer_verify --warmup 0 --reps 1 --queries 1

./build/bench_z2k_ring_switch_eval \
  --basis-mode provided \
  --alpha-basis '1,1;0,1' \
  --beta-basis '1,0;1,1' \
  --warmup 0 --reps 1 --queries 1
```

Recommended extra validation:

```bash
ctest --test-dir build --output-on-failure -R 'test_z2k_ring_switch_pcs|test_z2k_frobenius_pcs'
./build/bench_z2k_ring_switch_outer_prove \
  --basis-mode provided \
  --alpha-basis '1,1;0,1' \
  --beta-basis '1,0;1,1' \
  --warmup 0 --reps 1 --queries 1
./build/bench_z2k_ring_switch_outer_verify \
  --basis-mode provided \
  --alpha-basis '1,1;0,1' \
  --beta-basis '1,0;1,1' \
  --warmup 0 --reps 1 --queries 1
```

Stop rule:

- all five bench binaries compile and run in both default and provided modes
- `test_z2k_ring_switch_pcs` stays green
- provided-basis coverage includes both `GR(4,2)` and `GR(4,4)`
- README documents the new CLI surface

Estimated effort:

- small

## File-Level Change Map

### Files expected to change

- `bench/bench_z2k_ring_switch_common.hpp`
- `bench/bench_z2k_ring_switch_commit.cpp`
- `bench/bench_z2k_ring_switch_eval.cpp`
- `bench/bench_z2k_ring_switch_outer_commit.cpp`
- `bench/bench_z2k_ring_switch_outer_prove.cpp`
- `bench/bench_z2k_ring_switch_outer_verify.cpp`
- `tests/test_z2k_ring_switch_pcs.cpp`
- `README.md`
- optionally `README_zh.md`

### Files that should not need semantic changes

- `src/Compiler/Z2k/RingSwitchPCS.cpp`
- `src/Compiler/Z2k/RingSwitchProofSerialize.cpp`
- `include/Compiler/Z2k/RingSwitchPCS.hpp`
- `include/GaloisRing/Basis.hpp`
- `src/GaloisRing/Basis.cpp`

If implementation pressure starts pushing protocol semantics back into this plan, stop and split that work into a separate protocol change.

## Main Risks And Mitigations

### Risk 1. Bench CLI becomes awkward or inconsistent across binaries

Mitigation:

- centralize all basis parsing and help-text generation in `bench_z2k_ring_switch_common.hpp`
- do not let five binaries each invent their own flag semantics

### Risk 2. Commit/eval benches keep their duplicated parser logic

Mitigation:

- explicitly refactor commit/eval onto the common helper first
- treat that refactor as part of WP1, not as optional cleanup

### Risk 3. High-confidence coverage turns into a huge slow sweep

Mitigation:

- use named deterministic cases rather than broad random sweeps
- require a bounded matrix with explicit contexts and basis families

### Risk 4. New tests accidentally only prove one more hard-coded witness

Mitigation:

- derive multiple bases from deterministic invertible matrices
- cover both degree-2 and degree-4 contexts
- require at least one context beyond the original `GR(4,2)` pair

### Risk 5. Basis syntax is too opaque to use in practice

Mitigation:

- print a concise encoding note in `--help`
- add one README example
- report basis mode and basis dimensions in benchmark output

## Total Effort Estimate

If executed cleanly as written:

- bench CLI follow-up: about 2 to 3 days
- high-confidence test follow-up: about 2 to 4 days
- total: about 4 to 7 days

This is a medium follow-up, not a protocol rewrite.

## Recommended Execution Order

1. WP1 bench common-helper unification
2. WP2 bench CLI wiring and docs
3. WP3 reusable ring-switch test helper layer
4. WP4 positive high-confidence matrix
5. WP5 negative/boundary matrix
6. WP6 validation sweep and doc polish

Doing benches first keeps the CLI contract fixed before the test matrix starts growing.
