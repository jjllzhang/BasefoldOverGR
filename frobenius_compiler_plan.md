# Frobenius Compiler Implementation Plan

Last updated: 2026-03-13  
Working branch: `frobenius-compiler-feasibility`

## Goal

Implement the paper's Frobenius-maps-based compiler (`Protocol 2`) for multilinear PCS over `Z_{2^k}` on top of the existing PCS-over-`GR(2^k, r)` backend boundary already used by the ring-switch compiler.

The first production target is:

- correctness first,
- single-point evaluation proofs only,
- composed with the existing `Z2kPCSBackend` / BaseFold backend,
- no change to the backend API,
- no attempt to implement paper-wide extras such as multipoint openings, proof composition, or WHIR.

## Current Ground Truth

### Already true in the repo

- The generic backend layer is already split out and reusable:
  - `include/Compiler/Z2k/PCSBackend.hpp`
  - `src/Compiler/Z2k/PCSBackend.cpp`
  - `include/Compiler/Z2k/BaseFoldBackendAdapter.hpp`
  - `src/Compiler/Z2k/BaseFoldBackendAdapter.cpp`
- The current compiler line in production is ring-switch only:
  - `include/Compiler/Z2k/RingSwitchPCS.hpp`
  - `src/Compiler/Z2k/RingSwitchPCS.cpp`
- A correctness-oriented feasibility prototype now exists:
  - `tests/test_z2k_frobenius_feasibility.cpp`

### What the feasibility prototype proved

On small parameters `GR(4,2)` and `GR(4,4)`, the current NTL/Galois-ring stack can support the core algebra needed by `Protocol 2`:

- a usable Frobenius action `tau` / `sigma`,
- a normal basis `beta`,
- a dual basis `alpha`,
- the identity `Tr(alpha_u * beta_v) = delta_uv`,
- the paper's partial-recovery formula
  - `t(u || r_suffix) = sum_i tau^i(alpha_u * s_i)`,
- recombination of recovered partial evaluations back to the original claim.

This is enough to say the compiler is implementable in this codebase.

### What is still missing

The feasibility prototype is intentionally not production code:

- it is test-only,
- it uses brute-force enumeration for basis-coordinate recovery,
- it does not expose reusable algebra APIs,
- it does not implement any compiler proof objects, serializers, or benches,
- it does not preserve production performance constraints.

The next work is therefore not "port the test into the compiler"; it is "turn the validated algebra into a library-quality layer, then build the compiler on top of it."

## Non-Goals

The initial Frobenius compiler line should not try to solve all adjacent problems at once.

Out of scope for the first landing:

- multipoint opening (`Protocol 7`),
- Blaze/Orion style proof composition,
- WHIR / 3-foldable FRI,
- general basis-polymorphism across every existing compiler surface,
- replacing the ring-switch compiler,
- aggressive performance tuning before correctness is locked.

## High-Level Design

### Design rule 1: keep the backend boundary unchanged

The Frobenius compiler should keep using:

- `Z2kPCSBackendHandle`
- `Z2kPCSBackendCommitArtifacts`
- `Z2kPCSBackendEvalProof`
- `Z2kPCSBackendCommit(...)`
- `Z2kPCSBackendProveEval(...)`
- `Z2kPCSBackendVerifyEval(...)`

The backend still only needs to prove one evaluation:

- committed polynomial over `GR(2^k, r)`,
- point `r'` in `GR(2^k, r)^(ell - kappa)`,
- value `t_star = t'(r')`.

This is the same backend contract already used by ring-switch.

### Design rule 2: do not destabilize ring-switch while Frobenius is landing

Do not refactor `RingSwitchPCS` into a generic "basis-aware compiler" in the first pass.

Instead:

- add a separate Frobenius compiler surface,
- add a separate low-level algebra layer for Frobenius/normal-basis operations,
- only extract shared helpers later if duplication is clearly worth it.

This keeps the current ring-switch path stable and testable while Frobenius is being built.

### Design rule 3: no brute force in production paths

The feasibility test used enumeration because it was the fastest way to answer "is this mathematically viable here?"

That approach must not appear in production code.

Production coordinate recovery should use trace and explicit bases:

- if `x = sum_i c_i * beta_i`, then `c_i = Tr(alpha_i * x)` in `Z_{2^k}`,
- `tau(x)` can therefore be implemented as:
  1. recover normal-basis coordinates via `Tr(alpha_i * x)`,
  2. rotate the coordinates,
  3. recombine with `beta_i`.

This gives a clear, correct first implementation even before optimization.

## Target File Layout

### New algebra layer

Add a new low-level module under `GaloisRing/`:

- `include/GaloisRing/FrobeniusBasis.hpp`
- `src/GaloisRing/FrobeniusBasis.cpp`

Expected responsibilities:

- basis search and validation,
- dual basis construction,
- trace-based coordinate extraction,
- recombination from normal-basis coordinates,
- Frobenius maps `tau` and `sigma`,
- optional setup-time helper tables.

Suggested public data structures:

- `FrobeniusBasisParams`
- `FrobeniusBasisData`
- `NormalBasisData`

Suggested public functions:

- `FindNormalBasisOrThrow(...)`
- `BuildDualBasisOrThrow(...)`
- `RecoverNormalBasisCoords(...)`
- `ComposeFromNormalBasisCoords(...)`
- `ApplyFrobeniusTau(...)`
- `ApplyFrobeniusSigma(...)`
- `TraceToBaseRing(...)`

### New compiler layer

Add a new compiler surface under `Compiler/Z2k/`:

- `include/Compiler/Z2k/FrobeniusPCS.hpp`
- `src/Compiler/Z2k/FrobeniusPCS.cpp`
- `include/Compiler/Z2k/FrobeniusProofSerialize.hpp`
- `src/Compiler/Z2k/FrobeniusProofSerialize.cpp`

Suggested public structs:

- `FrobeniusPCSSetupInput`
- `FrobeniusPCSParams`
- `FrobeniusPCSCommitArtifacts`
- `FrobeniusPCSOuterCommitArtifacts`
- `FrobeniusPCSOuterEvalProof`
- `FrobeniusPCSEvalProof`
- `FrobeniusPCSProver`
- `FrobeniusPCSVerifier`

Proof shape should mirror the ring-switch split where practical:

- outer-only proof:
  - `s_by_i`
  - `h_by_level`
  - `t_star`
- composed proof:
  - outer proof fields
  - backend proof

### Tests

Keep the feasibility prototype and grow around it rather than deleting it.

Expected test layout:

- keep:
  - `tests/test_z2k_frobenius_feasibility.cpp`
- add:
  - `tests/test_z2k_frobenius_pcs.cpp`

The feasibility test remains the algebraic stop-rule regression.
The new PCS test owns setup/pack/prove/verify/serialize coverage.

### Bench and docs

After correctness is stable:

- add `bench/bench_z2k_frobenius_commit.cpp`
- add `bench/bench_z2k_frobenius_eval.cpp`
- update `README.md`
- update `README_zh.md`

## Implementation Strategy

## Phase 0: Feasibility and Planning

Status: `[x] Done`

Completed work:

- isolated branch for Frobenius work,
- created a standalone feasibility regression,
- validated Protocol 2 algebra on `GR(4,2)` and `GR(4,4)`.

Artifacts:

- `tests/test_z2k_frobenius_feasibility.cpp`
- this plan document

## Phase 1: Production Algebra Layer

Status: `[x] Done`

### Objective

Replace test-only brute-force logic with reusable, deterministic, production-quality algebra helpers.

### Tasks

- `[x]` Create `include/GaloisRing/FrobeniusBasis.hpp`
- `[x]` Create `src/GaloisRing/FrobeniusBasis.cpp`
- `[x]` Implement setup-time normal basis search
- `[x]` Implement setup-time dual basis construction
- `[x]` Implement trace-based coordinate extraction
- `[x]` Implement recombination from normal-basis coordinates
- `[x]` Implement `tau` and `sigma`
- `[x]` Add unit/shape validation helpers

Completed artifacts:

- `include/GaloisRing/FrobeniusBasis.hpp`
- `src/GaloisRing/FrobeniusBasis.cpp`
- `tests/test_galois_ring_basic.cpp`
- `tests/test_z2k_frobenius_feasibility.cpp`

### Recommended concrete approach

#### 1. Normal basis search

Start with setup-time search only. That is acceptable because setup is offline and much less latency-sensitive than proving.

Suggested search order:

1. powers of the Teichmuller generator,
2. small affine combinations of Teichmuller powers,
3. fail loudly if no normal basis is found.

Success criterion for a candidate `beta_0`:

- the orbit `beta_i = tau^i(beta_0)` has length `r`,
- the orbit spans the extension as a basis over `Z_{2^k}`.

#### 2. Dual basis construction

Do not use brute-force search in production.

Use a setup-time linear solve over `Z_{2^k}`:

- choose polynomial basis `omega_j = x^j`,
- build the trace matrix
  - `M[u,j] = Tr(beta_u * omega_j)`,
- solve `M * a = e_u` over `Z_{2^k}` for each target coordinate `u`,
- set `alpha_u = sum_j a_j * omega_j`.

This requires a small solver over base-ring scalars with unit pivots.

#### 3. Coordinate extraction

Recover coefficients in the normal basis by

- `c_i = Tr(alpha_i * x)`.

This should return base-ring scalars directly and avoids generic matrix inversion per query.

#### 4. Frobenius map

Implement

- `tau(x)` by extract -> rotate -> recombine,
- `sigma(x)` as `tau^(r-1)(x)`.

This is correctness-first and directly backed by the feasibility test.

### Acceptance criteria

- `tests/test_z2k_frobenius_feasibility.cpp` no longer needs private brute-force logic copied from production,
- `tau` fixes base-ring constants,
- `tau^r = id`,
- `Tr(alpha_u * beta_v) = delta_uv`,
- coordinate extraction round-trips.

### Stop rule

Do not proceed to Phase 2 if production coordinate extraction still depends on enumerating all ring elements or all coefficient vectors.

## Phase 2: Frobenius Setup and Packing

Status: `[x] Done`

### Objective

Create a compiler setup path and a packing path analogous to ring-switch, but specialized to the normal-basis Frobenius route.

### Tasks

- `[x]` Add `FrobeniusPCSSetupInput`
- `[x]` Add `FrobeniusPCSParams`
- `[x]` Validate:
  - active `ZZ_p` modulus,
  - active `ZZ_pE` modulus,
  - `r = 2^kappa`,
  - backend message length is `2^(ell-kappa)`,
  - normal basis and dual basis dimensions match `r`
- `[x]` Implement packing from `t_table` over `Z_{2^k}` to `t'_table` over `GR(2^k, r)` using `beta`
- `[x]` Add commit artifact structs and builder APIs

Completed artifacts:

- `include/Compiler/Z2k/FrobeniusPCS.hpp`
- `src/Compiler/Z2k/FrobeniusPCS.cpp`
- `tests/test_z2k_frobenius_pcs.cpp`

### Recommended concrete approach

Mirror the current ring-switch user-facing contract:

- compiler still accepts a Boolean-hypercube table `t_table`,
- compiler packs it into a `GR` table `t'_table`,
- packed table is converted to monomial coefficients only where the backend needs them.

Packing formula:

- `t'(w) = sum_v beta_v * t(v || w)`.

### Acceptance criteria

- packing round-trips via the dual basis:
  - `Tr(alpha_u * t'(w)) = t(u || w)`,
- setup rejects mismatched dimensions and contexts,
- commit path produces the same commitment as the backend run on the packed polynomial.

## Phase 3: Outer Proof, Correctness First

Status: `[x] Done`

### Objective

Implement the outer proof and verifier for the Frobenius compiler before optimizing anything.

### Tasks

- `[x]` Implement transcript config and domain separator
- `[x]` Implement outer prove path
- `[x]` Implement outer verify path
- `[x]` Implement composed prove/verify path
- `[x]` Add prover/verifier wrappers similar to ring-switch

Completed artifacts:

- `include/Compiler/Z2k/FrobeniusPCS.hpp`
- `src/Compiler/Z2k/FrobeniusPCS.cpp`
- `tests/test_z2k_frobenius_pcs.cpp`

### Protocol mapping

#### Step 1: send orbit evaluations

Compute:

- `s_i = t'(sigma^i(r_suffix))`.

Store these as `proof.s_by_i`.

#### Step 2: Equality Check 1

Recover:

- `partials[u] = sum_i tau^i(alpha_u * s_i)`.

Then recombine against `eq(z_prefix; u)` and compare to the claimed value.

#### Step 3: sample `r'_prefix`

Use the transcript after absorbing:

- commitment,
- query point `z`,
- claimed value,
- all `s_i`.

#### Step 4: batched sumcheck

Define:

- `lambda_i = eq(r'_prefix; v(i))`
- `h(X) = t'(X) * sum_i lambda_i * eq(X; sigma^i(r_suffix))`.

Correctness-first implementation:

- explicitly materialize `sigma^i(r_suffix)` for all `i`,
- build the batched `g_table[w] = sum_i lambda_i * eq(sigma^i(r_suffix); w)`,
- run the existing `ProductSumcheckProver`.

This is the simplest production-correct version.
Do not optimize this until the proof is stable and tested.

#### Step 5: final check + backend proof

Compute:

- `t_star = t'(r'_suffix)`
- `g_star = sum_i lambda_i * eq(r'_suffix; sigma^i(r_suffix))`

Check:

- final sumcheck claim equals `t_star * g_star`,
- then invoke the backend proof on `r'`.

### Acceptance criteria

- honest proof verifies,
- direct prove and artifact-based prove agree on outer messages,
- tampering in `s_by_i`, `h_by_level`, or `t_star` is rejected,
- outer-only verify passes independently of backend composition.

## Phase 4: Serialization, Proof Size, and Benches

Status: `[x] Done`

### Objective

Bring the Frobenius compiler up to the same repo standard as the ring-switch compiler.

### Tasks

- `[x]` Add fixed-width serializer
- `[x]` Add proof-size accounting through the serializer path
- `[x]` Add outer-only and composed proof-size helpers
- `[x]` Add commit/eval benches
- `[x]` Document proof-size semantics and bench semantics

Completed artifacts:

- `include/Compiler/Z2k/FrobeniusProofSerialize.hpp`
- `src/Compiler/Z2k/FrobeniusProofSerialize.cpp`
- `bench/bench_z2k_frobenius_common.hpp`
- `bench/bench_z2k_frobenius_commit.cpp`
- `bench/bench_z2k_frobenius_eval.cpp`
- `tests/test_z2k_frobenius_pcs.cpp`
- `README.md`
- `README_zh.md`

### Design rule

Do not hand-maintain a separate size formula once serialization exists.

Follow the same contract as ring-switch:

- outer size counts only outer proof fields,
- composed size appends backend proof bytes with an explicit length prefix,
- public inputs and commitment are excluded from proof-size reporting.

### Acceptance criteria

- serialized outer bytes equal outer size count,
- serialized composed bytes equal composed size count,
- bench output is wired to serializer-backed byte counts,
- README documents the scope clearly.

## Phase 5: Optimization and Refactor Cleanup

Status: `[ ] Not started`

### Objective

Reduce obvious overhead only after correctness, serialization, and docs are in place.

### Candidate optimizations

- `[ ]` Cache normal-basis coordinates of the queried suffix point
- `[ ]` Cache orbit points `sigma^i(r_suffix)`
- `[ ]` Avoid repeated recomputation of `eq(sigma^i(r_suffix); w)`
- `[ ]` Precompute setup-time tables for `tau`-rotation and basis recombination
- `[ ]` Factor common proof-serialize helpers if Frobenius and ring-switch drift too close

### Explicit non-goal for this phase

Do not collapse Frobenius and ring-switch into one generic compiler unless profiling shows the duplication is actually a maintenance problem.

## Testing Plan

## Algebraic regressions

- keep and extend:
  - `tests/test_z2k_frobenius_feasibility.cpp`

Required coverage:

- `GR(4,2)` and `GR(4,4)` still pass,
- `tau` and `trace` regressions,
- normal basis / dual basis regressions,
- coordinate round-trip regressions,
- partial recovery regressions.

## PCS regressions

Add `tests/test_z2k_frobenius_pcs.cpp` with:

- setup success and rejection cases,
- packing round-trip,
- commit agreement with backend,
- honest outer proof verify,
- honest composed proof verify,
- tamper rejection,
- dimension-zero or smallest supported dimension cases,
- serialize/size consistency.

## Validation commands

Phase 1:

```bash
cmake -S . -B build
cmake --build build -j 4 --target test_galois_ring test_z2k_frobenius_feasibility
ctest --test-dir build --output-on-failure -R "test_galois_ring|test_z2k_frobenius_feasibility"
```

Phase 2-4:

```bash
cmake --build build -j 4 --target test_z2k_frobenius_pcs
ctest --test-dir build --output-on-failure -R "test_z2k_frobenius_feasibility|test_z2k_frobenius_pcs"
```

Phase 4+:

```bash
cmake --build build -j 4 --target bench_z2k_frobenius_commit bench_z2k_frobenius_eval
./build/bench_z2k_frobenius_commit --warmup 0 --reps 1
./build/bench_z2k_frobenius_eval --warmup 0 --reps 1 --queries 2
```

## Risks and Mitigations

### Risk 1: normal basis search is fragile for some extension polynomials

Mitigation:

- keep the feasibility regression,
- make setup fail loudly with a clear error,
- start from the exact extension polynomials already used in tests and benches,
- search deterministically before adding randomness.

### Risk 2: dual basis construction needs a solver over `Z_{2^k}`

Mitigation:

- isolate a small base-ring linear solver in the algebra layer,
- require unit pivots only,
- test it independently before using it in compiler setup.

### Risk 3: `tau` via extract -> rotate -> recombine is too slow

Mitigation:

- accept this for the first correct version,
- profile before optimizing,
- add caches only after proof correctness is locked.

### Risk 4: premature generalization harms ring-switch stability

Mitigation:

- keep Frobenius codepaths separate,
- do not widen `RingSwitchBasisKind` in the first implementation pass,
- only extract shared utilities after both paths are green.

## Definition of Done

The Frobenius compiler is "done" for the first landing when all of the following are true:

- a production algebra layer exists under `GaloisRing/`,
- no production path relies on brute-force enumeration,
- setup/commit/prove/verify work end to end,
- proof serialization and proof-size accounting are live,
- honest/tampered tests pass,
- ring-switch tests still pass,
- README documents the scope and limitations,
- the feasibility regression remains green.

## Immediate Next Step

The next implementation step should be:

1. create `include/GaloisRing/FrobeniusBasis.hpp` and `src/GaloisRing/FrobeniusBasis.cpp`,
2. move the algebra that is now only in `tests/test_z2k_frobenius_feasibility.cpp` into that layer,
3. replace brute-force coordinate recovery with trace-based recovery plus setup-time dual-basis construction,
4. only then start `FrobeniusPCSSetup(...)`.
