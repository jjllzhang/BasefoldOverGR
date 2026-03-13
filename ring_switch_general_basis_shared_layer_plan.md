# Ring-Switch General Basis Shared-Layer Plan

Status: WP0-WP1 completed; WP2-WP6 pending
Last updated: 2026-03-14
Scope: make `RingSwitchPCS` support paper-aligned general `alpha/beta` Galois-ring bases through a clean shared basis layer

## Progress

- 2026-03-14: `WP0` landed on branch `feat/ring-switch-general-basis`.
- Public `RingSwitch` setup/params types now carry explicit `alpha/beta` basis
  data through `GaloisRingBasisData`.
- Default setup remains the active polynomial basis on both sides.
- Caller-provided `alpha/beta` basis mode is now expressible in the public API
  but intentionally rejected by `RingSwitchPCSSetup(...)` until `WP2`, so the
  repo does not pretend to support semantics that have not been implemented yet.
- 2026-03-14: `WP1` landed in the working tree.
- Generic basis algebra now lives in `include/GaloisRing/Basis.hpp` and
  `src/GaloisRing/Basis.cpp`.
- `FrobeniusBasis.*`, `FrobeniusPCS.cpp`, and `RingSwitchPCS.cpp` now consume
  the shared basis helpers instead of owning generic basis algebra privately.

## Goal

Implement general-basis support for the ring-switch compiler (`Protocol 1` / Appendix C.1 semantics), so that:

- `beta` is the packing basis for `t'(w) = sum_v t(v||w) * beta_v`
- `alpha` is the decomposition basis for `eq(r_suffix; w) = sum_u A_{u||w} * alpha_u`
- the prover/verifier can recover partial evaluations when `alpha` and `beta` are arbitrary valid `Z_{2^k}`-bases of `GR(2^k, 2^kappa)`
- the default public behavior remains the current polynomial-basis specialization

The plan intentionally keeps the backend PCS boundary unchanged and does not change the proof message shape.

## Non-Goals

- No multipoint opening work
- No Blaze-Orion composition work
- No WHIR/STIR work
- No BaseFold backend protocol changes
- No ring-switch proof format redesign

## Current Gaps

Current code is still specialized to the active polynomial basis:

- public API now carries basis data, but setup still only enables the default
  polynomial-basis path
- packing assumes the polynomial block layout
- `eq(r_suffix; w)` is decomposed using polynomial coefficients
- `s_u` is decomposed using polynomial coefficients
- partial-evaluation recovery implicitly relies on `alpha = beta = polynomial basis`

These assumptions currently live mainly in:

- `include/Compiler/Z2k/RingSwitchPCS.hpp`
- `src/Compiler/Z2k/RingSwitchPCS.cpp`
- `src/Compiler/Z2k/RingSwitchProofSerialize.cpp`
- `tests/test_z2k_ring_switch_pcs.cpp`
- `bench/bench_z2k_ring_switch_common.hpp`

## Recommended Design

### 1. Introduce a generic shared basis algebra layer

Add a new shared module under `GaloisRing/` that is not Frobenius-specific.

Recommended new files:

- `include/GaloisRing/Basis.hpp`
- `src/GaloisRing/Basis.cpp`

Recommended public types:

```cpp
struct GaloisRingBasisData {
  std::vector<NTL::ZZ_pE> basis;
  std::vector<NTL::ZZ_pE> dual_basis;
};
```

Recommended public helpers:

- `std::vector<NTL::ZZ_pE> BuildPolynomialBasis(long dimension);`
- `std::vector<NTL::ZZ_pE> BuildDualBasisOrThrow(const std::vector<NTL::ZZ_pE>& basis);`
- `void ValidateBasisShapeOrThrow(const std::vector<NTL::ZZ_pE>& basis, const char* label, const char* func_name);`
- `void ValidateBasisDataOrThrow(const GaloisRingBasisData& basis, const char* label, const char* func_name);`
- `std::vector<NTL::ZZ_p> RecoverBasisCoordsOrThrow(const GaloisRingBasisData& basis, const NTL::ZZ_pE& element, const char* func_name);`
- `NTL::ZZ_pE ComposeFromBasisCoordsOrThrow(const std::vector<NTL::ZZ_pE>& basis, const std::vector<NTL::ZZ_p>& coords, const char* func_name);`
- `NTL::ZZ_p TraceToBaseRing(const NTL::ZZ_pE& element);`
- `bool IsBaseRingConstant(const NTL::ZZ_pE& value);`

Implementation note:

- move the generic linear-algebra pieces currently buried in `src/GaloisRing/FrobeniusBasis.cpp` into this new shared file
- keep Frobenius-only logic such as normal-basis discovery, Teichmuller validation, and tau/sigma helpers inside `FrobeniusBasis.*`

### 2. Replace ring-switch basis descriptors with canonicalized basis data

Current ring-switch public setup API only carries `kind + dimension`, which is not enough for general paper-style bases.

Recommended ring-switch API shape:

```cpp
struct RingSwitchPCSProvidedBasisInput {
  bool has_alpha_basis = false;
  GaloisRingBasisData alpha_basis;
  bool has_beta_basis = false;
  GaloisRingBasisData beta_basis;
};

struct RingSwitchPCSSetupInput {
  long ell = 0;
  long kappa = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  bool use_provided_basis = false;
  RingSwitchPCSProvidedBasisInput provided_basis;
  Z2kPCSBackendHandle backend;
};

struct RingSwitchPCSParams {
  long ell = 0;
  long kappa = 0;
  long ell_prime = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  GaloisRingBasisData alpha_basis;
  GaloisRingBasisData beta_basis;
  Z2kPCSBackendHandle backend;
};
```

Recommended setup semantics:

- `use_provided_basis=false`:
  - build the active polynomial basis for both `alpha` and `beta`
  - derive their dual bases with shared algebra helpers
- `use_provided_basis=true`:
  - require both `alpha` and `beta`
  - if a provided basis omits its `dual_basis`, derive it in setup
  - validate that both are genuine `Z_{2^k}`-bases in the active ring context

Important semantic point:

- `alpha` and `beta` do not need to be equal
- `alpha` and `beta` do not need to be dual to each other
- each one should, however, be stored together with its own dual basis for coordinate recovery

### 3. Keep proof shape unchanged

The current ring-switch proof object already matches the paper-facing outer/composed message shape:

- `s_by_u`
- `h_by_level`
- `t_star`
- backend opening proof

General-basis support should not add basis data into the proof.

Reason:

- basis data is part of setup/params, not per-proof transcript traffic
- verifier already requires `RingSwitchPCSParams`

This means serializer changes should be mechanical only.

## Concrete Work Packages

### WP0. Freeze the public API shape

Status: completed on 2026-03-14

Files:

- `include/Compiler/Z2k/RingSwitchPCS.hpp`
- `README.md`

Tasks:

- delete `RingSwitchBasisKind` / `RingSwitchBasisDescriptor` from the ring-switch public API
- define the new setup/params basis-carrying structs
- decide that default setup remains polynomial-basis on both sides
- explicitly document that `alpha` and `beta` are independent bases

Acceptance criteria:

- no remaining ring-switch public type depends on `kind + dimension` only
- the header can express both the default polynomial mode and caller-provided basis mode

Estimated effort:

- small

### WP1. Extract shared basis algebra from `FrobeniusBasis`

Status: completed on 2026-03-14

Files:

- new `include/GaloisRing/Basis.hpp`
- new `src/GaloisRing/Basis.cpp`
- `include/GaloisRing/FrobeniusBasis.hpp`
- `src/GaloisRing/FrobeniusBasis.cpp`
- `CMakeLists.txt`

Tasks:

- move generic basis helpers out of `FrobeniusBasis.cpp`
- keep Frobenius-only functionality in `FrobeniusBasis.*`
- make `FrobeniusBasis.*` consume the shared helpers instead of owning them privately
- update build wiring in `galoisring_algebra`

Recommended extraction targets from current `FrobeniusBasis.cpp`:

- polynomial-basis construction
- coordinate-matrix construction
- invertibility check over the base ring
- linear-system solver with unit pivots
- dual-basis construction
- coordinate recovery from dual basis
- composition from basis coordinates
- `TraceToBaseRing`
- `IsBaseRingConstant`

Acceptance criteria:

- `test_z2k_frobenius_pcs` remains green without semantic change
- shared basis helpers are available without including `FrobeniusPCS` or normal-basis-specific APIs

Estimated effort:

- medium

### WP2. Refactor ring-switch setup and validation onto shared basis data

Files:

- `include/Compiler/Z2k/RingSwitchPCS.hpp`
- `src/Compiler/Z2k/RingSwitchPCS.cpp`
- `tests/test_z2k_ring_switch_pcs.cpp`
- `bench/bench_z2k_ring_switch_common.hpp`

Tasks:

- replace descriptor validation with basis-data validation
- add default polynomial-basis construction in setup
- add provided-basis path in setup
- derive missing dual bases during setup
- make all dimension checks come from `basis.size()`

Implementation note:

- serializer and proof-shape helpers currently read `params.beta_basis.dimension`; these will need mechanical updates after the params type changes

Acceptance criteria:

- existing polynomial-basis tests still pass through the new setup path
- setup rejects malformed basis vectors, malformed dual bases, and wrong-dimension bases

Estimated effort:

- medium

### WP3. Generalize the ring-switch core semantics

Files:

- `src/Compiler/Z2k/RingSwitchPCS.cpp`
- `tests/test_z2k_ring_switch_pcs.cpp`

Tasks:

1. Generalize packing.

- current function: `PackZ2kCoeffsToGREvals`
- current behavior: writes `t(v||w)` into polynomial coefficient slot `v`
- target behavior: compose `t'(w)` from base-ring coordinates against `beta_basis.basis`

2. Generalize decomposition for `eq(r_suffix; w)`.

- current function: `BuildRingSwitchComponentTensor`
- current behavior: decompose `eq_at_w` using polynomial coefficients
- target behavior: recover coordinates of `eq_at_w` with respect to `alpha_basis`
- output `A_{u||w}` remains a base-ring table

3. Generalize decomposition for `s_u`.

- current function: `RecoverPartialEvaluationsFromSByU`
- current behavior: decompose each `s_u` using polynomial coefficients, then recombine with polynomial-basis `alpha_u`
- target behavior:
  - decompose `s_u` using `beta_basis`
  - recombine the recovered coefficient matrix using `alpha_basis`

4. Keep the transcript/sumcheck/backend wiring unchanged.

- `ComputeSByU`, `BuildBatchedGTable`, transcript challenge flow, batched sumcheck, and backend opening all stay structurally the same once the basis-dependent tables are correct

Acceptance criteria:

- the following semantic equalities hold for non-polynomial `alpha/beta` test vectors:
  - packed `t'(w)` matches `sum_v t(v||w) * beta_v`
  - `A_{u||w}` reconstructs `eq(r_suffix; w)` in `alpha`
  - recovered partial evaluations equal the direct multilinear evaluation path
  - full `RingSwitchPCSCommit/ProveEval/VerifyEval` passes end to end

Estimated effort:

- medium to large

### WP4. Update proof-size and serializer plumbing

Files:

- `src/Compiler/Z2k/RingSwitchProofSerialize.cpp`
- possibly `include/Compiler/Z2k/RingSwitchProofSerialize.hpp`

Tasks:

- replace uses of `params.beta_basis.dimension` with `params.beta_basis.basis.size()`
- keep serialized proof bytes unchanged under default polynomial mode
- make proof-shape validation use the canonical basis size from params

Acceptance criteria:

- serializer output for current polynomial-basis proofs is unchanged or only changes where the current code depended on removed params fields
- proof-size helpers remain exact relative to the serializer

Estimated effort:

- small

### WP5. Expand tests to cover the actual design space

Files:

- `tests/test_z2k_ring_switch_pcs.cpp`
- optionally a new dedicated helper file under `tests/`

Tasks:

- keep all current polynomial-basis tests
- add a supplied non-polynomial `beta` basis test
- add a supplied non-polynomial `alpha` basis test
- add an `alpha != beta` round-trip test
- add setup rejection tests for:
  - non-basis `alpha`
  - non-basis `beta`
  - wrong-size basis
  - malformed dual basis
- add end-to-end prove/verify tests under a supplied non-polynomial pair

Recommended test structure:

- keep the current small `GR(4,2)` setting for semantic tests
- hand-construct one invertible non-polynomial basis pair over the active context
- verify direct formulas against the paper's Appendix C.1 equations before relying on end-to-end proofs

Acceptance criteria:

- `test_z2k_ring_switch_pcs` covers both default polynomial mode and provided-basis mode
- default mode remains green

Estimated effort:

- medium

### WP6. Update benches and docs

Files:

- `bench/bench_z2k_ring_switch_common.hpp`
- `bench/bench_z2k_ring_switch_commit.cpp`
- `bench/bench_z2k_ring_switch_eval.cpp`
- `bench/bench_z2k_ring_switch_outer_commit.cpp`
- `bench/bench_z2k_ring_switch_outer_prove.cpp`
- `bench/bench_z2k_ring_switch_outer_verify.cpp`
- `README.md`
- optionally `README_zh.md`

Tasks:

- keep polynomial basis as the default bench mode
- update bench helper construction to the new setup API
- optionally add a later CLI path for supplied basis data, but do not block the main implementation on CLI design
- document the two setup modes:
  - default polynomial basis
  - caller-provided `alpha/beta` basis data
- document that the proof shape is unchanged

Acceptance criteria:

- all ring-switch bench binaries compile and run in default mode
- README clearly states the new scope and setup options

Estimated effort:

- small to medium

## File-Level Change Map

### New files

- `include/GaloisRing/Basis.hpp`
- `src/GaloisRing/Basis.cpp`

### Existing files with semantic changes

- `include/Compiler/Z2k/RingSwitchPCS.hpp`
- `src/Compiler/Z2k/RingSwitchPCS.cpp`
- `tests/test_z2k_ring_switch_pcs.cpp`

### Existing files with shared-layer refactor changes

- `include/GaloisRing/FrobeniusBasis.hpp`
- `src/GaloisRing/FrobeniusBasis.cpp`
- `CMakeLists.txt`

### Existing files with mechanical follow-through

- `src/Compiler/Z2k/RingSwitchProofSerialize.cpp`
- `bench/bench_z2k_ring_switch_common.hpp`
- ring-switch bench binaries
- `README.md`
- optionally `README_zh.md`

## Validation Matrix

Run at minimum:

```bash
cmake -S . -B build
cmake --build build -j
ctest --test-dir build --output-on-failure -R 'test_z2k_ring_switch_pcs|test_z2k_frobenius_pcs'
```

Recommended extra validation:

```bash
ctest --test-dir build --output-on-failure -R 'test_pcs|test_z2k_ring_switch_pcs|test_z2k_frobenius_pcs'
./build/bench_z2k_ring_switch_outer_commit --ell 8 --kappa 2
./build/bench_z2k_ring_switch_outer_prove --ell 8 --kappa 2
./build/bench_z2k_ring_switch_outer_verify --ell 8 --kappa 2
```

Default-mode regression expectation:

- polynomial-basis proofs should still verify
- proof-size accounting should remain stable
- current bench binaries should still work without new flags

## Main Risks

### Risk 1. Mixing up `alpha` and `beta` recovery paths

Most likely semantic bug:

- decomposing `s_u` with `alpha` instead of `beta`
- reconstructing `A_{u||w}` in `beta` instead of `alpha`

Mitigation:

- add small hand-check tests for the Appendix C.1 equations before end-to-end tests

### Risk 2. Over-coupling ring-switch to Frobenius-specific APIs

Bad outcome:

- ring-switch starts depending on "normal basis" terminology or Teichmuller validation that has nothing to do with Protocol 1

Mitigation:

- keep shared basis algebra generic in a new `GaloisRing/Basis.*` layer
- make `FrobeniusBasis.*` consume that layer, not the other way around

### Risk 3. Accidental public API churn larger than needed

Bad outcome:

- touching backend APIs or proof shapes that do not need to change

Mitigation:

- keep changes isolated to setup/params/basis helpers plus the ring-switch arithmetic
- keep serializer/proof object shape stable

## Suggested Execution Order

Recommended implementation order:

1. WP1 shared basis layer extraction
2. WP2 ring-switch setup/params refactor
3. WP3 basis-dependent arithmetic refactor
4. WP5 semantic tests
5. WP4 serializer/proof-size follow-through
6. WP6 bench/docs cleanup

Reason:

- WP3 becomes much simpler once basis recovery/composition is already shared and tested
- serializer and bench changes are easier after the ring-switch params shape has stabilized

## Rough Size Estimate

Expected total churn for the clean shared-layer version:

- new code: about 250-450 lines in `GaloisRing/Basis.*`
- ring-switch semantic refactor: about 200-350 lines
- tests and bench/docs follow-through: about 250-500 lines
- total touched files: about 9-13

Overall effort:

- medium to large
- not a protocol rewrite
- the hard part is semantic discipline, not transcript or backend plumbing

## Stop Rule for "Paper-Aligned Enough"

This work should be considered complete when all of the following are true:

- ring-switch setup can run in both default polynomial mode and provided general-basis mode
- `PackZ2kCoeffsToGREvals`, `BuildRingSwitchComponentTensor`, and partial-evaluation recovery all respect independent `alpha/beta` bases
- `RingSwitchPCSCommit`, `RingSwitchPCSProveEval`, and `RingSwitchPCSVerifyEval` pass under a non-polynomial `alpha/beta` test pair
- Frobenius tests still pass after the shared-layer extraction
- README documents the new setup surface and clearly distinguishes the default specialization from the paper-general path
