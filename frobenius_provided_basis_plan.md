# Frobenius Provided-Basis Refactor Plan

## Goal

Make the Frobenius compiler accept paper-style basis parameters at setup time,
while keeping the current auto-search path as a compatible fallback.

## Current Status

- Implemented:
  - `FrobeniusPCSSetupInput` now supports an explicit provided-basis path.
  - `FrobeniusPCSSetup(...)` now branches between auto-discovery and
    caller-provided basis data.
  - the algebra layer validates provided basis data against an independently
    checked Frobenius context, computing a Teichmuller generator internally
    when the caller does not supply one.
  - focused Frobenius PCS tests now cover valid provided-basis setup, packing
    and commitment agreement with the auto path, end-to-end prove/verify, and
    rejection of malformed dual-basis, wrong-dimension, and broken-orbit
    inputs.
  - README/API-facing docs now describe the two setup modes and include a
    minimal caller-provided basis example.
- Remaining follow-up:
  - none for the scope of this refactor.

The target outcome is:

- callers may provide a normal basis `beta` and its dual basis `alpha` as setup
  inputs,
- setup validates the supplied basis data instead of silently trusting it,
- the existing prove/verify protocol remains unchanged,
- the current auto-discovery path still works for benches and existing tests.

## Approved Decisions

- Keep the current auto-search path for now.
- Add an explicit provided-basis setup path rather than replacing setup
  wholesale.
- In the provided-basis path, require both `beta` and `alpha` to be supplied.
- Allow an optional supplied Teichmuller generator, but still perform strong
  normality validation even when it is omitted by computing a valid generator
  during setup.
- Do not change the proof format or the backend interface in this refactor.

## Why This Refactor Is Narrow

The current Frobenius prover/verifier code only depends on:

- `params.basis_data.normal_basis`
- `params.precomputed`

Once setup can construct those two objects from caller-provided basis data,
the rest of the compiler remains the same.

## API Plan

### 1. Extend setup input

Add a new setup-side input object:

- `FrobeniusPCSProvidedBasisInput`

Planned fields:

- `NormalBasisData normal_basis`
- `bool has_teichmuller_generator`
- `NTL::ZZ_pE teichmuller_generator`

Then extend `FrobeniusPCSSetupInput` with:

- `bool use_provided_basis`
- `FrobeniusPCSProvidedBasisInput provided_basis`

Compatibility rule:

- if `use_provided_basis == false`, setup keeps the current auto-search logic,
- if `use_provided_basis == true`, setup uses the provided basis path.

### 2. Add a provided-basis constructor in the algebra layer

Add a new helper under `GaloisRing/FrobeniusBasis.*`:

- `BuildFrobeniusBasisFromProvidedNormalBasisOrThrow(...)`

This helper will:

- validate active ring contexts,
- validate basis shape and dual-basis identity,
- validate that `beta` spans the extension as a free base-ring module,
- validate normality against an independently checked Teichmuller/Frobenius
  context,
- return a `FrobeniusBasisData` object compatible with the existing compiler
  setup flow.

## Validation Rules

### Always required for provided basis

- `beta.size() == alpha.size() == 2^kappa`
- `beta` is base-ring-linearly independent and spans the extension
- `Tr(alpha_u * beta_v) = delta_uv`
- `basis_data` metadata matches the active `(base_modulus, extension_modulus,
  kappa)` context

### Strong validation for every provided-basis setup

Setup should additionally verify:

- the Teichmuller generator used for validation is valid in the active context,
- `beta` is a normal-basis Frobenius orbit under the true Frobenius action
  derived from that generator,
- `tau(beta_i) = beta_{i+1}` cyclically,
- `tau(alpha_i) = alpha_{i+1}` cyclically.

If the caller does not provide a Teichmuller generator, setup computes one and
still applies the same validation.

## Implementation Steps

### Phase A: Algebra-layer support

- extend `FrobeniusBasis.hpp` declarations,
- implement provided-basis constructor in `FrobeniusBasis.cpp`,
- factor out any shared validation helpers needed by both auto-search and
  provided-basis paths.

### Phase B: Compiler setup support

- extend `FrobeniusPCSSetupInput`,
- branch inside `FrobeniusPCSSetup(...)`,
- keep `BuildFrobeniusPCSPrecomputedTables(...)` unchanged,
- keep `ValidateFrobeniusPCSParamsOrThrow(...)` as the final common gate.

### Phase C: Tests

Add setup/proof regressions for:

- setup accepts a caller-provided valid basis,
- provided-basis setup matches the auto path on packing,
- provided-basis setup matches the auto path on commitment,
- provided-basis setup accepts an honest proof,
- setup rejects a malformed supplied dual basis,
- setup rejects a supplied basis with wrong dimension,
- setup rejects a supplied basis whose claimed Teichmuller/Frobenius orbit is
  inconsistent.

## Test Matrix

Primary test file:

- `tests/test_z2k_frobenius_pcs.cpp`

Planned new test groups:

- `TestFrobeniusPCSSetup_ProvidedBasisAcceptsValidBasis`
- `TestFrobeniusPCSSetup_ProvidedBasisRejectsBrokenDualBasis`
- `TestFrobeniusPCSSetup_ProvidedBasisRejectsWrongDimension`
- `TestFrobeniusPCSSetup_ProvidedBasisRejectsBrokenOrbitOrder`
- `TestFrobeniusPCSSetup_ProvidedBasisMatchesAutoPacking`
- `TestFrobeniusPCSSetup_ProvidedBasisMatchesAutoCommit`
- `TestFrobeniusPCSVerifyEval_AcceptsHonestProof_ProvidedBasis`

Keep green:

- `tests/test_z2k_frobenius_feasibility.cpp`
- existing Frobenius PCS tests

## Stop Rule

Do not remove the auto-search path in this refactor.

Do not weaken provided-basis setup into a trust-the-caller mode. The supplied
basis must always be checked against an independently validated Frobenius
context, whether the caller supplies the Teichmuller generator or setup derives
it internally.

## Definition of Done

This refactor is done when:

- setup accepts paper-style provided basis data,
- current auto-search setup still works unchanged,
- provided-basis and auto paths agree on packed tables and commitments for the
  same basis,
- an honest proof verifies through the provided-basis setup path,
- malformed supplied basis data is rejected by setup,
- existing Frobenius regressions stay green.
