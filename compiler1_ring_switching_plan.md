# Compiler I Ring-Switching Implementation Plan

Status: implementation tracking.

Current progress:

- `WP0` complete and committed.
- `WP1` complete and committed.
- `WP2` complete in the working tree for this turn.
- `WP3` complete in the working tree for this turn.

Branch: `feat/compiler1-ring-switching`

## Goal

Implement the paper's "Compiler I" (ring-switching based compiler) for
multilinear PCS over `Z_{2^k}` by treating the existing PCS over
`GR(2^k, r)` as a black-box backend.

Confirmed first-landing scope:

- Non-interactive single-point evaluation proofs only.
- Backend integration must be backend-agnostic from day one, even though the
  only concrete backend implementation in the first landing is the current
  BaseFold PCS over `GR(2^k, r)`.
- Keep the backend proof as an opaque subproof: the compiler layer should only
  call backend setup/commit/prove/verify entrypoints, not depend on backend
  internals.
- Proof serialization and benchmark entrypoints are part of the first landing,
  not follow-up work.

## Fixed Design Decisions

These decisions are now fixed and should not be reopened unless the paper or
the repo forces a contradiction.

1. Backend abstraction:
   - the compiler layer must expose an explicit backend-agnostic adapter
   - the initial adapter implementation may target only the current BaseFold PCS

2. Basis abstraction:
   - the API must preserve explicit basis abstraction from day one
   - the only first-landing basis implementation is
     `alpha = beta = current ZZ_pE polynomial basis`

3. Setup semantics:
   - the compiler layer must expose `Setup`
   - this first `Setup` is lightweight: the caller supplies extension-polynomial
     choice, ring representation, and backend config; `Setup` validates and
     packages them into compiler params

4. Deliverable boundary:
   - the first implementation includes library API, tests, proof serializer, and
     bench entrypoints

## Paper Alignment Anchors

The plan below is tied to the following paper locations:

- Packing from `Z_{2^k}` multilinears to `GR(2^k, r)` multilinears:
  `mineru_md/.../Polylogarithmic Proofs for Multilinears over Z2k.md:835-845`
- Compiler assumptions and the "treat Prove+Verify as an IOP" viewpoint:
  `mineru_md/.../Polylogarithmic Proofs for Multilinears over Z2k.md:883-887`
- Protocol 1 (ring-switching compiler):
  `mineru_md/.../Polylogarithmic Proofs for Multilinears over Z2k.md:264-300`
- Why the ring-switching final check can be seen as evaluating a verifier-known
  multilinear `r`:
  `mineru_md/.../Polylogarithmic Proofs for Multilinears over Z2k.md:246-262`
- Appendix C.1 proof that the `s_u` list encodes the desired partial
  evaluations:
  `mineru_md/.../Polylogarithmic Proofs for Multilinears over Z2k.md:1286-1318`
- Appendix C.1 explanation that Equality Check 3 is `g(r'_suffix) = r(r')`:
  `mineru_md/.../Polylogarithmic Proofs for Multilinears over Z2k.md:1318-1332`

The first implementation must preserve those data dependencies exactly:

1. Pack `t` into `t'`.
2. Prover sends the `s_u` list.
3. Verifier checks the claimed full evaluation `s` using recovered partial
   evaluations.
4. Prover and verifier run the batched sumcheck on
   `h(X) = tilde_g(X) * t'(X)`.
5. Prover gives one backend PCS proof for `t'(r'_suffix) = t_star`.
6. Verifier computes `g_star` from the verifier-known `r` polynomial and checks
   `s^(ell') = t_star * g_star`.

## Current Repo Mapping

What already exists:

- Backend PCS black-box API:
  [include/BaseFold/BaseFoldPCS.hpp](/home/zjl/BasefoldOverGR/include/BaseFold/BaseFoldPCS.hpp#L132)
- Multilinear evaluation and equality polynomial helpers:
  [include/BaseFold/Multilinear.hpp](/home/zjl/BasefoldOverGR/include/BaseFold/Multilinear.hpp#L13)
- Existing specialized sumcheck for `f(X) * eq_z(X)`:
  [include/BaseFold/Sumcheck.hpp](/home/zjl/BasefoldOverGR/include/BaseFold/Sumcheck.hpp#L17)
- Basic Galois-ring conversion utilities in the polynomial basis:
  [include/GaloisRing/utils.hpp](/home/zjl/BasefoldOverGR/include/GaloisRing/utils.hpp#L38)

What does not exist yet:

- A backend-agnostic compiler adapter surface for `Z_{2^k}` PCS.
- A compiler-layer proof type for `Z_{2^k}`.
- Packing helpers that explicitly build `t'` from `t` using the compiler's
  chosen basis.
- Verifier-side decomposition helpers that turn
  `eq(r_kappa, ..., r_{ell-1}; w)` into the coefficients `A_{u||w}`.
- A generic sumcheck for the product of two arbitrary multilinears.
- A transcript schedule, proof serializer, and benchmark entrypoints for the
  outer compiler proof.

## Basis Strategy For First Landing

For Compiler I, the Appendix C.1 proof only requires:

- one basis `alpha = {alpha_u}` used to decompose the suffix equality vector,
  and
- one basis `beta = {beta_v}` used for packing and for decomposing `s_u`.

Appendix C.1 explicitly allows `alpha` to be different from `beta`, but it does
not require the Frobenius/trace machinery that Compiler II needs.

That means the simplest paper-consistent first landing is:

- choose `alpha = beta =` the current polynomial basis induced by the active
  `ZZ_pE` modulus;
- do all coefficient extraction by polynomial-basis decomposition;
- postpone dual-basis / trace / normal-basis work to Compiler II.

This is the main reason Compiler I is implementable now while Compiler II is not
yet cheap to land.

API consequence:

- the public compiler params should carry basis handles or descriptors for
  `alpha` and `beta`
- the first concrete basis provider should simply bind both handles to the
  active polynomial basis exposed by the current `ZZ_pE` modulus

## Proposed Public Surface

New module set (proposed names, exact placement to be confirmed):

- `include/BaseFold/Z2kPCSBackend.hpp`
- `include/BaseFold/Z2kRingSwitchPCS.hpp`
- `include/BaseFold/Z2kRingSwitchProofSerialize.hpp`
- `src/BaseFold/Z2kRingSwitchPCS.cpp`
- `src/BaseFold/Z2kRingSwitchProofSerialize.cpp`
- `tests/test_z2k_ring_switch_pcs.cpp`
- `bench/bench_z2k_ring_switch_eval.cpp`
- `bench/bench_z2k_ring_switch_proof_size.cpp`

Proposed API shape, mirroring the current BaseFold style:

- `RingSwitchPCSSetup(...)`
- `RingSwitchPCSCommit(...)`
- `RingSwitchPCSBuildCommitArtifacts(...)`
- `RingSwitchPCSProveEval(...)`
- `RingSwitchPCSVerifyEval(...)`
- `RingSwitchPCSEvalProofSizeBytes(...)`

Backend adapter shape:

- `Z2kPCSBackend` interface with setup/commit/prove/verify hooks
- `BaseFoldZ2kPCSBackend` concrete adapter for the current BaseFold PCS

Internal note:

- The paper's abstract syntax includes `Setup`, `Commit`, `Open`, `Prove`,
  `Verify`.
- The repo's current BaseFold surface does not expose the paper's explicit
  `Open`; it uses `Commit` plus optional cached commit artifacts.
- To stay consistent with the repo, the compiler layer should follow the
  existing BaseFold style and treat "open hint" as an internal cached artifact,
  not as a public API.
- The compiler's public polynomial input should follow the paper and be a
  Boolean-hypercube / Lagrange-basis table.
- The current backend core should stay monomial-basis for now.
- Therefore the compiler/backend facade must explicitly convert Boolean tables
  into monomial coefficients before invoking the backend black box.

## Data Types To Add

### Compiler params

`RingSwitchPCSParams`

- `long ell`
- `long kappa`
- `long ell_prime`
- basis descriptors/handles for `alpha` and `beta`
- backend kind / adapter vtable / backend-owned params blob
- an optional backend challenge config hook if we want to support the existing
  extension-challenge path later

Validation invariants:

- `kappa >= 1`
- `ell >= kappa`
- backend point dimension is exactly `ell_prime = ell - kappa`
- backend message length is exactly `2^ell_prime`
- current `ZZ_p` modulus is `2^k`
- current `ZZ_pE` degree is exactly `r = 2^kappa`
- `alpha` and `beta` basis dimensions both equal `r`
- first implementation rejects any basis configuration other than the active
  polynomial basis for both `alpha` and `beta`

### Commit artifacts

`RingSwitchPCSCommitArtifacts`

- packed coefficients `t_packed_coeffs`
- backend-opaque commit artifact blob
- backend commitment object/root

### Proof object

`RingSwitchPCSEvalProof`

- `std::vector<NTL::ZZ_pE> s_by_u`
- `std::vector<basefold::QuadraticPoly> h_by_round`
- `NTL::ZZ_pE packed_eval_at_rprime`
- backend-opaque eval subproof

No extra hidden state should be required for verification.

## Work Packages

### WP0: Lightweight setup, parameter, and basis plumbing

Deliverables:

- backend-agnostic adapter interface and BaseFold adapter implementation
- `RingSwitchPCSSetup`
- compiler param validator
- helper that checks the current NTL contexts match the compiler assumptions
- helper that interprets the active `ZZ_pE` polynomial basis as the compiler
  basis

Implementation notes:

- `RingSwitchPCSSetup` should be lightweight, not fully automatic: the caller
  supplies the extension polynomial / active NTL context / backend config, and
  setup validates and packages them.
- The first landing should not add new algebraic setup code for choosing bases.
- It should reuse the active `ZZ_pE` polynomial basis as both `alpha` and
  `beta`, while still preserving basis abstraction in the params/API.

Exit criteria:

- construction fails loudly when `ZZ_p::modulus()` is not a power of two
- construction fails loudly when `deg(ZZ_pE::modulus()) != 2^kappa`
- construction fails loudly when the requested backend adapter is missing
- construction fails loudly when the caller requests any non-polynomial basis

### WP1: Packing and decomposition utilities

Deliverables:

- `PackZ2kCoeffsToGREvals(...)`
- `BooleanHypercubeTableToMonomialCoeffs(...)`
- `DecomposeGRElementToBaseCoeffsPolynomialBasis(...)`
- `BuildRingSwitchComponentTensor(...)`
- basis-aware wrappers whose first implementation dispatches only to the
  polynomial-basis helpers

Planned behavior:

1. Packing:
   - input: Boolean-hypercube / Lagrange-basis table of `t` over `Z_{2^k}` of length
     `2^ell`
   - output: packed Boolean-hypercube / Lagrange-basis table of `t'` over
     `GR(2^k, r)` of
     length `2^ell_prime`
   - formula must match paper Eq. (1) / packing definition

2. Suffix equality decomposition:
   - for each `w in {0,1}^{ell_prime}`, compute
     `eq(r_kappa, ..., r_{ell-1}; w)` in `GR`
   - decompose it in the polynomial basis into base-ring coefficients
     `A_{u||w} in Z_{2^k}`

3. Verifier-known polynomial `r`:
   - flatten the `A_{u||w}` table into the full `2^ell` Boolean-hypercube table
     of the verifier-known multilinear `r`
   - additionally compute the monomial-basis coefficients of `r` so existing
     backend-core and multilinear evaluation utilities can be reused

Exit criteria:

- packing round-trips on small examples
- computed `A_{u||w}` satisfy the Appendix C.1 coefficient identities
- Boolean-table to monomial conversion is covered by tests, because later
  backend-facade calls will depend on it

### WP2: Generic product sumcheck

Deliverables:

- a new prover for the product of two arbitrary multilinears
- a matching verifier-side chain checker

Reason:

- the current `SumcheckProver` only supports `f(X) * eq_z(X)`
- Compiler I needs `h(X) = tilde_g(X) * t'(X)` where `tilde_g` is verifier-known
  but otherwise arbitrary

Planned interface:

- `ProductSumcheckProver(f_coeffs, g_coeffs)`
- `CheckProductSumcheckChain(initial_claim, h_by_round, r_suffix)`

Planned algorithm:

- keep boolean evaluation tables for both multilinears
- at each round, for each remaining assignment pair `(f0, f1)` and `(g0, g1)`,
  accumulate the degree-2 polynomial
  `(f0 + (f1-f0) X) * (g0 + (g1-g0) X)`
- after a verifier challenge `r_i`, fold both tables with the same challenge

This utility should be written to be reusable by Compiler II later.

Exit criteria:

- sumcheck passes on self-consistent small examples
- tampering one `h_i` causes verification failure

### WP3: Compiler commit path

Deliverables:

- `RingSwitchPCSCommit`
- `RingSwitchPCSBuildCommitArtifacts`

Planned behavior:

1. Validate params and dimensions.
2. Pack `t` into `t'`.
3. Call backend adapter `Commit` or `BuildCommitArtifacts` on `t'`.
4. Return the backend commitment as the compiler commitment.

Paper alignment:

- this must exactly implement Protocol 1 commit step:
  "fix the packed polynomial `t'` and output `Commit'(t')`"

Exit criteria:

- compiler commitment matches direct backend commitment on the packed polynomial

### WP4: Compiler prove path

Deliverables:

- `RingSwitchPCSProveEval`
- optional unchecked/internal variant if needed by the implementation

Planned behavior:

1. Validate dimensions and recompute the claimed original evaluation `s`.
2. Pack the Boolean table of `t` into the Boolean table of `t'`.
3. Build the `A_{u||w}` table from the suffix part of the query point.
4. Compute the prover messages
   `s_u = sum_w A_{u||w} * t'(w)` for all `u`.
5. Recover the partial evaluations from the coefficient matrix exactly as in
   Appendix C.1 and check the paper's Equality Check 1 locally.
6. Use transcript-derived batching challenges
   `r'_0, ..., r'_{kappa-1}` to form `g(w)`.
7. Run the generic product sumcheck on `h(X) = tilde_g(X) * t'(X)`.
8. Let `r'_suffix = (r'_kappa, ..., r'_{ell-1})`.
9. Convert the packed Boolean table of `t'` into packed monomial coefficients.
10. Compute `t_star = t'(r'_suffix)`.
11. Call backend adapter `ProveEval` on the packed monomial coefficients and
    `(r'_suffix, t_star)`.
12. Output the outer proof carrying `s_by_u`, `h_by_round`, `t_star`, and the
    opaque backend proof.

Important implementation note:

- Step 5 is not optional bookkeeping. It is the easiest way to keep the code
  aligned with Appendix C.1 and to catch basis-indexing bugs early.

### WP5: Compiler verify path

Deliverables:

- `RingSwitchPCSVerifyEval`

Planned behavior:

1. Recompute the transcript-derived challenges from public input and proof
   messages.
2. Decompose each `s_u` in the `beta` basis and recover the partial evaluations
   exactly as in Appendix C.1 Eq. (14) and Eq. (15).
3. Run Equality Check 1:
   `s == sum_v partial_eval[v] * eq(prefix_query; v)`.
4. Recompute `s^(0)` from the `s_u` list and the batching scalars
   `eq(r'_0, ..., r'_{kappa-1}; u)`.
5. Verify the generic sumcheck chain.
6. Rebuild the verifier-known polynomial `r` from `A_{u||w}`.
7. Compute `g_star` as `r(r'_0, ..., r'_{ell-1})`, using the Appendix C.1
   observation that this equals `tilde_g(r'_suffix)`.
8. Check `s^(ell_prime) == t_star * g_star`.
9. Verify the backend proof on point `r'_suffix` and value `t_star` through the
   backend adapter.

Exit criteria:

- all checks pass on honest proofs
- each of the following tamperings fails independently:
  wrong `s`, wrong `s_u`, wrong `h_i`, wrong `t_star`, wrong backend proof

### WP6: Tests

Minimum required tests:

1. Packing test:
   - direct packing matches manual small examples

2. Component-table test:
   - `A_{u||w}` decomposition reconstructs the suffix equality vector

3. Matrix-recovery test:
   - `s_u` decomposed over `beta` recovers the expected partial evaluations
     exactly as in Appendix C.1

4. End-to-end compiler test over `GR(4, 2)`:
   - choose `kappa = 1`, small `ell`
   - prove and verify an honest claim

5. Negative end-to-end tests:
   - wrong original claim `s`
   - mutated `s_u`
   - mutated sumcheck polynomial
   - mutated backend subproof

Recommended first context:

- `Z_4` multilinears evaluated at points in `GR(4, 2)`
- this keeps `kappa = 1` and makes manual debugging manageable

### WP7: Serializer and benchmark follow-up

Deliverables:

- serializer for the outer compiler proof
- proof-size accounting for the composed proof
- benchmark entrypoints for:
  - packing time
  - outer prover time excluding backend prove
  - outer verifier time excluding backend verify
  - total composed proof size

Implementation notes:

- the serializer must treat the backend subproof as an opaque payload routed
  through the selected backend adapter
- benchmark binaries should report both outer-only and end-to-end composed
  costs, because the compiler layer is meant to wrap interchangeable backends

## Suggested Implementation Sequence

Implement in this order:

1. WP0
2. WP1
3. WP2
4. WP3
5. WP4
6. WP5
7. WP6
8. WP7

WP7 should start only after the end-to-end correctness tests in WP6 are green,
but it is still part of the first landing.

## Risks and Failure Modes

High-risk areas:

- bit ordering mismatch between the paper's `v || w` concatenation and the repo's
  multilinear coefficient ordering
- silently using the wrong basis when decomposing `s_u`
- transcript schedule drift between prover and verifier
- backend adapter becoming nominally generic but still leaking BaseFold-specific
  assumptions into compiler code
- attempting to reuse the current specialized BaseFold sumcheck instead of
  adding a generic product-sumcheck utility

Recommended defensive checks:

- after computing `s_u`, reconstruct the partial-evaluation matrix explicitly in
  debug/test paths
- in verifier code, assert every dimension derived from `ell`, `kappa`, and
  backend params before doing algebra
- keep transcript labels disjoint from the backend PCS transcript labels
- in adapter tests, run the same compiler tests through the public adapter
  surface rather than directly calling BaseFold helpers
