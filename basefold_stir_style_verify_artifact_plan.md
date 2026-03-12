# BaseFold STIR-Style Verifier Artifact Plan

## Status

- Owner: TBD
- Status: Phase 6 complete
- Last updated: 2026-03-12

## Current Progress

- [x] Phase 1: fixed-width `BaseFoldPCSEvalProof` deserializer and round-trip tests
- [x] Phase 2 helper layer: `bench/bench_pcs_artifact_common.hpp`
  - deterministic short hash `artifact_id`
  - human-readable `display_key`
  - `manifest.jsonl` append/load helpers
  - `meta.json` read/write helpers
  - combined `public_inputs.bin` binary round-trip helpers
  - metadata-driven verification-context restore helper
- [x] Phase 3: `dump_pcs_eval_artifact`
- [x] Phase 4: `bench_pcs_verify_artifact`
- [x] Phase 5: tests
- [x] README artifact-flow documentation

## Goal

Add a STIR-style artifact workflow for BaseFold verifier benchmarking:

- generate real BaseFold eval-proof artifacts once
- store multiple parameterized cases under one artifact root
- load one selected case into memory
- exclude file read and proof deserialization from `verifier mean`
- benchmark only `BaseFoldPCSVerifyEval(...)`

This plan must preserve the current benchmark split:

- `bench_pcs_commit` stays responsible for top commit timing
- `bench_pcs_prove` stays responsible for prove-only timing
- `bench_pcs_verify` stays as the self-contained verifier benchmark
- `bench_pcs_eval` stays as prove+verify benchmark

The new artifact path is additive rather than a replacement.

## Settled Design Decisions

- BaseFold proof deserialization is allowed for this work. The old ring-switch-only `serialize but not deserialize` choice is not treated as a repo-wide default.
- Artifact storage is one root directory containing many parameterized cases.
- Each artifact tool invocation operates on exactly one case. No `--mode both`.
- `bench_pcs_verify_artifact` excludes file read and deserialization from measured verifier time.
- Artifacts store only `commitment_root`, not full `BaseFoldPCSCommitArtifacts`.
- Artifact public inputs are stored as one combined `public_inputs.bin`.
- Artifact CLI supports explicit `--artifact-id`; if omitted, a canonical id is auto-generated.
- Artifact metadata stores the actual `zeta_coeffs` used, even if the case was created with `--auto-zeta teich`.

## Non-Goals

- Do not change the semantics of existing `bench_pcs_*` binaries.
- Do not add end-to-end "load from artifact and re-prove" flows.
- Do not store full top-level committed oracle or `MerkleTree` internals in artifacts.
- Do not mix multiple cases into one proof file.
- Do not count file IO or deserialization inside verifier timing.

## High-Level Architecture

### New protocol-layer capability

Add fixed-width deserialization support for `BaseFoldPCSEvalProof`, matching the existing fixed-width serializer contract.

### New benchmark-layer capability

Add a benchmark artifact container that stores:

- proof bytes
- verifier public inputs
- enough metadata to reconstruct the verification context

### New binaries

- `dump_pcs_eval_artifact`
- `bench_pcs_verify_artifact`

## Artifact Layout

Artifact root contains a manifest plus an object store.

Example:

```text
<artifact-root>/
  manifest.jsonl
  objects/
    bfv_a13c9e7f82d1f4a0/
      meta.json
      public_inputs.bin
      proof.bin
    bfv_83d6e1c1a4b2f915/
      meta.json
      public_inputs.bin
      proof.bin
```

The manifest is the primary index. Object paths use a short hash-based
`artifact_id`; the full human-readable context lives in `manifest.jsonl` and the
per-object `meta.json` rather than in directory names.

## Artifact Case Contents

### `meta.json`

Required fields:

- `artifact_id`
- `display_key`
- `context_id`
- `context_label`
- `mode`
- `c`
- `k0`
- `d`
- `poly_dim`
- `lambda`
- `gamma`
- `queries`
- `seed`
- `use_checked_prover_path`
- `use_extension_challenges`
- `scalar_modulus`
- `base_prime` for ring mode
- `F_coeffs`
- `zeta_coeffs`
- `zeta_source` as `explicit` or `auto_teich`
- `challenge_extension_coeffs` when extension challenges are enabled
- `hash_backend`
- `proof_encoding` and version info
- `proof_size_bytes`

### `public_inputs.bin`

Combined verifier public inputs:

- `commitment_root`
- point `z`
- claimed value `y`

### `proof.bin`

Fixed-width serialized `BaseFoldPCSEvalProof`.

## CLI Design

### `dump_pcs_eval_artifact`

Purpose:

- build a real BaseFold proof once
- serialize and persist one artifact case

Accepted modes:

- `--mode field`
- `--mode ring`

Explicitly rejected:

- `--mode both`

Suggested CLI:

```text
dump_pcs_eval_artifact
  --artifact-root <dir>
  [--artifact-id <id>]
  --mode field|ring
  [--c <int>] [--k0 <int>] [--d <int>] [--queries <int>]
  [--checked] [--seed <u64>]
  [--use-extension-challenges]
  [--field-mod ...] [--field-F ...] [--field-zeta ...]
  [--ring-mod ...] [--ring-p ...] [--ring-F ...] [--ring-zeta ...]
  [--field-challenge-ext ...] [--ring-challenge-ext ...]
  [--auto-zeta teich]
```

Behavior:

- derive deterministic `f_coeffs`, `z`, and `y` from parameters and seed
- build top commit artifacts in memory
- generate one real `BaseFoldPCSEvalProof`
- write one case directory under `artifact-root`
- print artifact id, case path, and proof size

### `bench_pcs_verify_artifact`

Purpose:

- load one artifact case
- deserialize outside the timed region
- benchmark pure verify

Suggested CLI:

```text
bench_pcs_verify_artifact
  --artifact-root <dir>
  --artifact-id <id>
  [--warmup <int>] [--reps <int>] [--profile]
  [--merkle-leaves-per-thread <int>] [--merkle-level-threshold <int>] [--merkle-max-threads <int>]
  [--verifier-query-per-thread <int>] [--verifier-query-threshold <int>] [--verifier-query-max-threads <int>]
```

Behavior:

- read `meta.json`
- restore `ZZ_p`, `ZZ_pE`, and challenge-extension context
- read and deserialize `public_inputs.bin`
- read and deserialize `proof.bin`
- start timing only after all of the above succeeds
- warm up and benchmark repeated `BaseFoldPCSVerifyEval(...)`
- print verifier timing plus artifact metadata summary

## File Plan

### New files

- `include/PCS/BaseFold/ProofDeserialize.hpp`
- `src/PCS/BaseFold/ProofDeserialize.cpp`
- `bench/bench_pcs_artifact_common.hpp`
- `bench/dump_pcs_eval_artifact.cpp`
- `bench/bench_pcs_verify_artifact.cpp`
- `tests/test_pcs_artifact.cpp`

### Expected modified files

- `CMakeLists.txt`
- `README.md`
- optionally `README_zh.md`

## Protocol-Layer Work

### Phase 1: proof deserialization

Add BaseFold proof deserialization matching the existing fixed-width serializer.

Required tasks:

- define reader utilities for fixed-width byte decoding
- deserialize `IOPPMerkleCommitments`
- deserialize `h_by_level`
- deserialize `pi0_codeword`
- deserialize base `query_multiproofs`
- deserialize optional extension payload
- validate structure lengths against encoded counts
- validate extension-width assumptions against metadata-provided options

Deliverables:

- `ProofDeserialize.hpp`
- `ProofDeserialize.cpp`
- unit tests covering base and extension cases

Acceptance criteria:

- serialize -> deserialize -> verify succeeds for base mode
- serialize -> deserialize -> verify succeeds for extension mode

## Benchmark Artifact Layer

### Phase 2: artifact schema and common helpers

Add benchmark-only artifact helpers in `bench/bench_pcs_artifact_common.hpp`.

Required tasks:

- define metadata struct
- define manifest-entry struct
- define public-inputs struct
- define canonical artifact-id generator
- define human-readable display-key generator
- add manifest append/load helpers for `manifest.jsonl`
- add JSON read/write helpers for metadata
- add binary read/write helpers for public inputs
- add helper to reconstruct verification context from metadata
- add helper to compute and print artifact summary

Canonical artifact id should be a stable hash over case-defining metadata rather
than a long path-encoded parameter string.

The hashed input should include at least:

- `context_id`
- `context_label`
- `mode`
- `c`
- `k0`
- `d`
- `poly_dim`
- `lambda`
- `gamma`
- `queries`
- `seed`
- extension-challenge marker

The human-readable `display_key` should carry the same core context in a compact
printable form for logs, manifest inspection, and CLI selection output.

Acceptance criteria:

- generated artifact ids are deterministic
- generated display keys are deterministic and human-readable
- `manifest.jsonl` is sufficient as the primary discovery/index surface
- metadata is sufficient to reconstruct verifier context without extra CLI protocol parameters

### Phase 3: artifact dump tool

Implement `dump_pcs_eval_artifact`.

Required tasks:

- reuse current deterministic input generation from benchmark helpers
- reuse current prove path selection including checked/unchecked and extension-challenge modes
- serialize proof to `proof.bin`
- serialize public inputs to `public_inputs.bin`
- write `meta.json`
- reject overwrite unless an explicit overwrite flag is added
- reject `--mode both`

Acceptance criteria:

- dumping one field case produces a complete case directory
- dumping one ring case produces a complete case directory
- `proof_size_bytes` in `meta.json` matches serializer output

Update 2026-03-12:

- `dump_pcs_eval_artifact` is now implemented in `bench/dump_pcs_eval_artifact.cpp`.
- The tool remains single-case only and rejects `--mode both`.
- Phase 3 currently rejects overwriting existing `artifact_id` directories; no `--overwrite` flag has been added.
- Targeted integration coverage now dumps one field case and one ring+extension case, then reloads `meta.json`, `public_inputs.bin`, and `proof.bin` to verify the stored proof.

### Phase 4: artifact-based verifier benchmark

Implement `bench_pcs_verify_artifact`.

Required tasks:

- load one selected case by `artifact-id`
- rebuild field/ring context from metadata
- deserialize public inputs and proof
- perform warmup runs
- measure only verify loop
- print `artifact load wall time` and `artifact deserialize wall time` as non-headline diagnostics
- print headline `verifier mean` and `input proof size`
- preserve existing anti-optimization sink pattern

Acceptance criteria:

- measured verifier timing excludes file IO and deserialization
- verifier succeeds on dumped artifacts
- output clearly distinguishes measured verifier time from non-measured load time

Update 2026-03-12:

- `bench_pcs_verify_artifact` is now implemented in `bench/bench_pcs_verify_artifact.cpp`.
- The tool loads one case by `artifact-id`, restores the verification context from metadata, deserializes `public_inputs.bin` and `proof.bin` outside the timed region, and times only repeated verify calls.
- Output now includes `artifact load wall time` and `artifact deserialize wall time` as explicitly excluded diagnostics, plus headline `verifier mean` and `input proof size`.
- `--merkle-*` flags are intentionally unsupported here because artifact verify does not build Merkle trees inside the measured path; only `--verifier-query-*` affects the timed verifier loop.

## Testing Plan

### Unit tests

- base proof fixed-width round-trip
- extension proof fixed-width round-trip
- public-inputs binary round-trip
- metadata JSON round-trip

### Integration tests

- dump one field artifact and verify it through artifact loader
- dump one ring artifact and verify it through artifact loader
- compare `bench_pcs_verify` and `bench_pcs_verify_artifact` under the same logical case

### Negative tests

- tamper one byte in `proof.bin` and verify fails
- tamper `commitment_root` and verify fails
- tamper point `z` and verify fails
- tamper claimed value `y` and verify fails
- tamper challenge extension modulus metadata and verify fails in extension mode

Update 2026-03-12:

- Phase 5 regression coverage now lives in `tests/test_pcs.cpp`.
- The direct-vs-artifact comparison test uses one shared logical field case and asserts stable equivalence of verifier inputs and proof bytes, rather than timing closeness.
- Negative coverage now exercises on-disk tampering of `proof.bin`, `commitment_root`, point `z`, and claimed `y`, plus extension-metadata tampering after load.

## Benchmark Validation

Run at least:

```bash
./build-release/dump_pcs_eval_artifact --artifact-root /tmp/basefold-artifacts --mode field --d 10 --queries 2 --warmup 0 --reps 1
./build-release/bench_pcs_verify_artifact --artifact-root /tmp/basefold-artifacts --artifact-id <id> --warmup 0 --reps 3
./build-release/bench_pcs_verify --mode field --d 10 --queries 2 --warmup 0 --reps 3
```

Validation expectations:

- artifact verify succeeds
- `input proof size` matches across artifact and non-artifact flows
- verifier mean is directionally consistent with current self-contained verify benchmark

## Documentation Updates

Update README benchmark section to explain two verifier-only modes:

- self-contained `bench_pcs_verify`
- artifact-driven `bench_pcs_verify_artifact`

Document clearly:

- artifact-driven verifier timing excludes file IO and deserialization
- artifact tools are single-case only
- artifact metadata stores actual `zeta_coeffs`
- artifacts store only `commitment_root`, not full commit artifacts

Update 2026-03-12:

- `README.md` and `README_zh.md` now document both verifier-only modes:
  self-contained `bench_pcs_verify` and artifact-driven
  `bench_pcs_verify_artifact`.
- The README benchmark section now includes the artifact dump/load workflow,
  its on-disk layout, exact timing semantics, and the distinction between
  `input proof size` reporting and headline verifier timing.
- The docs also make explicit that artifact tools are single-case only,
  store actual `zeta_coeffs` in metadata, and store only `commitment_root`
  rather than full top-level commit artifacts.

## Risks

- deserializer bugs can silently misparse proof layout if field-width assumptions are wrong
- extension-challenge cases are sensitive to modulus reconstruction and width checks
- metadata drift can cause artifact load failures if CLI defaults are silently reused instead of metadata truth
- overly broad artifact ids can cause collisions across cases

## Mitigations

- make round-trip tests mandatory before wiring benchmark CLI
- keep proof serializer and deserializer on one fixed-width contract
- reconstruct all verifier-critical parameters from metadata, not from current CLI defaults
- include enough case-defining fields in the canonical artifact id

## Implementation Order

- [x] Phase 1: BaseFold proof deserialization
- [x] Phase 2: artifact schema and shared helpers
- [x] Phase 3: `dump_pcs_eval_artifact`
- [x] Phase 4: `bench_pcs_verify_artifact`
- [x] Phase 5: tests
- [x] Phase 6: README updates

## Done Criteria

This plan is complete when all of the following are true:

- one command can dump a single BaseFold eval-proof artifact case
- one command can benchmark verifier time from a dumped artifact case
- artifact verifier timing excludes file read and deserialization
- base and extension proof round-trip tests pass
- README documents the artifact flow and its timing semantics
