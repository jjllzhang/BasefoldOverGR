# Prover Parallel Optimization Plan

Date: 2026-03-17

Status: in progress

## 0. Current Status Update

Completed in the current workspace:

- Workstream A is complete:
  - a dedicated profile bucket now exists for `Msg0CoeffsAtSuffixChallenges(...)`
  - on `ring-gr-2p16-128-ext`, `c=4`, `d=12`, `queries=2048`, `OMP_NUM_THREADS=1`, the old `Other` bucket was reduced from about `1072 ms` to about `55 ms`
- Workstream D foundation is complete:
  - prover commit-round parallel thresholds are now runtime-configurable
  - env vars and bench CLI overrides exist for both base and extension commit rounds
  - extension commit-round threshold tuning has now started, with the current default plus auto-tuned ring behavior adjusted toward realistic parallel activation
- Workstream B primary target is complete:
  - the prover no longer recomputes `msg0_coeffs` by refolding the original `f_coeffs`
  - it now derives `msg0_coeffs` from the remaining `ExtensionSumcheckProver` state by converting the surviving Boolean evaluation table back to monomial coefficients
- Workstream E is now partially complete:
  - `CurrentPolynomialSmall(...)` and generic `CurrentPolynomial()` now use deterministic thread-local reductions
  - `ReceiveChallengeSmall(...)` and generic `ReceiveChallengeGenericState(...)` now parallelize the independent half-table folds
  - the remaining init path now also parallelizes the table lift plus `prefix_eq_by_vars_` expansion
  - these changes only affect execution strategy; transcript order, challenge derivation, and proof layout remain unchanged

Measured result after the `msg0_coeffs` dedup change:

- `ring-gr-2p16-128-ext`, `c=4`, `d=12`, `queries=2048`, `OMP_NUM_THREADS=1`
  - `prove-phase mean`: about `6075.8 ms -> 5067.8 ms`
  - `Msg0CoeffsAtSuffixChallenges`: about `1030.5 ms -> 0.006 ms`
  - proof size remained unchanged

Not started yet:

- Workstream C (`ExtensionMerkleTree::Build(...)` parallelization)
- any further `ExtensionSumcheckProver` tuning beyond the current loop + init parallelization
- any follow-up tuning of default prover commit-round thresholds after real thread sweeps

Revised priority after the completed work:

1. Revisit the dominant ring context and continue with `ExtensionCommitRound` residual cost before `ExtensionMerkleTree::Build(...)`.
2. Keep `ExtensionMerkleTree::Build(...)` later unless a new profile shows it overtaking commit-round cost on the ring target or unless field-first latency becomes the goal.
3. Leave deeper sumcheck follow-up for later unless a new thread/threshold pass exposes another clean win.

Current tuning evidence:

- `ring-gr-2p16-128-ext`, `c=4`, `d=12`, `queries=2048`, `OMP_NUM_THREADS=8`
  - older sweep before the sumcheck parallelization had already shown that realistic thresholds matter:
    - `ext-per-thread=4096`: `prove ~4638 ms`, `ExtensionCommitRound ~2644 ms`
    - `ext-per-thread=2048`: `prove ~3986 ms`, `ExtensionCommitRound ~1895 ms`
    - `ext-per-thread=1024`: `prove ~3207 ms`, `ExtensionCommitRound ~1218 ms`
    - `ext-per-thread=512`: `prove ~2836 ms`, `ExtensionCommitRound ~825 ms`
  - refreshed sweep after the sumcheck work:
    - `ext-per-thread=256`: `prove ~1638 ms`, `ExtensionCommitRound ~612 ms`
    - `ext-per-thread=128`: `prove ~1330 ms`, `ExtensionCommitRound ~507 ms`
    - `ext-per-thread=64`: `prove ~1177 ms`, `ExtensionCommitRound ~475 ms`
  - current default path now auto-selects the ring-equivalent `64` threshold unless the user explicitly overrides it:
    - default config: `prove ~1142 ms`, `ExtensionCommitRound ~459 ms`
- `field-prime128-ext`, `c=4`, `d=12`, `queries=2048`, `OMP_NUM_THREADS=8`
  - refreshed sweep after the sumcheck work:
    - `ext-per-thread=256`: `prove ~44.3 ms`, `ExtensionCommitRound ~9.2 ms`
    - `ext-per-thread=128`: `prove ~39.4 ms`, `ExtensionCommitRound ~6.6 ms`
    - `ext-per-thread=64`: `prove ~40.9 ms`, `ExtensionCommitRound ~7.6 ms`
  - current default path keeps the field side on the `128` threshold; single-shot profile is still about:
    - default config: `prove ~39.9 ms`, `ExtensionCommitRound ~6.9 ms`

Current default choice:

- keep `base_elements_per_thread=4096`
- keep `ext_elements_per_thread=128` as the explicit cross-context default
- when that value is left unmodified, auto-tune the effective threshold down to `64` on composite-modulus ring contexts
- explicit env/CLI/process overrides disable that auto-tuning and use the requested threshold verbatim

Current sumcheck scaling evidence after the new parallelization:

- `ring-gr-2p16-128-ext`, `c=4`, `d=12`, `queries=2048`
  - `OMP_NUM_THREADS=1`: `prove ~4992 ms`, `ExtensionSumcheck ~1605 ms`
  - `OMP_NUM_THREADS=8`, current default path with ring auto-tuning: `prove ~1150 ms`, `ExtensionSumcheck ~342 ms`
  - inside the 8-thread profile:
    - `CurrentPolynomial ~129 ms`
    - `ReceiveChallenge ~152 ms`
    - `init ~62 ms`
    - `ExtensionCommitRound ~459 ms`
    - `ExtensionMerkleTree::Build ~208 ms`
- `field-prime128-ext`, `c=4`, `d=12`, `queries=2048`, `OMP_NUM_THREADS=8`
  - `prove ~39.4 ms`
  - `ExtensionSumcheck ~5.0 ms`
  - `ExtensionCommitRound ~6.6 ms`
  - `ExtensionMerkleTree::Build ~12.6 ms`

## 1. Goal And Scope

This plan targets lower `prove-phase mean` for the current BaseFold PCS prover hot path, with emphasis on the extension-challenge prover used by:

- `field-prime128-ext`
- `ring-gr-2p16-128-ext`
- other extension-challenge contexts that share the same `BaseFoldPCSProveEvalWithChallengeConfig...` code path

Primary metric:

- `bench_basefold_pcs_prove` `prove-phase mean`

Secondary metrics:

- `bench_basefold_pcs_eval` `prove-phase mean`
- relevant profile buckets from `--profile --warmup 0 --reps 1`

Out of scope for this plan:

- verifier-query parallelism
- top commit / encoder throughput as measured by `bench_basefold_pcs_commit`
- changing transcript order or proof semantics

Important scope reminder:

- current `bench_basefold_pcs_prove` starts timing after `BaseFoldPCSBuildCommitArtifactsUnchecked(...)`, so top-level `EncodeFoldableUnchecked`, base Merkle commit, base sumcheck precomputation, extension commit-round precomputation, and lifted top-oracle caching are already outside the timed prover window

## 2. Current Evidence

### 2.1 Timed prover really excludes top commit

`bench_basefold_pcs_prove` builds `commit_artifacts` before starting the timed loop, so current prover optimization work should focus on the prove path itself, not the top commit path.

Relevant code:

- `bench/bench_basefold_pcs_prove.cpp`
- `src/PCS/BaseFold/BaseFoldPCSCommit.cpp`

### 2.2 Current timed prove is dominated by extension-path work

Measured with `--profile --warmup 0 --reps 1`.

Field, `field-prime128-ext`, `d=12`, `queries=2048`, `OMP_NUM_THREADS=1`:

- `BaseFoldPCSProveEval total`: about `79.4 ms`
- `ExtensionSumcheck total`: about `18.3 ms`
- `ExtensionCommitRound`: about `20.4 ms`
- `ExtensionMerkleTree::Build`: about `12.2 ms`
- `ExtensionMerkleTree::Open`: about `5.1 ms`
- `Other`: about `18.7 ms`

Ring, `ring-gr-2p16-128-ext`, `d=12`, `queries=2048`, `OMP_NUM_THREADS=1`:

- `BaseFoldPCSProveEval total`: about `4411.2 ms`
- `ExtensionSumcheck total`: about `1627.8 ms`
- `ExtensionCommitRound`: about `1551.6 ms`
- `ExtensionMerkleTree::Build`: about `110.7 ms`
- `ExtensionMerkleTree::Open`: about `41.8 ms`
- `Other`: about `1072.3 ms`

The extension prover buckets dominate. Base-path prover buckets are effectively zero in this workload because the extension-challenge path is the hot path.

### 2.3 Current thread sweep shows essentially no prover speedup

For `ring-gr-2p16-128-ext`, `d=12`, `queries=2048`:

- `OMP_NUM_THREADS=1`: about `4351.7 ms`
- `OMP_NUM_THREADS=8`: about `4361.7 ms`

The timed prove path currently does not benefit from threads in any meaningful way.

### 2.4 Existing prover parallelism is present but mostly inactive at common sizes

`ProverCommitRoundNoValidate(...)` and `ProverCommitRoundExtensionNoValidate(...)` already call `ForEachIndexMaybeParallel(...)`, but:

- they hard-code `kParallelThreshold = 4096`
- `ForEachIndexMaybeParallel(...)` only spawns threads when `threads_to_use = work_items / parallel_threshold >= 2`

This means:

- `work_items = 4096` still runs single-threaded
- for `c=2`, `k0=1`, `d=12`, the largest prover-commit round has `n_i = 4096`, so none of the commit rounds go parallel
- even at somewhat larger `d`, only the largest one or two rounds go parallel

So the existing OpenMP hooks are real, but the current threshold policy suppresses them for many relevant prover workloads.

### 2.5 Extension Merkle build is still serial

Base Merkle build already has an OpenMP implementation in `src/PCS/Common/Merkle.cpp`.

Extension Merkle build in `src/PCS/BaseFold/BaseFoldPCSExtension.cpp` is still fully serial:

- serial leaf hashing loop
- serial internal-node loop

This is a clean candidate for low-risk parallelization.

### 2.6 Extension sumcheck is still serial

The main hot loops in `ExtensionSumcheckProver` are serial:

- `CurrentPolynomialSmall(...)`
- `ReceiveChallengeSmall(...)`
- `ExtensionSumcheckProver::CurrentPolynomial()`
- `ExtensionSumcheckProver::ReceiveChallengeGenericState(...)`

These loops are the most promising high-impact parallel targets.

### 2.7 The large `Other` bucket likely hides an unprofiled serial fold

The strongest current candidate is:

- `Msg0CoeffsAtSuffixChallenges(...)`

This function performs a full extension-domain fold over the coefficient table after the main round loop, and it is not currently broken out as a dedicated profile bucket.

This suggests two likely wins:

- add visibility first
- then either remove duplicated work or parallelize it

## 3. Constraints And Guardrails

Any prover parallelization must preserve the following:

- Transcript schedule remains strictly serial.
- Challenge derivation order remains unchanged.
- Proof bytes remain bit-for-bit identical for the same witness, parameters, and seed.
- `ZZ_p` / `ZZ_pE` / extension contexts are reinitialized correctly inside parallel regions.
- No nested-OpenMP blowup between prover loops and Merkle loops.
- Profiling should continue to report wall-clock style bucket timings, not per-thread accumulated time.

Do not violate these repo-specific facts:

- benchmark conclusions must come from `build-release/`
- `bench_basefold_pcs_prove` and `bench_basefold_pcs_commit` measure different scopes and must not be conflated
- top commit work is already outside the timed prove window, so moving work from timed prove into commit artifacts only counts as a prover win if it matches the intended benchmark contract

## 4. Optimization Workstreams

## Workstream A: Close Profiling Gaps First

Goal:

- shrink the current `Other` bucket into named, attributable buckets before making larger structural changes

Changes:

- add a dedicated profile bucket for `Msg0CoeffsAtSuffixChallenges(...)`
- optionally add a small bucket for the top-oracle extension fallback lift if it ever appears in timed paths
- keep bucket boundaries at coarse wall-clock function level, not inside per-thread loops

Why first:

- current ring profile has over `1 second` in `Other`
- parallelizing the wrong function first is high effort and low confidence

Exit criterion:

- `Other` is reduced to a small residual bucket, or we can name the remaining dominant contributors precisely

## Workstream B: Remove Duplicated Serial Work

Goal:

- remove unnecessary serial work before parallelizing it

Primary target:

- avoid recomputing `msg0_coeffs` from scratch after the round loop if the same suffix-folded state is already available from `ExtensionSumcheckProver`

Candidate direction:

- expose a safe accessor from `ExtensionSumcheckProver` for the final suffix-folded low-dimensional state
- derive `msg0_coeffs` from that state instead of running `Msg0CoeffsAtSuffixChallenges(...)` as a second full fold

Why this matters:

- this likely accounts for a large part of current `Other`
- it reduces both single-thread time and the amount of work later parallel code would need to cover

Risks:

- do not couple proof construction to mutable sumcheck internals in a way that is fragile across paths
- keep exact coefficient ordering and proof serialization unchanged

Exit criterion:

- measurable single-thread win in `bench_basefold_pcs_prove`
- no proof-byte regression

## Workstream C: Parallelize Extension Merkle Build

Current status:

- not started
- deprioritized after the `msg0_coeffs` dedup landed, because current ring profile shows `ExtensionMerkleTree::Build` at about `221 ms`, far below `ExtensionCommitRound`

Goal:

- give the extension prover the same style of Merkle parallelism already available on the base path

Changes:

- refactor `ExtensionMerkleTree::Build(...)` to use the same parallel strategy pattern as `MerkleTree::Build(...)`
- parallelize leaf hashing
- parallelize large internal levels
- reuse or mirror the current Merkle tuning knobs where practical

Expected payoff:

- modest but reliable
- lower risk than sumcheck changes
- especially useful once larger extension oracles appear at bigger `d`

Validation:

- root equality for serial vs threaded build
- proof equality for `OpenMany(...)` over the resulting tree

Exit criterion:

- `ExtensionMerkleTree::Build` bucket scales down with threads
- no digest or proof-format drift

## Workstream D: Make Prover Commit-Round Parallelism Actually Reachable

Current status:

- substantially complete for the current phase
- runtime knobs have landed
- extension default has now settled into `128` plus an automatic ring-only fallback to `64` when no explicit override is present
- thread sweeps showed the existing per-index structure already scales once the threshold is realistic
- remaining work is optional follow-up tuning, not the next blocker

Goal:

- turn the existing element-parallel prover-commit loops into something that activates on realistic workloads

Current blocker:

- hard-coded `kParallelThreshold = 4096`
- effective policy requires at least `8192` work items to get two threads

Changes:

- replace the hard-coded threshold with configurable knobs
- use separate knobs for base and extension commit rounds if needed
- prefer a chunk-size style policy over the current floor division behavior

Suggested knobs:

- `BASEFOLD_PROVER_COMMIT_PARALLEL_THRESHOLD`
- `BASEFOLD_EXT_PROVER_COMMIT_PARALLEL_THRESHOLD`
- optional bench CLI overrides if the repo wants the same style as verifier-query tuning

Design rule:

- preserve serial transcript order
- parallelize only the independent per-index folding inside each round

Expected payoff:

- medium payoff on larger `d`
- small payoff on small `d`
- necessary foundation before judging commit-round scalability fairly

Exit criterion:

- `ExtensionCommitRound` bucket shows thread scaling on `d` where `n_i` is large enough
- no correctness drift

## Workstream E: Parallelize Extension Sumcheck

Current status:

- partially complete
- the hot `ReceiveChallenge*` and `CurrentPolynomial*` loops now parallelize on both the degree-2 small-extension path and the generic extension path
- the remaining visible sumcheck cost is now concentrated mostly in `ReceiveChallenge` / `CurrentPolynomial`; init has been reduced materially but is not zero
- after the new parallelization, `ExtensionSumcheck total` is again the largest bucket in the dominant ring profile, so this workstream remains active rather than closed

Goal:

- attack the largest current prover bucket

Priority order:

1. degree-2 small-extension path
2. generic extension path

Reason for degree-2 first:

- `field-prime128-ext` and `ring-gr-2p16-128-ext` both use degree-2 extension challenges
- these are the highest-value measured contexts

Target loops:

- `CurrentPolynomialSmall(...)`
- `ReceiveChallengeSmall(...)`
- `ExtensionSumcheckProver::CurrentPolynomial()`
- `ExtensionSumcheckProver::ReceiveChallengeGenericState(...)`

Parallelization pattern:

- `CurrentPolynomial*`: thread-local partial accumulators, then deterministic reduction
- `ReceiveChallenge*`: parallel in-place fold over the first half of the table
- keep per-thread scratch local; do not share mutable extension temporaries across threads

Important constraints:

- deterministic output ordering
- safe NTL context initialization in worker threads
- avoid over-counting profile time

Expected payoff:

- high
- this is the largest named bucket in current ring profiles

Exit criterion:

- clear reduction in `ExtensionSumcheck total`
- headline prover time moves meaningfully on ring and field extension-challenge benches

## Workstream F: Optional Query-Proof Construction Parallelism

Goal:

- parallelize smaller remaining independent work only if the earlier workstreams stop moving the total enough

Candidates:

- parallel value-copy loop inside base `MerkleTree::OpenMany(...)`
- parallel value-copy loop inside `ExtensionMerkleTree::OpenMany(...)`
- possibly parallel proof opening across tree levels after query plans are fully derived

Why this is later:

- current profile share is much smaller than extension sumcheck or commit-round folding
- extra complexity is not justified until larger wins are exhausted

Exit criterion:

- only pursue if earlier stages leave a visible post-round query-opening bottleneck

## 5. Recommended Implementation Order

Phase 1:

- add missing profiler bucket(s), especially `Msg0CoeffsAtSuffixChallenges(...)`
- add prover commit-round threshold knobs

Status:

- complete

Phase 2:

- remove duplicated `msg0_coeffs` work if possible
- parallelize `ExtensionMerkleTree::Build(...)`

Status:

- `msg0_coeffs` dedup is complete
- `ExtensionMerkleTree::Build(...)` parallelization is intentionally deferred

Phase 3:

- parallelize degree-2 `ExtensionSumcheckProver`
- parallelize extension commit-round loops with better activation thresholds

Revised status / order:

- split this phase
- extension commit-round scaling/tuning is complete enough for now
- the first round of `ExtensionSumcheckProver` parallelization is now in, including init-path parallelization
- after refreshing the ring profile, the next non-sumcheck target is `ExtensionCommitRound`, not `ExtensionMerkleTree::Build(...)`

Phase 4:

- generalize sumcheck parallelization to generic extension degree
- tune threshold defaults using real thread sweeps

Revised status:

- generic extension `ReceiveChallenge` / `CurrentPolynomial` parallelization is already in
- threshold refinement is now complete enough for the current phase
- the remaining value here is mainly any follow-up on the sumcheck init path plus future degree-generic cleanup if new contexts need it

Phase 5:

- only if still needed, optimize query-proof construction and other residual buckets

## 6. Validation Matrix

## Functional Validation

- run existing PCS tests in `build-release/`
- add targeted regression tests that compare serial vs threaded proof bytes for the same seed
- compare verifier acceptance for serial-produced and threaded-produced proofs
- cover both:
  - degree-2 extension contexts
  - at least one generic extension context

Minimum contexts:

- `field-prime128-ext`
- `ring-gr-2p16-128-ext`
- one degree-3 extension context such as `field-prime64-ext` or `ring-gr-2p16-64-ext`

## Performance Validation

Use two benchmark modes:

- profiling mode for bucket attribution:
  - `bench_basefold_pcs_prove --profile --warmup 0 --reps 1`
- timing mode for headline speedup:
  - `bench_basefold_pcs_prove --warmup 1 --reps 3`

Thread sweep:

- `1, 2, 4, 8, 16, 32`

Primary benchmark contexts:

- `field-prime128-ext`, `d=12`, `queries=2048`
- `ring-gr-2p16-128-ext`, `d=12`, `queries=2048`

Additional large-`d` checks:

- one larger `d` point where prover commit rounds are definitely above the parallel threshold

Success criteria:

- single-thread prover time does not regress materially
- `8` threads gives visible prover improvement in at least the ring degree-2 extension context
- profile buckets confirm the win comes from named hotspots, not measurement noise

## 7. Stop Rules

Stop and reassess before continuing if any of the following happens:

- proof bytes change unexpectedly
- verifier accepts only one of serial/threaded proofs
- `Other` remains dominant even after instrumentation
- sumcheck parallelization increases total prover time because of NTL or thread overhead
- nested parallel regions create oversubscription that hides local wins

If a phase fails, revert to the previous stable phase and keep the measurement evidence.

## 8. Current Next Patch Recommendation

The first patch recommendation in the original draft has already been overtaken by landed work.

The next patch should still focus on `ExtensionCommitRound`, and no longer on threshold knobs.

Recommended order:

1. Keep the new `ext_elements_per_thread=128` default; treat the thread/threshold sweep as complete enough for now.
2. If commit-round optimization remains the goal, improve the remaining inner structure or data layout of the degree-2 extension kernel rather than tuning the knob further.
3. Only after commit-round no longer dominates the ring degree-2 profile, move focus to the next-largest residual bucket.

Reason:

- latest 8-thread ring profile with the new default is now roughly:
  - `ExtensionCommitRound`: about `459 ms`
  - `ExtensionSumcheck total`: about `342 ms`
  - `ExtensionMerkleTree::Build`: about `208 ms`
- so the commit-round knob work has paid off, but commit-round still edges out the other named buckets on the dominant ring workload
