#ifndef BASEFOLD_IOPP_HPP_
#define BASEFOLD_IOPP_HPP_

#include <NTL/ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

#include "PCS/BaseFold/FoldableCode.hpp"
#include "PCS/Common/Merkle.hpp"
#include "PCS/Common/Transcript.hpp"

namespace basefold {

// This header declares a minimal set of interfaces for implementing the
// BaseFold IOPP (Protocols 2 and 3 in:
//   Zeilberger et al., 2023, "BaseFold: Efficient Field-Agnostic Polynomial
//   Commitment Schemes from Foldable Codes")
// specialized to the ζ-foldable recursion already encoded by FoldableCodeParams.
//
// -----------------------------------------------------------------------------
// Indexing conventions (paper vs. code)
// -----------------------------------------------------------------------------
// The paper uses 1-based indices for oracle access. In this code we use 0-based:
// - An oracle π_i is represented as a dense vector of length n_i, indexed by
//   j ∈ [0, n_i).
// - For a folding round i (0 <= i < d), π_{i+1} has length n_{i+1} = 2*n_i and
//   the consistency check reads π_{i+1}[µ] and π_{i+1}[µ + n_i], for µ ∈ [0, n_i).
//
// -----------------------------------------------------------------------------
// ζ-specialization (FoldableCodeParams)
// -----------------------------------------------------------------------------
// The generic BaseFold IOPP assumes diagonal matrices T_i and T'_i with
// diag(T_i)[j] != diag(T'_i)[j]. Your FoldableCodeParams stores diag_T[i] =
// diag(T_i) and a fixed ζ such that:
//   diag(T'_i)[j] := ζ * diag(T_i)[j].
// Under FoldableCodeParams' validation, ζ != 0,1 and diag(T_i)[j] is non-zero
// (and intended to be a unit), hence the two folding points are distinct and the
// denominator (diag(T'_i)[j] - diag(T_i)[j]) is a unit.
//
// This is exactly the additional condition needed to run the same BaseFold IOPP
// over a Galois ring: the verifier/prover only ever "divide" by these unit
// denominators (see Basefold_over_GR.pdf, Protocol 2).

using FieldElement = NTL::ZZ_pE;
using Oracle = NTL::vec_ZZ_pE;

// Returns k_i = k0 * 2^i (message length of C_i).
long MessageLengthAtLevel(const FoldableCodeParams &params, long level);

// Returns n_i = c * k0 * 2^i (codeword length of C_i).
long CodewordLengthAtLevel(const FoldableCodeParams &params, long level);

// For folding round i and coordinate j, returns the two interpolation x-points:
//   x_left  = diag(T_i)[j]
//   x_right = diag(T'_i)[j] = ζ * diag(T_i)[j]
// Preconditions:
// - 0 <= level_i < params.d
// - 0 <= j < n_i
void FoldingPoints(FieldElement &x_left, FieldElement &x_right,
                   const FoldableCodeParams &params, long level_i, long j);

// Evaluates the unique degree-1 polynomial p(X) such that:
//   p(x1) = y1, p(x2) = y2
// at the point x.
// Precondition: x1 != x2 (and x2-x1 is invertible in the ambient field/ring).
FieldElement EvalLineAt(const FieldElement &x, const FieldElement &x1,
                        const FieldElement &y1, const FieldElement &x2,
                        const FieldElement &y2);

// Verifier challenges α_i for i=0..d-1 used to fold π_{i+1} -> π_i.
// In Protocol 2, the verifier samples α_{d-1},...,α_0 in that order, but it is
// convenient to store them by index: alphas[i] == α_i.
struct IOPPChallenges {
  std::vector<FieldElement> alphas;  // size == params.d
};

// A collection of oracles π_0..π_d.
// Convention: pi[i] stores π_i, so pi.size() == params.d + 1 and pi[d] is the
// input oracle π_d.
struct IOPPOracles {
  std::vector<Oracle> oracles_by_level;
};

// -----------------------------------------------------------------------------
// Prover side (Protocol 2: IOPP.commit)
// -----------------------------------------------------------------------------

// Computes one folding round:
//   input:  π_{i+1} (length 2*n_i), challenge α_i
//   output: π_i     (length n_i)
// with coordinate-wise degree-1 interpolation at points (diag(T_i)[j], ...) and
// (diag(T'_i)[j], ...).
void ProverCommitRound(Oracle &pi_i, const Oracle &pi_ip1,
                       const FieldElement &alpha_i, long level_i,
                       const FoldableCodeParams &params);

// Computes the full commit phase for a fixed π_d and challenges α_0..α_{d-1}.
// Expected output layout:
// - oracles.oracles_by_level has size params.d + 1
// - oracles.oracles_by_level[d] == pi_d
// - oracles.oracles_by_level[i] is the prover's derived oracle π_i for all i<d.
void ProverCommitAll(IOPPOracles &oracles, const Oracle &pi_d,
                     const IOPPChallenges &challenges,
                     const FoldableCodeParams &params);

// -----------------------------------------------------------------------------
// Verifier side (Protocol 3: IOPP.query)
// -----------------------------------------------------------------------------

// Query-plan for one repetition of IOPP.query.
struct IOPPQueryPlan {
  // Initial µ sampled uniformly from [0, n_{d-1}).
  long initial_mu = 0;

  // µ used at each folding level i (0 <= i < d) when checking π_{i+1} -> π_i.
  // Size must equal params.d.
  //
  // In particular, the verifier reads:
  //   π_{i+1}[µ_i], π_{i+1}[µ_i + n_i], and π_i[µ_i]
  // where n_i = CodewordLengthAtLevel(params, i).
  std::vector<long> mu_by_level;
};

// Computes the µ update schedule from Protocol 3:
// - Sample µ_{d-1} ∈ [0, n_{d-1})
// - For i = d-1 down to 0:
//     use µ_i to check π_{i+1} -> π_i
//     if i>0 and µ_i >= n_{i-1}, set µ_{i-1} = µ_i - n_{i-1}
IOPPQueryPlan MakeQueryPlan(long initial_mu, const FoldableCodeParams &params);

// The explicit openings needed to verify one IOPP.query repetition without
// random-access to the oracles (useful when you later Merkle-commit to oracles).
struct IOPPQueryOpenings {
  // For each i in [0, d):
  //   upper_left_by_level[i]  = π_{i+1}[µ_i]
  //   upper_right_by_level[i] = π_{i+1}[µ_i + n_i]
  //   folded_by_level[i]      = π_i[µ_i]
  std::vector<FieldElement> upper_left_by_level;
  std::vector<FieldElement> upper_right_by_level;
  std::vector<FieldElement> folded_by_level;

  // Final oracle π_0 (typically small enough that the verifier reads all entries
  // to check π_0 ∈ C_0).
  Oracle pi0_codeword;
};

// Verifies a single IOPP.query repetition from explicit openings.
// Returns true iff:
// - All folding consistency checks pass (Protocol 3, steps 1-4), and
// - π_0 is a valid codeword w.r.t. generator matrix G0 (Protocol 3 final step).
bool VerifyQueryFromOpenings(const IOPPQueryPlan &plan,
                             const IOPPChallenges &challenges,
                             const IOPPQueryOpenings &openings,
                             const FoldableCodeParams &params);

// Convenience verifier for the oracle model: directly indexes into concrete
// π_0..π_d vectors according to the query plan.
bool VerifyQueryFromOracles(const IOPPQueryPlan &plan,
                            const IOPPChallenges &challenges,
                            const IOPPOracles &oracles,
                            const FoldableCodeParams &params);

// Merkle roots for π_0..π_d (roots_by_level[i] commits to π_i).
struct IOPPMerkleCommitments {
  std::vector<MerkleRoot> roots_by_level;  // size == params.d + 1
};

// Derives IOPP challenges α_0..α_{d-1} using Fiat-Shamir and the Merkle roots.
//
// Recommended transcript schedule (mirrors Protocol 2 adaptivity):
//   absorb(root_d)
//   for i = d-1 downto 0:
//     alpha_i = challenge_field("alpha/<i>")
//     absorb(root_i)
//
// Preconditions:
// - commitments.roots_by_level has size params.d + 1.
IOPPChallenges FiatShamirDeriveChallenges(
    FiatShamirTranscript &transcript,
    const IOPPMerkleCommitments &commitments,
    const FoldableCodeParams &params);

// Derives `num_queries` independent query plans from the transcript.
// Expected schedule:
// - Call after FiatShamirDeriveChallenges (so the transcript has absorbed all
//   oracle roots in the intended order).
// - Each query draws µ ∈ [0, n_{d-1}) via ChallengeIndex and expands it via
//   MakeQueryPlan.
std::vector<IOPPQueryPlan> FiatShamirDeriveQueryPlans(
    FiatShamirTranscript &transcript, long num_queries,
    const FoldableCodeParams &params);

// -----------------------------------------------------------------------------
// C0 codeword check (Protocol 3 final step)
// -----------------------------------------------------------------------------

// Returns true iff pi0 is in the linear code C0 generated by params.G0.
//
// For a finite-field implementation, a typical approach is to solve for m0:
//   m0 * G0 == pi0
// and accept iff a solution exists.
bool IsCodewordC0(const Oracle &pi0, const FoldableCodeParams &params);

// Optional helper: attempts to recover some m0 such that m0 * G0 == pi0.
// Returns false iff pi0 is not in C0.
bool DecodeC0(NTL::vec_ZZ_pE &msg0_out, const Oracle &pi0,
              const FoldableCodeParams &params);

}  // namespace basefold

#endif  // BASEFOLD_IOPP_HPP_
