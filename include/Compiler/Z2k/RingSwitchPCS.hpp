#ifndef BASEFOLD_Z2K_RING_SWITCH_PCS_HPP_
#define BASEFOLD_Z2K_RING_SWITCH_PCS_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

#include "Compiler/Z2k/PCSBackend.hpp"
#include "PCS/Common/Sumcheck.hpp"

namespace basefold {

enum class RingSwitchBasisKind {
  kPolynomial = 0,
};

struct RingSwitchBasisDescriptor {
  RingSwitchBasisKind kind = RingSwitchBasisKind::kPolynomial;
  long dimension = 0;
};

struct RingSwitchPCSSetupInput {
  long ell = 0;
  long kappa = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  RingSwitchBasisDescriptor alpha_basis;
  RingSwitchBasisDescriptor beta_basis;
  Z2kPCSBackendHandle backend;
};

struct RingSwitchPCSParams {
  long ell = 0;
  long kappa = 0;
  long ell_prime = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  RingSwitchBasisDescriptor alpha_basis;
  RingSwitchBasisDescriptor beta_basis;
  Z2kPCSBackendHandle backend;
};

struct RingSwitchComponentTensor {
  long basis_dimension = 0;
  long ell_prime = 0;
  // Row-major storage of A_{u||w}: row `u`, then column `w`.
  NTL::vec_ZZ_pE a_by_u_then_w;
  // Boolean-hypercube value table of the verifier-known polynomial r, with
  // linear index index(u || w) = u + 2^kappa * w.
  NTL::vec_ZZ_pE r_table;
  // Monomial-basis coefficients of r, cached for reuse with the current
  // backend core and existing multilinear helpers.
  NTL::vec_ZZ_pE r_monomial_coeffs;
};

struct RingSwitchPCSCommitArtifacts {
  NTL::vec_ZZ_pE t_packed_table;
  NTL::vec_ZZ_pE t_packed_monomial_coeffs;
  Z2kPCSBackendCommitArtifacts backend_commit_artifacts;
  MerkleRoot commitment;
};

struct RingSwitchPCSOuterCommitArtifacts {
  NTL::vec_ZZ_pE t_packed_table;
  NTL::vec_ZZ_pE t_packed_monomial_coeffs;
};

struct RingSwitchPCSOuterEvalProof {
  std::vector<NTL::ZZ_pE> s_by_u;
  std::vector<QuadraticPoly> h_by_level;
  NTL::ZZ_pE t_star;
};

struct RingSwitchPCSEvalProof {
  std::vector<NTL::ZZ_pE> s_by_u;
  std::vector<QuadraticPoly> h_by_level;
  NTL::ZZ_pE t_star;
  Z2kPCSBackendEvalProof backend_proof;
};

RingSwitchBasisDescriptor ActivePolynomialBasisDescriptor();

void ValidateCurrentZ2kRingContextOrThrow(const NTL::ZZ &base_modulus,
                                          const NTL::ZZ_pX &extension_modulus,
                                          long kappa);

void ValidateRingSwitchPCSParamsOrThrow(const RingSwitchPCSParams &params);

RingSwitchPCSParams RingSwitchPCSSetup(const RingSwitchPCSSetupInput &input);

// Converts a length-2^d Boolean-hypercube value table (Lagrange-basis
// coefficients over {0,1}^d) into the corresponding monomial-basis
// coefficients.
NTL::vec_ZZ_pE BooleanHypercubeTableToMonomialCoeffs(
    const NTL::vec_ZZ_pE &table_values);

// Packs the paper-style Boolean-hypercube table of t over Z_{2^k} into the
// paper-style Boolean-hypercube table of t' over GR(2^k, 2^kappa).
NTL::vec_ZZ_pE PackZ2kCoeffsToGREvals(const RingSwitchPCSParams &params,
                                      const NTL::vec_ZZ_pE &t_table);

MerkleRoot RingSwitchPCSCommit(const RingSwitchPCSParams &params,
                               const NTL::vec_ZZ_pE &t_table);

RingSwitchPCSOuterCommitArtifacts RingSwitchPCSBuildOuterCommitArtifacts(
    const RingSwitchPCSParams &params, const NTL::vec_ZZ_pE &t_table);

RingSwitchPCSCommitArtifacts RingSwitchPCSBuildCommitArtifacts(
    const RingSwitchPCSParams &params, const NTL::vec_ZZ_pE &t_table);

RingSwitchPCSOuterEvalProof RingSwitchPCSProveOuterEval(
    const RingSwitchPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries);

RingSwitchPCSOuterEvalProof RingSwitchPCSProveOuterEvalFromCommitArtifacts(
    const RingSwitchPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries,
    const RingSwitchPCSOuterCommitArtifacts &commit_artifacts);

RingSwitchPCSEvalProof RingSwitchPCSProveEval(
    const RingSwitchPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries);

RingSwitchPCSEvalProof RingSwitchPCSProveEvalFromCommitArtifacts(
    const RingSwitchPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const RingSwitchPCSCommitArtifacts &commit_artifacts);

bool RingSwitchPCSVerifyEval(const RingSwitchPCSParams &params,
                             const MerkleRoot &commitment,
                             const std::vector<FieldElement> &z,
                             const FieldElement &claimed_s, long num_queries,
                             const RingSwitchPCSEvalProof &proof);

bool RingSwitchPCSVerifyOuterEval(const RingSwitchPCSParams &params,
                                  const MerkleRoot &commitment,
                                  const std::vector<FieldElement> &z,
                                  const FieldElement &claimed_s,
                                  long num_queries,
                                  const RingSwitchPCSOuterEvalProof &proof);

// Prover-local state produced by Commit and consumed by Prove.
struct RingSwitchPCSCommittedWitness {
  MerkleRoot commitment;
  NTL::vec_ZZ_pE t_table;
  RingSwitchPCSCommitArtifacts commit_artifacts;
};

struct RingSwitchPCSOuterCommittedWitness {
  MerkleRoot commitment;
  NTL::vec_ZZ_pE t_table;
  RingSwitchPCSOuterCommitArtifacts commit_artifacts;
};

class RingSwitchPCSProver {
 public:
  explicit RingSwitchPCSProver(const RingSwitchPCSParams &params)
      : params_(params) {
    ValidateRingSwitchPCSParamsOrThrow(params_);
  }

  explicit RingSwitchPCSProver(const RingSwitchPCSSetupInput &input)
      : params_(RingSwitchPCSSetup(input)) {}

  RingSwitchPCSCommittedWitness Commit(const NTL::vec_ZZ_pE &t_table) const {
    RingSwitchPCSCommittedWitness committed;
    committed.t_table = t_table;
    committed.commit_artifacts = RingSwitchPCSBuildCommitArtifacts(params_, t_table);
    committed.commitment = committed.commit_artifacts.commitment;
    return committed;
  }

  RingSwitchPCSEvalProof Prove(const RingSwitchPCSCommittedWitness &committed,
                               const std::vector<FieldElement> &z,
                               const FieldElement &claimed_s,
                               long num_queries) const {
    if (committed.commitment != committed.commit_artifacts.commitment) {
      NTL::LogicError(
          "RingSwitchPCSProver::Prove: committed witness has inconsistent commitment");
    }
    return RingSwitchPCSProveEvalFromCommitArtifacts(
        params_, committed.t_table, z, claimed_s, num_queries,
        committed.commit_artifacts);
  }

  const RingSwitchPCSParams &params() const { return params_; }

 private:
  RingSwitchPCSParams params_;
};

class RingSwitchPCSOuterProver {
 public:
  explicit RingSwitchPCSOuterProver(const RingSwitchPCSParams &params)
      : params_(params) {
    ValidateRingSwitchPCSParamsOrThrow(params_);
  }

  explicit RingSwitchPCSOuterProver(const RingSwitchPCSSetupInput &input)
      : params_(RingSwitchPCSSetup(input)) {}

  RingSwitchPCSOuterCommittedWitness Commit(
      const NTL::vec_ZZ_pE &t_table, const MerkleRoot &commitment) const {
    RingSwitchPCSOuterCommittedWitness committed;
    committed.commitment = commitment;
    committed.t_table = t_table;
    committed.commit_artifacts =
        RingSwitchPCSBuildOuterCommitArtifacts(params_, t_table);
    return committed;
  }

  RingSwitchPCSOuterEvalProof Prove(
      const RingSwitchPCSOuterCommittedWitness &committed,
      const std::vector<FieldElement> &z, const FieldElement &claimed_s,
      long num_queries) const {
    return RingSwitchPCSProveOuterEvalFromCommitArtifacts(
        params_, committed.t_table, committed.commitment, z, claimed_s,
        num_queries, committed.commit_artifacts);
  }

  const RingSwitchPCSParams &params() const { return params_; }

 private:
  RingSwitchPCSParams params_;
};

class RingSwitchPCSVerifier {
 public:
  explicit RingSwitchPCSVerifier(const RingSwitchPCSParams &params)
      : params_(params) {
    ValidateRingSwitchPCSParamsOrThrow(params_);
  }

  explicit RingSwitchPCSVerifier(const RingSwitchPCSSetupInput &input)
      : params_(RingSwitchPCSSetup(input)) {}

  bool Verify(const MerkleRoot &commitment,
              const std::vector<FieldElement> &z,
              const FieldElement &claimed_s, long num_queries,
              const RingSwitchPCSEvalProof &proof) const {
    return RingSwitchPCSVerifyEval(params_, commitment, z, claimed_s,
                                   num_queries, proof);
  }

  const RingSwitchPCSParams &params() const { return params_; }

 private:
  RingSwitchPCSParams params_;
};

class RingSwitchPCSOuterVerifier {
 public:
  explicit RingSwitchPCSOuterVerifier(const RingSwitchPCSParams &params)
      : params_(params) {
    ValidateRingSwitchPCSParamsOrThrow(params_);
  }

  explicit RingSwitchPCSOuterVerifier(const RingSwitchPCSSetupInput &input)
      : params_(RingSwitchPCSSetup(input)) {}

  bool Verify(const MerkleRoot &commitment,
              const std::vector<FieldElement> &z,
              const FieldElement &claimed_s, long num_queries,
              const RingSwitchPCSOuterEvalProof &proof) const {
    return RingSwitchPCSVerifyOuterEval(params_, commitment, z, claimed_s,
                                        num_queries, proof);
  }

  const RingSwitchPCSParams &params() const { return params_; }

 private:
  RingSwitchPCSParams params_;
};

struct RingSwitchPCSSetupOutput {
  RingSwitchPCSParams params;
  RingSwitchPCSProver prover;
  RingSwitchPCSVerifier verifier;
};

inline RingSwitchPCSSetupOutput RingSwitchPCSSetupProtocol(
    const RingSwitchPCSSetupInput &input) {
  const RingSwitchPCSParams params = RingSwitchPCSSetup(input);
  return {params, RingSwitchPCSProver(params), RingSwitchPCSVerifier(params)};
}

struct RingSwitchPCSOuterSetupOutput {
  RingSwitchPCSParams params;
  RingSwitchPCSOuterProver prover;
  RingSwitchPCSOuterVerifier verifier;
};

inline RingSwitchPCSOuterSetupOutput RingSwitchPCSSetupOuterProtocol(
    const RingSwitchPCSSetupInput &input) {
  const RingSwitchPCSParams params = RingSwitchPCSSetup(input);
  return {params, RingSwitchPCSOuterProver(params),
          RingSwitchPCSOuterVerifier(params)};
}

NTL::vec_ZZ_pE DecomposeGRElementToBaseCoeffsPolynomialBasis(
    const RingSwitchPCSParams &params, const NTL::ZZ_pE &element);

NTL::vec_ZZ_pE DecomposeGRElementToBaseCoeffs(
    const RingSwitchPCSParams &params, const NTL::ZZ_pE &element,
    const RingSwitchBasisDescriptor &basis);

RingSwitchComponentTensor BuildRingSwitchComponentTensor(
    const RingSwitchPCSParams &params,
    const std::vector<NTL::ZZ_pE> &r_suffix);

}  // namespace basefold

#endif  // BASEFOLD_Z2K_RING_SWITCH_PCS_HPP_
