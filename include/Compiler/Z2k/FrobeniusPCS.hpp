#ifndef BASEFOLD_Z2K_FROBENIUS_PCS_HPP_
#define BASEFOLD_Z2K_FROBENIUS_PCS_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

#include "Compiler/Z2k/PCSBackend.hpp"
#include "GaloisRing/FrobeniusBasis.hpp"
#include "PCS/Common/Sumcheck.hpp"

namespace basefold {

// Optional paper-style setup input for Protocol 2. When enabled through
// FrobeniusPCSSetupInput::use_provided_basis, setup validates that normal_basis
// is a genuine Frobenius-ordered normal basis with the supplied dual basis.
// If has_teichmuller_generator is false, setup derives a generator internally
// and still performs the same strong validation.
struct FrobeniusPCSProvidedBasisInput {
  NormalBasisData normal_basis;
  bool has_teichmuller_generator = false;
  NTL::ZZ_pE teichmuller_generator;
};

struct FrobeniusPCSSetupInput {
  long ell = 0;
  long kappa = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  // false: discover the normal/dual basis during setup.
  // true: use provided_basis instead and validate it against the active context.
  bool use_provided_basis = false;
  FrobeniusPCSProvidedBasisInput provided_basis;
  long teichmuller_generator_max_trials = 1024;
  long affine_search_limit = 4;
  Z2kPCSBackendHandle backend;
};

struct FrobeniusPCSPrecomputedTables {
  // tau_basis_rows[i][j] = tau^i(beta_j), sigma_basis_rows[i][j] = sigma^i(beta_j)
  // for the normal basis beta.
  std::vector<std::vector<NTL::ZZ_pE>> tau_basis_rows;
  std::vector<std::vector<NTL::ZZ_pE>> sigma_basis_rows;
  // tau_alpha_by_u_then_i[u][i] = tau^i(alpha_u) for the dual basis alpha.
  std::vector<std::vector<NTL::ZZ_pE>> tau_alpha_by_u_then_i;
};

struct FrobeniusPCSParams {
  long ell = 0;
  long kappa = 0;
  long ell_prime = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  FrobeniusBasisData basis_data;
  FrobeniusPCSPrecomputedTables precomputed;
  Z2kPCSBackendHandle backend;
};

struct FrobeniusPCSOuterCommitArtifacts {
  NTL::vec_ZZ_pE t_packed_table;
  NTL::vec_ZZ_pE t_packed_monomial_coeffs;
};

struct FrobeniusPCSCommitArtifacts {
  NTL::vec_ZZ_pE t_packed_table;
  NTL::vec_ZZ_pE t_packed_monomial_coeffs;
  Z2kPCSBackendCommitArtifacts backend_commit_artifacts;
  MerkleRoot commitment;
};

struct FrobeniusPCSOuterEvalProof {
  std::vector<NTL::ZZ_pE> s_by_i;
  std::vector<QuadraticPoly> h_by_level;
  NTL::ZZ_pE t_star;
};

struct FrobeniusPCSEvalProof {
  std::vector<NTL::ZZ_pE> s_by_i;
  std::vector<QuadraticPoly> h_by_level;
  NTL::ZZ_pE t_star;
  Z2kPCSBackendEvalProof backend_proof;
};

void ValidateCurrentZ2kFrobeniusContextOrThrow(
    const NTL::ZZ &base_modulus, const NTL::ZZ_pX &extension_modulus,
    long kappa);

void ValidateFrobeniusPCSParamsOrThrow(const FrobeniusPCSParams &params);

FrobeniusPCSParams FrobeniusPCSSetup(const FrobeniusPCSSetupInput &input);

NTL::vec_ZZ_pE PackZ2kTableToFrobeniusGREvals(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table);

NTL::vec_ZZ_pE DecomposeGRElementToBaseCoeffsFrobeniusBasis(
    const FrobeniusPCSParams &params, const NTL::ZZ_pE &element);

MerkleRoot FrobeniusPCSCommit(const FrobeniusPCSParams &params,
                              const NTL::vec_ZZ_pE &t_table);

FrobeniusPCSOuterCommitArtifacts FrobeniusPCSBuildOuterCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table);

FrobeniusPCSCommitArtifacts FrobeniusPCSBuildCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table);

FrobeniusPCSOuterEvalProof FrobeniusPCSProveOuterEval(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries);

// Unchecked variant intended for benchmarking / hot paths.
//
// Differences vs FrobeniusPCSProveOuterEval:
// - Does NOT validate prove inputs.
// - Does NOT validate the claimed evaluation by recomputing t(z).
// - Does NOT run prover-side honest-witness reconstruction / equality checks.
FrobeniusPCSOuterEvalProof FrobeniusPCSProveOuterEvalUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries);

FrobeniusPCSOuterEvalProof FrobeniusPCSProveOuterEvalFromCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries,
    const FrobeniusPCSOuterCommitArtifacts &commit_artifacts);

// Unchecked variant intended for benchmarking / hot paths.
//
// Differences vs FrobeniusPCSProveOuterEvalFromCommitArtifacts:
// - Does NOT validate params or commit_artifacts shape.
// - Does NOT validate prove inputs.
// - Does NOT validate the claimed evaluation by recomputing t(z).
// - Does NOT run prover-side honest-witness reconstruction / equality checks.
FrobeniusPCSOuterEvalProof FrobeniusPCSProveOuterEvalFromCommitArtifactsUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries,
    const FrobeniusPCSOuterCommitArtifacts &commit_artifacts);

FrobeniusPCSEvalProof FrobeniusPCSProveEval(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries);

// Unchecked variant intended for benchmarking / hot paths.
//
// Differences vs FrobeniusPCSProveEval:
// - Does NOT validate prove inputs.
// - Does NOT validate the claimed evaluation by recomputing t(z).
// - Does NOT run prover-side honest-witness reconstruction / equality checks
//   in the outer Frobenius layer.
FrobeniusPCSEvalProof FrobeniusPCSProveEvalUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries);

FrobeniusPCSEvalProof FrobeniusPCSProveEvalFromCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const FrobeniusPCSCommitArtifacts &commit_artifacts);

// Unchecked variant intended for benchmarking / hot paths.
//
// Differences vs FrobeniusPCSProveEvalFromCommitArtifacts:
// - Does NOT validate params or commit_artifacts shape in the Frobenius layer.
// - Does NOT validate prove inputs.
// - Does NOT validate the claimed evaluation by recomputing t(z).
// - Does NOT run prover-side honest-witness reconstruction / equality checks
//   in the outer Frobenius layer.
FrobeniusPCSEvalProof FrobeniusPCSProveEvalFromCommitArtifactsUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const FrobeniusPCSCommitArtifacts &commit_artifacts);

bool FrobeniusPCSVerifyEval(const FrobeniusPCSParams &params,
                            const MerkleRoot &commitment,
                            const std::vector<FieldElement> &z,
                            const FieldElement &claimed_s, long num_queries,
                            const FrobeniusPCSEvalProof &proof);

bool FrobeniusPCSVerifyOuterEval(const FrobeniusPCSParams &params,
                                 const MerkleRoot &commitment,
                                 const std::vector<FieldElement> &z,
                                 const FieldElement &claimed_s,
                                 long num_queries,
                                 const FrobeniusPCSOuterEvalProof &proof);

struct FrobeniusPCSCommittedWitness {
  MerkleRoot commitment;
  NTL::vec_ZZ_pE t_table;
  FrobeniusPCSCommitArtifacts commit_artifacts;
};

struct FrobeniusPCSOuterCommittedWitness {
  MerkleRoot commitment;
  NTL::vec_ZZ_pE t_table;
  FrobeniusPCSOuterCommitArtifacts commit_artifacts;
};

class FrobeniusPCSProver {
 public:
  explicit FrobeniusPCSProver(const FrobeniusPCSParams &params)
      : params_(params) {
    ValidateFrobeniusPCSParamsOrThrow(params_);
  }

  explicit FrobeniusPCSProver(const FrobeniusPCSSetupInput &input)
      : params_(FrobeniusPCSSetup(input)) {}

  FrobeniusPCSCommittedWitness Commit(const NTL::vec_ZZ_pE &t_table) const {
    FrobeniusPCSCommittedWitness committed;
    committed.t_table = t_table;
    committed.commit_artifacts = FrobeniusPCSBuildCommitArtifacts(params_, t_table);
    committed.commitment = committed.commit_artifacts.commitment;
    return committed;
  }

  FrobeniusPCSEvalProof Prove(const FrobeniusPCSCommittedWitness &committed,
                              const std::vector<FieldElement> &z,
                              const FieldElement &claimed_s,
                              long num_queries) const {
    if (committed.commitment != committed.commit_artifacts.commitment) {
      NTL::LogicError(
          "FrobeniusPCSProver::Prove: committed witness has inconsistent commitment");
    }
    return FrobeniusPCSProveEvalFromCommitArtifacts(
        params_, committed.t_table, z, claimed_s, num_queries,
        committed.commit_artifacts);
  }

  const FrobeniusPCSParams &params() const { return params_; }

 private:
  FrobeniusPCSParams params_;
};

class FrobeniusPCSOuterProver {
 public:
  explicit FrobeniusPCSOuterProver(const FrobeniusPCSParams &params)
      : params_(params) {
    ValidateFrobeniusPCSParamsOrThrow(params_);
  }

  explicit FrobeniusPCSOuterProver(const FrobeniusPCSSetupInput &input)
      : params_(FrobeniusPCSSetup(input)) {}

  FrobeniusPCSOuterCommittedWitness Commit(
      const NTL::vec_ZZ_pE &t_table, const MerkleRoot &commitment) const {
    FrobeniusPCSOuterCommittedWitness committed;
    committed.commitment = commitment;
    committed.t_table = t_table;
    committed.commit_artifacts =
        FrobeniusPCSBuildOuterCommitArtifacts(params_, t_table);
    return committed;
  }

  FrobeniusPCSOuterEvalProof Prove(
      const FrobeniusPCSOuterCommittedWitness &committed,
      const std::vector<FieldElement> &z, const FieldElement &claimed_s,
      long num_queries) const {
    return FrobeniusPCSProveOuterEvalFromCommitArtifacts(
        params_, committed.t_table, committed.commitment, z, claimed_s,
        num_queries, committed.commit_artifacts);
  }

  const FrobeniusPCSParams &params() const { return params_; }

 private:
  FrobeniusPCSParams params_;
};

class FrobeniusPCSVerifier {
 public:
  explicit FrobeniusPCSVerifier(const FrobeniusPCSParams &params)
      : params_(params) {
    ValidateFrobeniusPCSParamsOrThrow(params_);
  }

  explicit FrobeniusPCSVerifier(const FrobeniusPCSSetupInput &input)
      : params_(FrobeniusPCSSetup(input)) {}

  bool Verify(const MerkleRoot &commitment,
              const std::vector<FieldElement> &z,
              const FieldElement &claimed_s, long num_queries,
              const FrobeniusPCSEvalProof &proof) const {
    return FrobeniusPCSVerifyEval(params_, commitment, z, claimed_s,
                                  num_queries, proof);
  }

  const FrobeniusPCSParams &params() const { return params_; }

 private:
  FrobeniusPCSParams params_;
};

class FrobeniusPCSOuterVerifier {
 public:
  explicit FrobeniusPCSOuterVerifier(const FrobeniusPCSParams &params)
      : params_(params) {
    ValidateFrobeniusPCSParamsOrThrow(params_);
  }

  explicit FrobeniusPCSOuterVerifier(const FrobeniusPCSSetupInput &input)
      : params_(FrobeniusPCSSetup(input)) {}

  bool Verify(const MerkleRoot &commitment,
              const std::vector<FieldElement> &z,
              const FieldElement &claimed_s, long num_queries,
              const FrobeniusPCSOuterEvalProof &proof) const {
    return FrobeniusPCSVerifyOuterEval(params_, commitment, z, claimed_s,
                                       num_queries, proof);
  }

  const FrobeniusPCSParams &params() const { return params_; }

 private:
  FrobeniusPCSParams params_;
};

struct FrobeniusPCSSetupOutput {
  FrobeniusPCSParams params;
  FrobeniusPCSProver prover;
  FrobeniusPCSVerifier verifier;
};

inline FrobeniusPCSSetupOutput FrobeniusPCSSetupProtocol(
    const FrobeniusPCSSetupInput &input) {
  const FrobeniusPCSParams params = FrobeniusPCSSetup(input);
  return {params, FrobeniusPCSProver(params), FrobeniusPCSVerifier(params)};
}

struct FrobeniusPCSOuterSetupOutput {
  FrobeniusPCSParams params;
  FrobeniusPCSOuterProver prover;
  FrobeniusPCSOuterVerifier verifier;
};

inline FrobeniusPCSOuterSetupOutput FrobeniusPCSSetupOuterProtocol(
    const FrobeniusPCSSetupInput &input) {
  const FrobeniusPCSParams params = FrobeniusPCSSetup(input);
  return {params, FrobeniusPCSOuterProver(params),
          FrobeniusPCSOuterVerifier(params)};
}

}  // namespace basefold

#endif  // BASEFOLD_Z2K_FROBENIUS_PCS_HPP_
