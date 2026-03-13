#ifndef BASEFOLD_Z2K_FROBENIUS_PCS_HPP_
#define BASEFOLD_Z2K_FROBENIUS_PCS_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_pE.h>

#include "Compiler/Z2k/PCSBackend.hpp"
#include "GaloisRing/FrobeniusBasis.hpp"

namespace basefold {

struct FrobeniusPCSSetupInput {
  long ell = 0;
  long kappa = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  long teichmuller_generator_max_trials = 1024;
  long affine_search_limit = 4;
  Z2kPCSBackendHandle backend;
};

struct FrobeniusPCSParams {
  long ell = 0;
  long kappa = 0;
  long ell_prime = 0;
  NTL::ZZ base_modulus;
  NTL::ZZ_pX extension_modulus;
  FrobeniusBasisData basis_data;
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

}  // namespace basefold

#endif  // BASEFOLD_Z2K_FROBENIUS_PCS_HPP_
