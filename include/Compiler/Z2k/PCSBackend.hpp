#ifndef BASEFOLD_Z2KPCS_BACKEND_HPP_
#define BASEFOLD_Z2KPCS_BACKEND_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_pE.h>

#include <cstdint>
#include <memory>
#include <vector>

#include "PCS/Common/Hash.hpp"
#include "PCS/Common/Merkle.hpp"

namespace basefold {

struct Z2kPCSBackendProofSizeOptions {
  bool include_version_byte = true;
  long challenge_ext_degree = 0;
};

struct Z2kPCSBackendProofEncodingOptions {
  bool include_version_byte = true;
  long challenge_ext_degree = 0;
};

using Z2kPCSBackendOpaquePtr = std::shared_ptr<const void>;

struct Z2kPCSBackendVTable {
  const char *backend_name = nullptr;

  void (*validate_params_or_throw)(const void *backend_params) = nullptr;
  long (*message_length)(const void *backend_params) = nullptr;
  long (*point_dimension)(const void *backend_params) = nullptr;

  MerkleRoot (*commit)(const NTL::vec_ZZ_pE &f_coeffs,
                       const void *backend_params) = nullptr;
  Z2kPCSBackendOpaquePtr (*build_commit_artifacts)(
      const NTL::vec_ZZ_pE &f_coeffs, const void *backend_params) = nullptr;
  MerkleRoot (*commitment_from_artifacts)(
      const void *commit_artifacts,
      const void *backend_params) = nullptr;
  Z2kPCSBackendOpaquePtr (*prove_eval)(
      const NTL::vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
      const FieldElement &claimed_y, long num_queries,
      const void *backend_params, const void *commit_artifacts_or_null) = nullptr;
  bool (*verify_eval)(const MerkleRoot &commitment,
                      const std::vector<FieldElement> &z,
                      const FieldElement &claimed_y, long num_queries,
                      const void *backend_proof,
                      const void *backend_params) = nullptr;
  bool (*verify_eval_unchecked)(const MerkleRoot &commitment,
                                const std::vector<FieldElement> &z,
                                const FieldElement &claimed_y,
                                long num_queries,
                                const void *backend_proof,
                                const void *backend_params) = nullptr;
  Bytes (*serialize_eval_proof)(
      const void *backend_proof,
      const Z2kPCSBackendProofEncodingOptions &options) = nullptr;
  std::uint64_t (*eval_proof_size_bytes)(
      const void *backend_proof,
      const Z2kPCSBackendProofSizeOptions &options) = nullptr;
};

struct Z2kPCSBackendHandle {
  const Z2kPCSBackendVTable *vtable = nullptr;
  Z2kPCSBackendOpaquePtr params;
  NTL::ZZ context_base_modulus;
  NTL::ZZ_pX context_extension_modulus;
  long context_extension_degree = 0;
};

struct Z2kPCSBackendCommitArtifacts {
  const Z2kPCSBackendVTable *vtable = nullptr;
  Z2kPCSBackendOpaquePtr payload;
  Z2kPCSBackendOpaquePtr params_owner;
  MerkleRoot commitment;
};

struct Z2kPCSBackendEvalProof {
  const Z2kPCSBackendVTable *vtable = nullptr;
  Z2kPCSBackendOpaquePtr payload;
  Z2kPCSBackendOpaquePtr params_owner;
};

const char *Z2kPCSBackendName(const Z2kPCSBackendHandle &backend);

void Z2kPCSBackendValidateParamsOrThrow(const Z2kPCSBackendHandle &backend);

long Z2kPCSBackendMessageLength(const Z2kPCSBackendHandle &backend);

long Z2kPCSBackendPointDimension(const Z2kPCSBackendHandle &backend);

MerkleRoot Z2kPCSBackendCommit(const Z2kPCSBackendHandle &backend,
                               const NTL::vec_ZZ_pE &f_coeffs);

Z2kPCSBackendCommitArtifacts Z2kPCSBackendBuildCommitArtifacts(
    const Z2kPCSBackendHandle &backend, const NTL::vec_ZZ_pE &f_coeffs);

Z2kPCSBackendEvalProof Z2kPCSBackendProveEval(
    const Z2kPCSBackendHandle &backend, const NTL::vec_ZZ_pE &f_coeffs,
    const std::vector<FieldElement> &z, const FieldElement &claimed_y,
    long num_queries,
    const Z2kPCSBackendCommitArtifacts *commit_artifacts = nullptr);

bool Z2kPCSBackendVerifyEval(const Z2kPCSBackendHandle &backend,
                             const MerkleRoot &commitment,
                             const std::vector<FieldElement> &z,
                             const FieldElement &claimed_y, long num_queries,
                             const Z2kPCSBackendEvalProof &proof);

bool Z2kPCSBackendVerifyEvalUnchecked(const Z2kPCSBackendHandle &backend,
                                      const MerkleRoot &commitment,
                                      const std::vector<FieldElement> &z,
                                      const FieldElement &claimed_y,
                                      long num_queries,
                                      const Z2kPCSBackendEvalProof &proof);

Bytes Z2kPCSBackendSerializeEvalProof(
    const Z2kPCSBackendHandle &backend, const Z2kPCSBackendEvalProof &proof,
    const Z2kPCSBackendProofEncodingOptions &options = {});

std::uint64_t Z2kPCSBackendEvalProofSizeBytes(
    const Z2kPCSBackendHandle &backend, const Z2kPCSBackendEvalProof &proof,
    const Z2kPCSBackendProofSizeOptions &options = {});

}  // namespace basefold

#endif  // BASEFOLD_Z2KPCS_BACKEND_HPP_
