#include "Compiler/Z2k/BaseFoldBackendAdapter.hpp"

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>

#include <memory>
#include <string>

#include "PCS/BaseFold/ProofSerialize.hpp"
#include "PCS/BaseFold/ProofSize.hpp"

using NTL::LogicError;
using NTL::ZZ_p;
using NTL::ZZ_pE;

namespace basefold {
namespace {

long Log2ExactPowerOfTwoLongOrThrow(long value, const char *what) {
  if (value <= 0 || (value & (value - 1)) != 0) {
    LogicError(what);
  }
  long out = 0;
  while (value > 1) {
    value >>= 1;
    ++out;
  }
  return out;
}

const FoldableCodeParams &AsBaseFoldParams(const void *backend_params,
                                           const char *func_name) {
  if (backend_params == nullptr) {
    LogicError((std::string(func_name) + ": backend params must be present")
                   .c_str());
  }
  return *static_cast<const FoldableCodeParams *>(backend_params);
}

const BaseFoldPCSCommitArtifacts &AsBaseFoldCommitArtifacts(
    const void *commit_artifacts, const char *func_name) {
  if (commit_artifacts == nullptr) {
    LogicError((std::string(func_name) +
                ": commit artifacts payload must be present")
                   .c_str());
  }
  return *static_cast<const BaseFoldPCSCommitArtifacts *>(commit_artifacts);
}

const BaseFoldPCSEvalProof &AsBaseFoldEvalProof(const void *backend_proof,
                                                const char *func_name) {
  if (backend_proof == nullptr) {
    LogicError((std::string(func_name) + ": backend proof must be present")
                   .c_str());
  }
  return *static_cast<const BaseFoldPCSEvalProof *>(backend_proof);
}

void BaseFoldBackendValidateParamsOrThrow(const void *backend_params) {
  const FoldableCodeParams &params =
      AsBaseFoldParams(backend_params, "BaseFoldBackendValidateParamsOrThrow");
  (void)MessageLength(params);
}

long BaseFoldBackendMessageLength(const void *backend_params) {
  const FoldableCodeParams &params =
      AsBaseFoldParams(backend_params, "BaseFoldBackendMessageLength");
  return MessageLength(params);
}

long BaseFoldBackendPointDimension(const void *backend_params) {
  const FoldableCodeParams &params =
      AsBaseFoldParams(backend_params, "BaseFoldBackendPointDimension");
  const long log_k0 = Log2ExactPowerOfTwoLongOrThrow(
      params.k0, "BaseFoldBackendPointDimension: params.k0 must be a power of two");
  const long point_dim = params.d + log_k0;
  if (point_dim < params.d) {
    LogicError("BaseFoldBackendPointDimension: dimension overflow");
  }
  return point_dim;
}

MerkleRoot BaseFoldBackendCommit(const NTL::vec_ZZ_pE &f_coeffs,
                                 const void *backend_params) {
  const FoldableCodeParams &params =
      AsBaseFoldParams(backend_params, "BaseFoldBackendCommit");
  return BaseFoldPCSCommit(f_coeffs, params);
}

Z2kPCSBackendOpaquePtr BaseFoldBackendBuildCommitArtifacts(
    const NTL::vec_ZZ_pE &f_coeffs, const void *backend_params) {
  const FoldableCodeParams &params =
      AsBaseFoldParams(backend_params, "BaseFoldBackendBuildCommitArtifacts");
  return std::static_pointer_cast<const void>(
      std::make_shared<BaseFoldPCSCommitArtifacts>(
          BaseFoldPCSBuildCommitArtifacts(f_coeffs, params)));
}

MerkleRoot BaseFoldBackendCommitmentFromArtifacts(const void *commit_artifacts,
                                                  const void *backend_params) {
  (void)backend_params;
  const BaseFoldPCSCommitArtifacts &artifacts = AsBaseFoldCommitArtifacts(
      commit_artifacts, "BaseFoldBackendCommitmentFromArtifacts");
  return artifacts.root_d;
}

Z2kPCSBackendOpaquePtr BaseFoldBackendProveEval(
    const NTL::vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries, const void *backend_params,
    const void *commit_artifacts_or_null) {
  const FoldableCodeParams &params =
      AsBaseFoldParams(backend_params, "BaseFoldBackendProveEval");
  BaseFoldPCSEvalProof proof;
  if (commit_artifacts_or_null == nullptr) {
    proof = BaseFoldPCSProveEval(f_coeffs, z, claimed_y, num_queries, params);
  } else {
    const BaseFoldPCSCommitArtifacts &commit_artifacts =
        AsBaseFoldCommitArtifacts(commit_artifacts_or_null,
                                  "BaseFoldBackendProveEval");
    proof = BaseFoldPCSProveEvalFromCommittedTopOracle(
        f_coeffs, z, claimed_y, num_queries, params, commit_artifacts);
  }
  return std::static_pointer_cast<const void>(
      std::make_shared<BaseFoldPCSEvalProof>(std::move(proof)));
}

bool BaseFoldBackendVerifyEval(const MerkleRoot &commitment,
                               const std::vector<FieldElement> &z,
                               const FieldElement &claimed_y, long num_queries,
                               const void *backend_proof,
                               const void *backend_params) {
  const FoldableCodeParams &params =
      AsBaseFoldParams(backend_params, "BaseFoldBackendVerifyEval");
  const BaseFoldPCSEvalProof &proof =
      AsBaseFoldEvalProof(backend_proof, "BaseFoldBackendVerifyEval");
  return BaseFoldPCSVerifyEval(commitment, z, claimed_y, num_queries, proof,
                               params);
}

Bytes BaseFoldBackendSerializeEvalProof(
    const void *backend_proof,
    const Z2kPCSBackendProofEncodingOptions &options) {
  const BaseFoldPCSEvalProof &proof =
      AsBaseFoldEvalProof(backend_proof, "BaseFoldBackendSerializeEvalProof");
  FixedProofEncodingOptions basefold_options;
  basefold_options.include_version_byte = options.include_version_byte;
  basefold_options.challenge_ext_degree = options.challenge_ext_degree;
  return SerializeBaseFoldPCSEvalProofFixedBytes(proof, basefold_options);
}

std::uint64_t BaseFoldBackendEvalProofSizeBytes(
    const void *backend_proof, const Z2kPCSBackendProofSizeOptions &options) {
  const BaseFoldPCSEvalProof &proof =
      AsBaseFoldEvalProof(backend_proof, "BaseFoldBackendEvalProofSizeBytes");
  BaseFoldProofSizeOptions basefold_options;
  basefold_options.include_version_byte = options.include_version_byte;
  basefold_options.challenge_ext_degree = options.challenge_ext_degree;
  return BaseFoldPCSEvalProofSizeBytes(proof, basefold_options);
}

const Z2kPCSBackendVTable kBaseFoldBackendVTable = {
    "basefold",
    &BaseFoldBackendValidateParamsOrThrow,
    &BaseFoldBackendMessageLength,
    &BaseFoldBackendPointDimension,
    &BaseFoldBackendCommit,
    &BaseFoldBackendBuildCommitArtifacts,
    &BaseFoldBackendCommitmentFromArtifacts,
    &BaseFoldBackendProveEval,
    &BaseFoldBackendVerifyEval,
    &BaseFoldBackendSerializeEvalProof,
    &BaseFoldBackendEvalProofSizeBytes,
};

}  // namespace

Z2kPCSBackendHandle MakeBaseFoldZ2kPCSBackend(
    const FoldableCodeParams &params) {
  if (ZZ_p::modulus() <= 1) {
    LogicError("MakeBaseFoldZ2kPCSBackend: ZZ_p context must be initialized");
  }
  if (!ZZ_pE::initialized()) {
    LogicError("MakeBaseFoldZ2kPCSBackend: ZZ_pE context must be initialized");
  }
  if (ZZ_pE::degree() <= 0) {
    LogicError("MakeBaseFoldZ2kPCSBackend: current ZZ_pE degree must be positive");
  }
  Z2kPCSBackendHandle backend;
  backend.vtable = &kBaseFoldBackendVTable;
  backend.params =
      std::static_pointer_cast<const void>(
          std::make_shared<FoldableCodeParams>(params));
  backend.context_base_modulus = ZZ_p::modulus();
  backend.context_extension_modulus = ZZ_pE::modulus().val();
  backend.context_extension_degree = ZZ_pE::degree();
  return backend;
}

}  // namespace basefold
