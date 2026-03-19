#include "Compiler/Z2k/PCSBackend.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>

#include <string>

using NTL::LogicError;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pX;

namespace basefold {
namespace {

void ValidateBackendVTableOrThrow(const Z2kPCSBackendVTable &vtable,
                                  const char *func_name) {
  if (vtable.backend_name == nullptr) {
    LogicError((std::string(func_name) + ": backend_name must be set").c_str());
  }
  if (vtable.validate_params_or_throw == nullptr ||
      vtable.message_length == nullptr || vtable.point_dimension == nullptr ||
      vtable.commit == nullptr || vtable.build_commit_artifacts == nullptr ||
      vtable.commitment_from_artifacts == nullptr ||
      vtable.prove_eval == nullptr || vtable.verify_eval == nullptr ||
      vtable.serialize_eval_proof == nullptr ||
      vtable.eval_proof_size_bytes == nullptr) {
    LogicError((std::string(func_name) +
                ": backend vtable is missing required callbacks")
                   .c_str());
  }
}

void ValidateBackendHandleMetadataOrThrow(const Z2kPCSBackendHandle &backend,
                                          const char *func_name) {
  if (backend.vtable == nullptr) {
    LogicError((std::string(func_name) + ": backend handle is missing vtable")
                   .c_str());
  }
  ValidateBackendVTableOrThrow(*backend.vtable, func_name);
  if (!backend.params) {
    LogicError((std::string(func_name) + ": backend handle is missing params")
                   .c_str());
  }
}

void ValidateBackendHandleContextOrThrow(const Z2kPCSBackendHandle &backend,
                                         const char *func_name) {
  if (backend.context_base_modulus <= 1) {
    LogicError(
        (std::string(func_name) + ": backend handle is missing context snapshot")
            .c_str());
  }
  if (ZZ_p::modulus() != backend.context_base_modulus) {
    LogicError((std::string(func_name) +
                ": current ZZ_p modulus does not match backend context")
                   .c_str());
  }
  if (!ZZ_pE::initialized()) {
    LogicError((std::string(func_name) +
                ": ZZ_pE context must be initialized for this backend")
                   .c_str());
  }
  if (backend.context_extension_degree <= 0) {
    LogicError(
        (std::string(func_name) + ": backend handle has invalid context degree")
            .c_str());
  }
  if (ZZ_pE::degree() != backend.context_extension_degree) {
    LogicError((std::string(func_name) +
                ": current ZZ_pE degree does not match backend context")
                   .c_str());
  }
  if (ZZ_pE::modulus().val() != backend.context_extension_modulus) {
    LogicError((std::string(func_name) +
                ": current ZZ_pE modulus does not match backend context")
                   .c_str());
  }
}

void ValidateBackendHandleOrThrow(const Z2kPCSBackendHandle &backend,
                                  const char *func_name) {
  ValidateBackendHandleMetadataOrThrow(backend, func_name);
  ValidateBackendHandleContextOrThrow(backend, func_name);
}

void ValidateCommitArtifactsOrThrow(
    const Z2kPCSBackendHandle &backend,
    const Z2kPCSBackendCommitArtifacts &commit_artifacts,
    const char *func_name) {
  ValidateBackendHandleOrThrow(backend, func_name);
  if (commit_artifacts.vtable != backend.vtable) {
    LogicError((std::string(func_name) +
                ": commit artifacts backend does not match handle")
                   .c_str());
  }
  if (!commit_artifacts.payload) {
    LogicError((std::string(func_name) +
                ": commit artifacts payload must be present")
                   .c_str());
  }
  if (!commit_artifacts.params_owner) {
    LogicError((std::string(func_name) +
                ": commit artifacts params owner must be present")
                   .c_str());
  }
  if (commit_artifacts.params_owner.get() != backend.params.get()) {
    LogicError((std::string(func_name) +
                ": commit artifacts must belong to the same backend params instance")
                   .c_str());
  }
}

void ValidateEvalProofOrThrow(const Z2kPCSBackendHandle &backend,
                              const Z2kPCSBackendEvalProof &proof,
                              const char *func_name) {
  ValidateBackendHandleOrThrow(backend, func_name);
  if (proof.vtable != backend.vtable) {
    LogicError(
        (std::string(func_name) + ": proof backend does not match handle")
                   .c_str());
  }
  if (!proof.payload) {
    LogicError((std::string(func_name) + ": proof payload must be present")
                   .c_str());
  }
  if (!proof.params_owner) {
    LogicError((std::string(func_name) + ": proof params owner must be present")
                   .c_str());
  }
  if (proof.params_owner.get() != backend.params.get()) {
    LogicError((std::string(func_name) +
                ": proof must belong to the same backend params instance")
                   .c_str());
  }
}

}  // namespace

const char *Z2kPCSBackendName(const Z2kPCSBackendHandle &backend) {
  ValidateBackendHandleMetadataOrThrow(backend, "Z2kPCSBackendName");
  return backend.vtable->backend_name;
}

void Z2kPCSBackendValidateParamsOrThrow(const Z2kPCSBackendHandle &backend) {
  ValidateBackendHandleOrThrow(backend, "Z2kPCSBackendValidateParamsOrThrow");
  backend.vtable->validate_params_or_throw(backend.params.get());
}

long Z2kPCSBackendMessageLength(const Z2kPCSBackendHandle &backend) {
  ValidateBackendHandleOrThrow(backend, "Z2kPCSBackendMessageLength");
  return backend.vtable->message_length(backend.params.get());
}

long Z2kPCSBackendPointDimension(const Z2kPCSBackendHandle &backend) {
  ValidateBackendHandleOrThrow(backend, "Z2kPCSBackendPointDimension");
  return backend.vtable->point_dimension(backend.params.get());
}

MerkleRoot Z2kPCSBackendCommit(const Z2kPCSBackendHandle &backend,
                               const NTL::vec_ZZ_pE &f_coeffs) {
  ValidateBackendHandleOrThrow(backend, "Z2kPCSBackendCommit");
  return backend.vtable->commit(f_coeffs, backend.params.get());
}

Z2kPCSBackendCommitArtifacts Z2kPCSBackendBuildCommitArtifacts(
    const Z2kPCSBackendHandle &backend, const NTL::vec_ZZ_pE &f_coeffs) {
  ValidateBackendHandleOrThrow(backend, "Z2kPCSBackendBuildCommitArtifacts");
  Z2kPCSBackendCommitArtifacts out;
  out.vtable = backend.vtable;
  out.payload =
      backend.vtable->build_commit_artifacts(f_coeffs, backend.params.get());
  out.params_owner = backend.params;
  if (!out.payload) {
    LogicError(
        "Z2kPCSBackendBuildCommitArtifacts: backend returned empty payload");
  }
  out.commitment =
      backend.vtable->commitment_from_artifacts(out.payload.get(),
                                                backend.params.get());
  return out;
}

Z2kPCSBackendEvalProof Z2kPCSBackendProveEval(
    const Z2kPCSBackendHandle &backend, const NTL::vec_ZZ_pE &f_coeffs,
    const std::vector<FieldElement> &z, const FieldElement &claimed_y,
    long num_queries, const Z2kPCSBackendCommitArtifacts *commit_artifacts) {
  ValidateBackendHandleOrThrow(backend, "Z2kPCSBackendProveEval");
  const void *commit_payload = nullptr;
  if (commit_artifacts != nullptr) {
    ValidateCommitArtifactsOrThrow(backend, *commit_artifacts,
                                   "Z2kPCSBackendProveEval");
    commit_payload = commit_artifacts->payload.get();
  }

  Z2kPCSBackendEvalProof out;
  out.vtable = backend.vtable;
  out.payload = backend.vtable->prove_eval(
      f_coeffs, z, claimed_y, num_queries, backend.params.get(), commit_payload);
  out.params_owner = backend.params;
  if (!out.payload) {
    LogicError("Z2kPCSBackendProveEval: backend returned empty proof payload");
  }
  return out;
}

bool Z2kPCSBackendVerifyEval(const Z2kPCSBackendHandle &backend,
                             const MerkleRoot &commitment,
                             const std::vector<FieldElement> &z,
                             const FieldElement &claimed_y, long num_queries,
                             const Z2kPCSBackendEvalProof &proof) {
  ValidateEvalProofOrThrow(backend, proof, "Z2kPCSBackendVerifyEval");
  return backend.vtable->verify_eval(commitment, z, claimed_y, num_queries,
                                     proof.payload.get(),
                                     backend.params.get());
}

bool Z2kPCSBackendVerifyEvalUnchecked(const Z2kPCSBackendHandle &backend,
                                      const MerkleRoot &commitment,
                                      const std::vector<FieldElement> &z,
                                      const FieldElement &claimed_y,
                                      long num_queries,
                                      const Z2kPCSBackendEvalProof &proof) {
  ValidateEvalProofOrThrow(backend, proof, "Z2kPCSBackendVerifyEvalUnchecked");
  if (backend.vtable->verify_eval_unchecked != nullptr) {
    return backend.vtable->verify_eval_unchecked(
        commitment, z, claimed_y, num_queries, proof.payload.get(),
        backend.params.get());
  }
  return backend.vtable->verify_eval(commitment, z, claimed_y, num_queries,
                                     proof.payload.get(),
                                     backend.params.get());
}

Bytes Z2kPCSBackendSerializeEvalProof(
    const Z2kPCSBackendHandle &backend, const Z2kPCSBackendEvalProof &proof,
    const Z2kPCSBackendProofEncodingOptions &options) {
  ValidateEvalProofOrThrow(backend, proof, "Z2kPCSBackendSerializeEvalProof");
  return backend.vtable->serialize_eval_proof(proof.payload.get(), options);
}

std::uint64_t Z2kPCSBackendEvalProofSizeBytes(
    const Z2kPCSBackendHandle &backend, const Z2kPCSBackendEvalProof &proof,
    const Z2kPCSBackendProofSizeOptions &options) {
  ValidateEvalProofOrThrow(backend, proof, "Z2kPCSBackendEvalProofSizeBytes");
  return backend.vtable->eval_proof_size_bytes(proof.payload.get(), options);
}

}  // namespace basefold
