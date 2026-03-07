#ifndef BASEFOLD_PROOFSERIALIZE_HPP_
#define BASEFOLD_PROOFSERIALIZE_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>

#include <cstdint>
#include <limits>
#include <vector>

#include "BaseFold/BaseFoldPCS.hpp"

namespace basefold {

struct FixedProofEncodingOptions {
  bool include_version_byte = true;
  long challenge_ext_degree = 0;
};

struct FixedProofEncodingContext {
  std::uint8_t format_version = 1;
  std::uint64_t hash_bytes = 32;
  std::uint64_t u64_bytes = 8;
  std::uint64_t coeff_bytes = 0;
  std::uint64_t field_elem_bytes = 0;
  std::uint64_t extension_elem_bytes = 0;
  long base_ext_degree = 0;
  long challenge_ext_degree = 0;
  bool include_version_byte = true;
};

class CountingSink {
 public:
  inline void WriteU8(std::uint8_t value) {
    (void)value;
    WriteFixedBytes(1);
  }

  inline void WriteU64(std::uint64_t value) {
    (void)value;
    WriteFixedBytes(8);
  }

  inline void WriteFixedBytes(std::uint64_t size) {
    if (size > std::numeric_limits<std::uint64_t>::max() - total_) {
      NTL::LogicError("CountingSink::WriteFixedBytes: byte count overflow");
    }
    total_ += size;
  }

  inline std::uint64_t bytes_written() const { return total_; }

 private:
  std::uint64_t total_ = 0;
};

namespace fixed_proof_serialize_detail {

inline std::uint64_t SizeToU64OrThrow(std::size_t value,
                                      const char *what) {
  if (value > static_cast<std::size_t>(std::numeric_limits<std::uint64_t>::max())) {
    NTL::LogicError(what);
  }
  return static_cast<std::uint64_t>(value);
}

inline std::uint64_t LongToU64OrThrow(long value, const char *what) {
  if (value < 0) {
    NTL::LogicError(what);
  }
  return static_cast<std::uint64_t>(value);
}

inline std::uint64_t MulU64OrThrow(std::uint64_t lhs, std::uint64_t rhs,
                                   const char *what) {
  if (lhs == 0 || rhs == 0) {
    return 0;
  }
  if (lhs > std::numeric_limits<std::uint64_t>::max() / rhs) {
    NTL::LogicError(what);
  }
  return lhs * rhs;
}

inline std::uint64_t ComputeFixedCoeffBytesOrThrow() {
  const NTL::ZZ mod_minus_one = NTL::ZZ_p::modulus() - 1;
  const long bits = NTL::NumBits(mod_minus_one);
  if (bits <= 0) {
    NTL::LogicError(
        "ComputeFixedCoeffBytesOrThrow: invalid base modulus bit width");
  }
  return static_cast<std::uint64_t>((bits + 7) / 8);
}

inline void ValidateExtensionWidthOrThrow(
    const BaseFoldPCSEvalProof &proof,
    const FixedProofEncodingOptions &options) {
  if (!proof.extension.enabled) {
    if (options.challenge_ext_degree < 0) {
      NTL::LogicError(
          "ValidateExtensionWidthOrThrow: challenge_ext_degree must be >= 0");
    }
    return;
  }
  if (options.challenge_ext_degree <= 0) {
    NTL::LogicError(
        "ValidateExtensionWidthOrThrow: extension proof requires "
        "challenge_ext_degree > 0");
  }
}

inline void SerializeVersion(CountingSink &sink,
                             const FixedProofEncodingContext &ctx) {
  if (ctx.include_version_byte) {
    sink.WriteU8(ctx.format_version);
  }
}

inline void SerializeVecHeader(CountingSink &sink, std::uint64_t count) {
  sink.WriteU64(count);
}

inline void SerializeDigestFixed(CountingSink &sink,
                                 const FixedProofEncodingContext &ctx) {
  sink.WriteFixedBytes(ctx.hash_bytes);
}

inline void SerializeMerkleRootFixed(CountingSink &sink,
                                     const MerkleRoot &root,
                                     const FixedProofEncodingContext &ctx) {
  (void)root;
  SerializeDigestFixed(sink, ctx);
}

inline void SerializeFieldElementFixed(CountingSink &sink,
                                       const FieldElement &value,
                                       const FixedProofEncodingContext &ctx) {
  const long deg = NTL::deg(NTL::rep(value));
  if (deg >= ctx.base_ext_degree) {
    NTL::LogicError(
        "SerializeFieldElementFixed: value degree exceeds base extension "
        "degree");
  }
  sink.WriteFixedBytes(ctx.field_elem_bytes);
}

inline void SerializeExtensionElementFixed(CountingSink &sink,
                                           const NTL::ZZ_pEX &value,
                                           const FixedProofEncodingContext &ctx) {
  if (ctx.challenge_ext_degree <= 0) {
    NTL::LogicError(
        "SerializeExtensionElementFixed: challenge_ext_degree must be > 0");
  }
  const long deg = NTL::deg(value);
  if (deg >= ctx.challenge_ext_degree) {
    NTL::LogicError(
        "SerializeExtensionElementFixed: extension value degree exceeds "
        "challenge extension degree");
  }
  sink.WriteFixedBytes(ctx.extension_elem_bytes);
}

inline void SerializeQuadraticPolyFixed(CountingSink &sink,
                                        const QuadraticPoly &poly,
                                        const FixedProofEncodingContext &ctx) {
  SerializeFieldElementFixed(sink, poly.a0, ctx);
  SerializeFieldElementFixed(sink, poly.a1, ctx);
  SerializeFieldElementFixed(sink, poly.a2, ctx);
}

inline void SerializeExtensionQuadraticPolyFixed(
    CountingSink &sink, const ExtensionQuadraticPoly &poly,
    const FixedProofEncodingContext &ctx) {
  SerializeExtensionElementFixed(sink, poly.a0, ctx);
  SerializeExtensionElementFixed(sink, poly.a1, ctx);
  SerializeExtensionElementFixed(sink, poly.a2, ctx);
}

inline void SerializeOracleFixed(CountingSink &sink, const Oracle &oracle,
                                 const FixedProofEncodingContext &ctx) {
  const std::uint64_t count = LongToU64OrThrow(
      oracle.length(), "SerializeOracleFixed: oracle length must be >= 0");
  SerializeVecHeader(sink, count);
  for (long i = 0; i < oracle.length(); ++i) {
    SerializeFieldElementFixed(sink, oracle[i], ctx);
  }
}

inline void SerializeMerkleAuthPathFixed(CountingSink &sink,
                                         const MerkleAuthPath &auth_path,
                                         const FixedProofEncodingContext &ctx) {
  const std::uint64_t count = SizeToU64OrThrow(
      auth_path.sibling_hashes.size(),
      "SerializeMerkleAuthPathFixed: sibling hash count overflow");
  SerializeVecHeader(sink, count);
  for (const Digest &digest : auth_path.sibling_hashes) {
    (void)digest;
    SerializeDigestFixed(sink, ctx);
  }
}

inline void SerializeMerkleOpeningFixed(CountingSink &sink,
                                        const MerkleOpening &opening,
                                        const FixedProofEncodingContext &ctx) {
  sink.WriteU64(LongToU64OrThrow(
      opening.index, "SerializeMerkleOpeningFixed: index must be >= 0"));
  SerializeFieldElementFixed(sink, opening.value, ctx);
  SerializeMerkleAuthPathFixed(sink, opening.auth_path, ctx);
}

inline void SerializeExtensionMerkleOpeningFixed(
    CountingSink &sink, const ExtensionMerkleOpening &opening,
    const FixedProofEncodingContext &ctx) {
  sink.WriteU64(LongToU64OrThrow(
      opening.index,
      "SerializeExtensionMerkleOpeningFixed: index must be >= 0"));
  SerializeExtensionElementFixed(sink, opening.value, ctx);
  SerializeMerkleAuthPathFixed(sink, opening.auth_path, ctx);
}

inline void SerializeMerkleMultiproofFixed(CountingSink &sink,
                                           const MerkleMultiproof &proof,
                                           const FixedProofEncodingContext &ctx) {
  const std::uint64_t query_count =
      SizeToU64OrThrow(proof.queried_indices.size(),
                       "SerializeMerkleMultiproofFixed: queried index count overflow");
  SerializeVecHeader(sink, query_count);
  for (long index : proof.queried_indices) {
    sink.WriteU64(LongToU64OrThrow(
        index, "SerializeMerkleMultiproofFixed: queried index must be >= 0"));
  }

  if (proof.values.length() < 0) {
    NTL::LogicError(
        "SerializeMerkleMultiproofFixed: values length must be >= 0");
  }
  const std::uint64_t value_count = static_cast<std::uint64_t>(proof.values.length());
  if (value_count != query_count) {
    NTL::LogicError(
        "SerializeMerkleMultiproofFixed: values length must match queried "
        "indices count");
  }
  for (long i = 0; i < proof.values.length(); ++i) {
    SerializeFieldElementFixed(sink, proof.values[i], ctx);
  }

  const std::uint64_t sibling_count = SizeToU64OrThrow(
      proof.sibling_hashes.size(),
      "SerializeMerkleMultiproofFixed: sibling hash count overflow");
  SerializeVecHeader(sink, sibling_count);
  for (const Digest &digest : proof.sibling_hashes) {
    (void)digest;
    SerializeDigestFixed(sink, ctx);
  }
}

inline void SerializeLegacyQueryProofFixed(
    CountingSink &sink, const BaseFoldPCSQueryProof &query_proof,
    const FixedProofEncodingContext &ctx) {
  const std::uint64_t left_count =
      SizeToU64OrThrow(query_proof.left.size(),
                       "SerializeLegacyQueryProofFixed: left count overflow");
  SerializeVecHeader(sink, left_count);
  for (const MerkleOpening &opening : query_proof.left) {
    SerializeMerkleOpeningFixed(sink, opening, ctx);
  }

  const std::uint64_t right_count =
      SizeToU64OrThrow(query_proof.right.size(),
                       "SerializeLegacyQueryProofFixed: right count overflow");
  SerializeVecHeader(sink, right_count);
  for (const MerkleOpening &opening : query_proof.right) {
    SerializeMerkleOpeningFixed(sink, opening, ctx);
  }

  const std::uint64_t folded_count =
      SizeToU64OrThrow(query_proof.folded.size(),
                       "SerializeLegacyQueryProofFixed: folded count overflow");
  SerializeVecHeader(sink, folded_count);
  for (const MerkleOpening &opening : query_proof.folded) {
    SerializeMerkleOpeningFixed(sink, opening, ctx);
  }
}

inline void SerializeExtensionQueryProofFixed(
    CountingSink &sink, const BaseFoldPCSQueryProofExtension &query_proof,
    const FixedProofEncodingContext &ctx) {
  const std::uint64_t left_count = SizeToU64OrThrow(
      query_proof.left.size(),
      "SerializeExtensionQueryProofFixed: left count overflow");
  SerializeVecHeader(sink, left_count);
  for (const ExtensionMerkleOpening &opening : query_proof.left) {
    SerializeExtensionMerkleOpeningFixed(sink, opening, ctx);
  }

  const std::uint64_t right_count = SizeToU64OrThrow(
      query_proof.right.size(),
      "SerializeExtensionQueryProofFixed: right count overflow");
  SerializeVecHeader(sink, right_count);
  for (const ExtensionMerkleOpening &opening : query_proof.right) {
    SerializeExtensionMerkleOpeningFixed(sink, opening, ctx);
  }

  const std::uint64_t folded_count = SizeToU64OrThrow(
      query_proof.folded.size(),
      "SerializeExtensionQueryProofFixed: folded count overflow");
  SerializeVecHeader(sink, folded_count);
  for (const ExtensionMerkleOpening &opening : query_proof.folded) {
    SerializeExtensionMerkleOpeningFixed(sink, opening, ctx);
  }
}

inline void SerializeCommitmentsFixed(CountingSink &sink,
                                      const IOPPMerkleCommitments &commitments,
                                      const FixedProofEncodingContext &ctx) {
  const std::uint64_t count =
      SizeToU64OrThrow(commitments.roots_by_level.size(),
                       "SerializeCommitmentsFixed: commitment count overflow");
  SerializeVecHeader(sink, count);
  for (const MerkleRoot &root : commitments.roots_by_level) {
    SerializeMerkleRootFixed(sink, root, ctx);
  }
}

inline void SerializeExtensionProofDataFixed(
    CountingSink &sink, const BaseFoldPCSExtensionProofData &extension,
    const FixedProofEncodingContext &ctx) {
  const std::uint64_t roots_count =
      SizeToU64OrThrow(extension.roots_by_level.size(),
                       "SerializeExtensionProofDataFixed: extension roots count overflow");
  SerializeVecHeader(sink, roots_count);
  for (const MerkleRoot &root : extension.roots_by_level) {
    SerializeMerkleRootFixed(sink, root, ctx);
  }

  const std::uint64_t h_count =
      SizeToU64OrThrow(extension.h_by_level.size(),
                       "SerializeExtensionProofDataFixed: extension h count overflow");
  SerializeVecHeader(sink, h_count);
  for (const ExtensionQuadraticPoly &h : extension.h_by_level) {
    SerializeExtensionQuadraticPolyFixed(sink, h, ctx);
  }

  const std::uint64_t r_count =
      SizeToU64OrThrow(extension.r_by_level.size(),
                       "SerializeExtensionProofDataFixed: extension r count overflow");
  SerializeVecHeader(sink, r_count);
  for (const NTL::ZZ_pEX &r : extension.r_by_level) {
    SerializeExtensionElementFixed(sink, r, ctx);
  }

  const std::uint64_t msg0_count =
      SizeToU64OrThrow(extension.msg0_coeffs.size(),
                       "SerializeExtensionProofDataFixed: extension msg0 count overflow");
  SerializeVecHeader(sink, msg0_count);
  for (const NTL::ZZ_pEX &coeff : extension.msg0_coeffs) {
    SerializeExtensionElementFixed(sink, coeff, ctx);
  }

  const std::uint64_t pi0_count =
      SizeToU64OrThrow(extension.pi0_full.size(),
                       "SerializeExtensionProofDataFixed: extension pi0 count overflow");
  SerializeVecHeader(sink, pi0_count);
  for (const NTL::ZZ_pEX &value : extension.pi0_full) {
    SerializeExtensionElementFixed(sink, value, ctx);
  }

  const std::uint64_t query_count =
      SizeToU64OrThrow(extension.query_proofs.size(),
                       "SerializeExtensionProofDataFixed: extension query count overflow");
  SerializeVecHeader(sink, query_count);
  for (const BaseFoldPCSQueryProofExtension &query : extension.query_proofs) {
    SerializeExtensionQueryProofFixed(sink, query, ctx);
  }
}

}  // namespace fixed_proof_serialize_detail

inline FixedProofEncodingContext BuildFixedProofEncodingContext(
    const BaseFoldPCSEvalProof &proof,
    const FixedProofEncodingOptions &options) {
  fixed_proof_serialize_detail::ValidateExtensionWidthOrThrow(proof, options);

  FixedProofEncodingContext ctx;
  ctx.include_version_byte = options.include_version_byte;
  ctx.base_ext_degree = NTL::ZZ_pE::degree();
  if (ctx.base_ext_degree <= 0) {
    NTL::LogicError(
        "BuildFixedProofEncodingContext: invalid base extension degree");
  }
  ctx.challenge_ext_degree = options.challenge_ext_degree;
  ctx.coeff_bytes = fixed_proof_serialize_detail::ComputeFixedCoeffBytesOrThrow();
  ctx.field_elem_bytes = fixed_proof_serialize_detail::MulU64OrThrow(
      ctx.coeff_bytes,
      fixed_proof_serialize_detail::LongToU64OrThrow(
          ctx.base_ext_degree,
          "BuildFixedProofEncodingContext: base extension degree must be >= 0"),
      "BuildFixedProofEncodingContext: field element byte width overflow");

  if (ctx.challenge_ext_degree > 0) {
    ctx.extension_elem_bytes = fixed_proof_serialize_detail::MulU64OrThrow(
        ctx.field_elem_bytes,
        fixed_proof_serialize_detail::LongToU64OrThrow(
            ctx.challenge_ext_degree,
            "BuildFixedProofEncodingContext: challenge extension degree must be "
            ">= 0"),
        "BuildFixedProofEncodingContext: extension element byte width overflow");
  }
  return ctx;
}

inline void SerializeBaseFoldPCSEvalProofFixed(
    CountingSink &sink, const BaseFoldPCSEvalProof &proof,
    const FixedProofEncodingContext &ctx) {
  fixed_proof_serialize_detail::SerializeVersion(sink, ctx);

  fixed_proof_serialize_detail::SerializeCommitmentsFixed(sink, proof.commitments,
                                                          ctx);

  const std::uint64_t h_count = fixed_proof_serialize_detail::SizeToU64OrThrow(
      proof.h_by_level.size(),
      "SerializeBaseFoldPCSEvalProofFixed: h_by_level count overflow");
  fixed_proof_serialize_detail::SerializeVecHeader(sink, h_count);
  for (const QuadraticPoly &h : proof.h_by_level) {
    fixed_proof_serialize_detail::SerializeQuadraticPolyFixed(sink, h, ctx);
  }

  fixed_proof_serialize_detail::SerializeOracleFixed(sink, proof.pi0_full, ctx);

  const std::uint64_t multiproof_count =
      fixed_proof_serialize_detail::SizeToU64OrThrow(
          proof.query_multiproofs.size(),
          "SerializeBaseFoldPCSEvalProofFixed: query_multiproofs count overflow");
  fixed_proof_serialize_detail::SerializeVecHeader(sink, multiproof_count);
  for (const MerkleMultiproof &multiproof : proof.query_multiproofs) {
    fixed_proof_serialize_detail::SerializeMerkleMultiproofFixed(sink, multiproof,
                                                                 ctx);
  }

  const std::uint64_t query_count = fixed_proof_serialize_detail::SizeToU64OrThrow(
      proof.query_proofs.size(),
      "SerializeBaseFoldPCSEvalProofFixed: query_proofs count overflow");
  fixed_proof_serialize_detail::SerializeVecHeader(sink, query_count);
  for (const BaseFoldPCSQueryProof &query : proof.query_proofs) {
    fixed_proof_serialize_detail::SerializeLegacyQueryProofFixed(sink, query,
                                                                 ctx);
  }

  sink.WriteU8(proof.extension.enabled ? 1 : 0);
  if (proof.extension.enabled) {
    fixed_proof_serialize_detail::SerializeExtensionProofDataFixed(
        sink, proof.extension, ctx);
  }
}

inline std::uint64_t CountSerializedBaseFoldPCSEvalProofFixedBytes(
    const BaseFoldPCSEvalProof &proof,
    const FixedProofEncodingOptions &options) {
  const FixedProofEncodingContext ctx =
      BuildFixedProofEncodingContext(proof, options);
  CountingSink sink;
  SerializeBaseFoldPCSEvalProofFixed(sink, proof, ctx);
  return sink.bytes_written();
}

}  // namespace basefold

#endif  // BASEFOLD_PROOFSERIALIZE_HPP_
