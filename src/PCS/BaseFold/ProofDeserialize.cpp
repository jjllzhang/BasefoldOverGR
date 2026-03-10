#include "PCS/BaseFold/ProofDeserialize.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pX.h>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace basefold {
namespace {

using fixed_proof_serialize_detail::LongToU64OrThrow;
using fixed_proof_serialize_detail::MulU64OrThrow;

class ByteReader {
 public:
  ByteReader(const Byte *data, std::size_t size) : data_(data), size_(size) {}

  std::size_t remaining() const { return size_ - pos_; }

  std::uint8_t ReadU8(const char *what) {
    RequireBytes(1, what);
    return data_[pos_++];
  }

  std::uint64_t ReadU64(const char *what) {
    RequireBytes(8, what);
    std::uint64_t out = 0;
    for (int i = 0; i < 8; ++i) {
      out |= static_cast<std::uint64_t>(data_[pos_++]) << (8 * i);
    }
    return out;
  }

  void ReadBytes(Byte *dst, std::size_t size, const char *what) {
    if (size == 0) {
      return;
    }
    if (dst == nullptr) {
      NTL::LogicError("ByteReader::ReadBytes: destination must be present");
    }
    RequireBytes(size, what);
    std::copy(data_ + pos_, data_ + pos_ + size, dst);
    pos_ += size;
  }

 private:
  void RequireBytes(std::size_t need, const char *what) const {
    if (need > remaining()) {
      NTL::LogicError(what);
    }
  }

  const Byte *data_ = nullptr;
  std::size_t size_ = 0;
  std::size_t pos_ = 0;
};

long U64ToLongOrThrow(std::uint64_t value, const char *what) {
  if (value > static_cast<std::uint64_t>(std::numeric_limits<long>::max())) {
    NTL::LogicError(what);
  }
  return static_cast<long>(value);
}

std::size_t U64ToSizeOrThrow(std::uint64_t value, const char *what) {
  if (value > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    NTL::LogicError(what);
  }
  return static_cast<std::size_t>(value);
}

std::uint64_t ComputeExpectedFieldElementBytesOrThrow(
    const FixedProofEncodingContext &ctx) {
  return MulU64OrThrow(
      ctx.coeff_bytes,
      LongToU64OrThrow(ctx.base_ext_degree,
                       "ComputeExpectedFieldElementBytesOrThrow: invalid base extension degree"),
      "ComputeExpectedFieldElementBytesOrThrow: byte width overflow");
}

std::uint64_t ComputeExpectedExtensionElementBytesOrThrow(
    const FixedProofEncodingContext &ctx) {
  return MulU64OrThrow(
      ComputeExpectedFieldElementBytesOrThrow(ctx),
      LongToU64OrThrow(ctx.challenge_ext_degree,
                       "ComputeExpectedExtensionElementBytesOrThrow: invalid challenge extension degree"),
      "ComputeExpectedExtensionElementBytesOrThrow: byte width overflow");
}

void ValidateReadWidthsOrThrow(const FixedProofEncodingContext &ctx) {
  if (ctx.base_ext_degree <= 0) {
    NTL::LogicError("ValidateReadWidthsOrThrow: invalid base extension degree");
  }
  if (ctx.coeff_bytes == 0) {
    NTL::LogicError("ValidateReadWidthsOrThrow: coefficient byte width must be positive");
  }
  if (ctx.field_elem_bytes != ComputeExpectedFieldElementBytesOrThrow(ctx)) {
    NTL::LogicError("ValidateReadWidthsOrThrow: field element byte width mismatch");
  }
  if (ctx.challenge_ext_degree > 0 &&
      ctx.extension_elem_bytes != ComputeExpectedExtensionElementBytesOrThrow(ctx)) {
    NTL::LogicError("ValidateReadWidthsOrThrow: extension element byte width mismatch");
  }
}

NTL::ZZ ReadCoeffFixed(ByteReader &reader, std::uint64_t coeff_bytes,
                       const char *what) {
  std::vector<unsigned char> tmp(U64ToSizeOrThrow(coeff_bytes, what), 0);
  reader.ReadBytes(tmp.data(), tmp.size(), what);
  return NTL::ZZFromBytes(tmp.data(), static_cast<long>(tmp.size()));
}

FieldElement ReadFieldElementFixed(ByteReader &reader,
                                   const FixedProofEncodingContext &ctx) {
  ValidateReadWidthsOrThrow(ctx);
  NTL::ZZ_pX poly;
  for (long i = 0; i < ctx.base_ext_degree; ++i) {
    const NTL::ZZ coeff = ReadCoeffFixed(
        reader, ctx.coeff_bytes,
        "ReadFieldElementFixed: insufficient bytes for coefficient");
    const NTL::ZZ reduced = coeff % NTL::ZZ_p::modulus();
    NTL::SetCoeff(poly, i, NTL::conv<NTL::ZZ_p>(reduced));
  }
  FieldElement value;
  NTL::conv(value, poly);
  return value;
}

NTL::ZZ_pEX ReadExtensionElementFixed(ByteReader &reader,
                                      const FixedProofEncodingContext &ctx) {
  ValidateReadWidthsOrThrow(ctx);
  if (ctx.challenge_ext_degree <= 0) {
    NTL::LogicError(
        "ReadExtensionElementFixed: challenge_ext_degree must be > 0");
  }
  NTL::ZZ_pEX poly;
  for (long i = 0; i < ctx.challenge_ext_degree; ++i) {
    NTL::SetCoeff(poly, i, ReadFieldElementFixed(reader, ctx));
  }
  return poly;
}

Digest ReadDigestFixed(ByteReader &reader, const FixedProofEncodingContext &ctx) {
  Digest digest{};
  if (ctx.hash_bytes != digest.size()) {
    NTL::LogicError("ReadDigestFixed: digest width mismatch");
  }
  reader.ReadBytes(digest.data(), digest.size(),
                   "ReadDigestFixed: insufficient bytes for digest");
  return digest;
}

QuadraticPoly ReadQuadraticPolyFixed(ByteReader &reader,
                                     const FixedProofEncodingContext &ctx) {
  QuadraticPoly poly;
  poly.a0 = ReadFieldElementFixed(reader, ctx);
  poly.a1 = ReadFieldElementFixed(reader, ctx);
  poly.a2 = ReadFieldElementFixed(reader, ctx);
  return poly;
}

ExtensionQuadraticPoly ReadExtensionQuadraticPolyFixed(
    ByteReader &reader, const FixedProofEncodingContext &ctx) {
  ExtensionQuadraticPoly poly;
  poly.a0 = ReadExtensionElementFixed(reader, ctx);
  poly.a1 = ReadExtensionElementFixed(reader, ctx);
  poly.a2 = ReadExtensionElementFixed(reader, ctx);
  return poly;
}

Oracle ReadOracleFixed(ByteReader &reader, const FixedProofEncodingContext &ctx,
                       const char *what) {
  const std::uint64_t count_u64 = reader.ReadU64(what);
  const long count = U64ToLongOrThrow(count_u64, what);
  Oracle oracle;
  oracle.SetLength(count);
  for (long i = 0; i < count; ++i) {
    oracle[i] = ReadFieldElementFixed(reader, ctx);
  }
  return oracle;
}

MerkleMultiproof ReadMerkleMultiproofFixed(ByteReader &reader,
                                           const FixedProofEncodingContext &ctx,
                                           const char *what) {
  MerkleMultiproof proof;
  const std::uint64_t value_count_u64 = reader.ReadU64(what);
  const long value_count = U64ToLongOrThrow(value_count_u64, what);
  proof.values.SetLength(value_count);
  for (long i = 0; i < value_count; ++i) {
    proof.values[i] = ReadFieldElementFixed(reader, ctx);
  }

  const std::uint64_t sibling_count_u64 = reader.ReadU64(what);
  const std::size_t sibling_count = U64ToSizeOrThrow(sibling_count_u64, what);
  proof.sibling_hashes.resize(sibling_count);
  for (Digest &digest : proof.sibling_hashes) {
    digest = ReadDigestFixed(reader, ctx);
  }
  return proof;
}

ExtensionMerkleMultiproof ReadExtensionMerkleMultiproofFixed(
    ByteReader &reader, const FixedProofEncodingContext &ctx, const char *what) {
  ExtensionMerkleMultiproof proof;
  const std::uint64_t value_count_u64 = reader.ReadU64(what);
  const std::size_t value_count = U64ToSizeOrThrow(value_count_u64, what);
  proof.values.resize(value_count);
  for (NTL::ZZ_pEX &value : proof.values) {
    value = ReadExtensionElementFixed(reader, ctx);
  }

  const std::uint64_t sibling_count_u64 = reader.ReadU64(what);
  const std::size_t sibling_count = U64ToSizeOrThrow(sibling_count_u64, what);
  proof.sibling_hashes.resize(sibling_count);
  for (Digest &digest : proof.sibling_hashes) {
    digest = ReadDigestFixed(reader, ctx);
  }
  return proof;
}

IOPPMerkleCommitments ReadCommitmentsFixed(ByteReader &reader,
                                           const FixedProofEncodingContext &ctx) {
  IOPPMerkleCommitments commitments;
  const std::uint64_t count_u64 = reader.ReadU64(
      "ReadCommitmentsFixed: insufficient bytes for commitment count");
  const std::size_t count = U64ToSizeOrThrow(
      count_u64, "ReadCommitmentsFixed: commitment count overflow");
  commitments.roots_by_level.resize(count);
  for (MerkleRoot &root : commitments.roots_by_level) {
    root = ReadDigestFixed(reader, ctx);
  }
  return commitments;
}

BaseFoldPCSExtensionProofData ReadExtensionProofDataFixed(
    ByteReader &reader, const FixedProofEncodingContext &ctx) {
  BaseFoldPCSExtensionProofData extension;
  extension.has_extension_payload = true;

  const std::uint64_t roots_count_u64 = reader.ReadU64(
      "ReadExtensionProofDataFixed: insufficient bytes for extension root count");
  const std::size_t roots_count = U64ToSizeOrThrow(
      roots_count_u64,
      "ReadExtensionProofDataFixed: extension root count overflow");
  extension.roots_by_level.resize(roots_count);
  for (MerkleRoot &root : extension.roots_by_level) {
    root = ReadDigestFixed(reader, ctx);
  }

  const std::uint64_t h_count_u64 = reader.ReadU64(
      "ReadExtensionProofDataFixed: insufficient bytes for extension h count");
  const std::size_t h_count = U64ToSizeOrThrow(
      h_count_u64, "ReadExtensionProofDataFixed: extension h count overflow");
  extension.h_by_level.resize(h_count);
  for (ExtensionQuadraticPoly &h : extension.h_by_level) {
    h = ReadExtensionQuadraticPolyFixed(reader, ctx);
  }

  const std::uint64_t msg0_count_u64 = reader.ReadU64(
      "ReadExtensionProofDataFixed: insufficient bytes for extension msg0 count");
  const std::size_t msg0_count = U64ToSizeOrThrow(
      msg0_count_u64,
      "ReadExtensionProofDataFixed: extension msg0 count overflow");
  extension.msg0_coeffs.resize(msg0_count);
  for (NTL::ZZ_pEX &coeff : extension.msg0_coeffs) {
    coeff = ReadExtensionElementFixed(reader, ctx);
  }

  const std::uint64_t pi0_count_u64 = reader.ReadU64(
      "ReadExtensionProofDataFixed: insufficient bytes for extension pi0 count");
  const std::size_t pi0_count = U64ToSizeOrThrow(
      pi0_count_u64,
      "ReadExtensionProofDataFixed: extension pi0 count overflow");
  extension.pi0_codeword.resize(pi0_count);
  for (NTL::ZZ_pEX &value : extension.pi0_codeword) {
    value = ReadExtensionElementFixed(reader, ctx);
  }

  extension.base_top_query_multiproof = ReadMerkleMultiproofFixed(
      reader, ctx,
      "ReadExtensionProofDataFixed: insufficient bytes for base top multiproof");

  const std::uint64_t multiproof_count_u64 = reader.ReadU64(
      "ReadExtensionProofDataFixed: insufficient bytes for extension multiproof count");
  const std::size_t multiproof_count = U64ToSizeOrThrow(
      multiproof_count_u64,
      "ReadExtensionProofDataFixed: extension multiproof count overflow");
  extension.query_multiproofs.resize(multiproof_count);
  for (ExtensionMerkleMultiproof &multiproof : extension.query_multiproofs) {
    multiproof = ReadExtensionMerkleMultiproofFixed(
        reader, ctx,
        "ReadExtensionProofDataFixed: insufficient bytes for extension multiproof");
  }

  // Current fixed-width proof contract omits verifier-reconstructable metadata.
  extension.r_by_level.clear();
  return extension;
}

FixedProofEncodingContext BuildReadContextOrThrow(
    const FixedProofEncodingOptions &options) {
  FixedProofEncodingContext ctx;
  ctx.include_version_byte = options.include_version_byte;
  ctx.base_ext_degree = NTL::ZZ_pE::degree();
  if (ctx.base_ext_degree <= 0) {
    NTL::LogicError(
        "BuildReadContextOrThrow: invalid base extension degree");
  }
  ctx.challenge_ext_degree = options.challenge_ext_degree;
  ctx.coeff_bytes =
      fixed_proof_serialize_detail::ComputeFixedCoeffBytesOrThrow();
  ctx.field_elem_bytes = ComputeExpectedFieldElementBytesOrThrow(ctx);
  if (ctx.challenge_ext_degree > 0) {
    ctx.extension_elem_bytes = ComputeExpectedExtensionElementBytesOrThrow(ctx);
  }
  return ctx;
}

}  // namespace

BaseFoldPCSEvalProof DeserializeBaseFoldPCSEvalProofFixedBytes(
    const Bytes &bytes, const FixedProofEncodingOptions &options) {
  return DeserializeBaseFoldPCSEvalProofFixedBytes(bytes.data(), bytes.size(),
                                                   options);
}

BaseFoldPCSEvalProof DeserializeBaseFoldPCSEvalProofFixedBytes(
    const Byte *data, std::size_t size, const FixedProofEncodingOptions &options) {
  if (size > 0 && data == nullptr) {
    NTL::LogicError(
        "DeserializeBaseFoldPCSEvalProofFixedBytes: non-empty input requires data");
  }

  const FixedProofEncodingContext ctx = BuildReadContextOrThrow(options);
  ByteReader reader(data, size);

  if (ctx.include_version_byte) {
    const std::uint8_t version = reader.ReadU8(
        "DeserializeBaseFoldPCSEvalProofFixedBytes: missing version byte");
    if (version != ctx.format_version) {
      NTL::LogicError(
          "DeserializeBaseFoldPCSEvalProofFixedBytes: unsupported format version");
    }
  }

  BaseFoldPCSEvalProof proof;
  proof.commitments = ReadCommitmentsFixed(reader, ctx);

  const std::uint64_t h_count_u64 = reader.ReadU64(
      "DeserializeBaseFoldPCSEvalProofFixedBytes: missing h count");
  const std::size_t h_count = U64ToSizeOrThrow(
      h_count_u64,
      "DeserializeBaseFoldPCSEvalProofFixedBytes: h count overflow");
  proof.h_by_level.resize(h_count);
  for (QuadraticPoly &h : proof.h_by_level) {
    h = ReadQuadraticPolyFixed(reader, ctx);
  }

  proof.pi0_codeword = ReadOracleFixed(
      reader, ctx,
      "DeserializeBaseFoldPCSEvalProofFixedBytes: missing pi0 codeword");

  const std::uint64_t multiproof_count_u64 = reader.ReadU64(
      "DeserializeBaseFoldPCSEvalProofFixedBytes: missing multiproof count");
  const std::size_t multiproof_count = U64ToSizeOrThrow(
      multiproof_count_u64,
      "DeserializeBaseFoldPCSEvalProofFixedBytes: multiproof count overflow");
  proof.query_multiproofs.resize(multiproof_count);
  for (MerkleMultiproof &multiproof : proof.query_multiproofs) {
    multiproof = ReadMerkleMultiproofFixed(
        reader, ctx,
        "DeserializeBaseFoldPCSEvalProofFixedBytes: missing multiproof payload");
  }

  const std::uint8_t has_extension = reader.ReadU8(
      "DeserializeBaseFoldPCSEvalProofFixedBytes: missing extension flag");
  if (has_extension > 1) {
    NTL::LogicError(
        "DeserializeBaseFoldPCSEvalProofFixedBytes: invalid extension flag");
  }
  proof.extension.has_extension_payload = (has_extension != 0);
  if (proof.extension.has_extension_payload) {
    if (ctx.challenge_ext_degree <= 0) {
      NTL::LogicError(
          "DeserializeBaseFoldPCSEvalProofFixedBytes: extension proof requires challenge_ext_degree > 0");
    }
    proof.extension = ReadExtensionProofDataFixed(reader, ctx);
  }

  if (reader.remaining() != 0) {
    NTL::LogicError(
        "DeserializeBaseFoldPCSEvalProofFixedBytes: trailing bytes remain");
  }
  return proof;
}

}  // namespace basefold
