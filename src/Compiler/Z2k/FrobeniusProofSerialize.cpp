#include "Compiler/Z2k/FrobeniusProofSerialize.hpp"

#include <NTL/ZZ_pE.h>

#include <limits>
#include <string>

#include "PCS/BaseFold/ProofSerialize.hpp"

using NTL::LogicError;

namespace basefold {
namespace {

FixedProofEncodingContext BuildFrobeniusOuterEncodingContext(
    const FrobeniusProofEncodingOptions &options) {
  FixedProofEncodingContext ctx;
  ctx.include_version_byte = options.include_version_byte;
  ctx.base_ext_degree = NTL::ZZ_pE::degree();
  if (ctx.base_ext_degree <= 0) {
    LogicError(
        "BuildFrobeniusOuterEncodingContext: invalid base extension degree");
  }
  ctx.coeff_bytes =
      fixed_proof_serialize_detail::ComputeFixedCoeffBytesOrThrow();
  ctx.field_elem_bytes = fixed_proof_serialize_detail::MulU64OrThrow(
      ctx.coeff_bytes,
      fixed_proof_serialize_detail::LongToU64OrThrow(
          ctx.base_ext_degree,
          "BuildFrobeniusOuterEncodingContext: base extension degree must be "
          ">= 0"),
      "BuildFrobeniusOuterEncodingContext: field element byte width overflow");
  return ctx;
}

template <typename OuterProofLike>
void ValidateOuterProofShapeOrThrow(const FrobeniusPCSParams &params,
                                    const OuterProofLike &proof,
                                    const char *func_name) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const long expected_s_count = static_cast<long>(
      params.basis_data.normal_basis.beta.size());
  if (static_cast<long>(proof.s_by_i.size()) != expected_s_count) {
    LogicError((std::string(func_name) +
                ": s_by_i count must equal Frobenius basis dimension")
                   .c_str());
  }
  if (static_cast<long>(proof.h_by_level.size()) != params.ell_prime) {
    LogicError((std::string(func_name) +
                ": h_by_level count must equal ell_prime")
                   .c_str());
  }
}

template <typename Sink, typename OuterProofLike>
void SerializeOuterProofToSink(Sink &sink, const FrobeniusPCSParams &params,
                               const OuterProofLike &proof,
                               const FixedProofEncodingContext &ctx) {
  (void)params;
  fixed_proof_serialize_detail::SerializeVersion(sink, ctx);

  const std::uint64_t s_count = fixed_proof_serialize_detail::SizeToU64OrThrow(
      proof.s_by_i.size(),
      "SerializeOuterProofToSink: s_by_i count overflow");
  fixed_proof_serialize_detail::SerializeVecHeader(sink, s_count);
  for (const FieldElement &s_i : proof.s_by_i) {
    fixed_proof_serialize_detail::SerializeFieldElementFixed(sink, s_i, ctx);
  }

  const std::uint64_t h_count = fixed_proof_serialize_detail::SizeToU64OrThrow(
      proof.h_by_level.size(),
      "SerializeOuterProofToSink: h_by_level count overflow");
  fixed_proof_serialize_detail::SerializeVecHeader(sink, h_count);
  for (const QuadraticPoly &h : proof.h_by_level) {
    fixed_proof_serialize_detail::SerializeQuadraticPolyFixed(sink, h, ctx);
  }

  fixed_proof_serialize_detail::SerializeFieldElementFixed(sink, proof.t_star,
                                                           ctx);
}

std::uint64_t AddU64OrThrow(std::uint64_t lhs, std::uint64_t rhs,
                            const char *what) {
  if (lhs > std::numeric_limits<std::uint64_t>::max() - rhs) {
    LogicError(what);
  }
  return lhs + rhs;
}

}  // namespace

Bytes SerializeFrobeniusPCSOuterProofFixedBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSOuterEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "SerializeFrobeniusPCSOuterProofFixedBytes");
  const FixedProofEncodingContext ctx =
      BuildFrobeniusOuterEncodingContext(options);
  ByteBufferSink sink;
  SerializeOuterProofToSink(sink, params, proof, ctx);
  return sink.bytes();
}

Bytes SerializeFrobeniusPCSOuterProofFixedBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "SerializeFrobeniusPCSOuterProofFixedBytes");
  const FixedProofEncodingContext ctx =
      BuildFrobeniusOuterEncodingContext(options);
  ByteBufferSink sink;
  SerializeOuterProofToSink(sink, params, proof, ctx);
  return sink.bytes();
}

Bytes SerializeFrobeniusPCSEvalProofFixedBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "SerializeFrobeniusPCSEvalProofFixedBytes");
  const FixedProofEncodingContext ctx =
      BuildFrobeniusOuterEncodingContext(options);
  ByteBufferSink sink;
  SerializeOuterProofToSink(sink, params, proof, ctx);

  const Bytes backend_bytes = Z2kPCSBackendSerializeEvalProof(
      params.backend, proof.backend_proof, options.backend_proof_options);
  sink.WriteU64(static_cast<std::uint64_t>(backend_bytes.size()));
  if (!backend_bytes.empty()) {
    sink.WriteBytes(backend_bytes.data(),
                    static_cast<std::uint64_t>(backend_bytes.size()));
  }
  return sink.bytes();
}

std::uint64_t FrobeniusPCSOuterProofSizeBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSOuterEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "FrobeniusPCSOuterProofSizeBytes");
  const FixedProofEncodingContext ctx =
      BuildFrobeniusOuterEncodingContext(options);
  CountingSink sink;
  SerializeOuterProofToSink(sink, params, proof, ctx);
  return sink.bytes_written();
}

std::uint64_t FrobeniusPCSOuterProofSizeBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "FrobeniusPCSOuterProofSizeBytes");
  const FixedProofEncodingContext ctx =
      BuildFrobeniusOuterEncodingContext(options);
  CountingSink sink;
  SerializeOuterProofToSink(sink, params, proof, ctx);
  return sink.bytes_written();
}

double FrobeniusPCSOuterProofSizeKB(
    const FrobeniusPCSParams &params, const FrobeniusPCSOuterEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  return static_cast<double>(
             FrobeniusPCSOuterProofSizeBytes(params, proof, options)) /
         1024.0;
}

double FrobeniusPCSOuterProofSizeKB(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  return static_cast<double>(
             FrobeniusPCSOuterProofSizeBytes(params, proof, options)) /
         1024.0;
}

std::uint64_t FrobeniusPCSEvalProofSizeBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  const std::uint64_t outer_bytes =
      FrobeniusPCSOuterProofSizeBytes(params, proof, options);
  Z2kPCSBackendProofSizeOptions backend_options;
  backend_options.include_version_byte =
      options.backend_proof_options.include_version_byte;
  backend_options.challenge_ext_degree =
      options.backend_proof_options.challenge_ext_degree;
  const std::uint64_t backend_bytes = Z2kPCSBackendEvalProofSizeBytes(
      params.backend, proof.backend_proof, backend_options);
  return AddU64OrThrow(
      AddU64OrThrow(outer_bytes, static_cast<std::uint64_t>(8),
                    "FrobeniusPCSEvalProofSizeBytes: byte count overflow"),
      backend_bytes,
      "FrobeniusPCSEvalProofSizeBytes: byte count overflow");
}

double FrobeniusPCSEvalProofSizeKB(
    const FrobeniusPCSParams &params, const FrobeniusPCSEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  return static_cast<double>(
             FrobeniusPCSEvalProofSizeBytes(params, proof, options)) /
         1024.0;
}

}  // namespace basefold
