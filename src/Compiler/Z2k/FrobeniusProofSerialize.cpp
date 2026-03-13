#include "Compiler/Z2k/FrobeniusProofSerialize.hpp"

#include "Compiler/Z2k/ProofSerializeCommon.hpp"

using NTL::LogicError;

namespace basefold {
namespace {

template <typename OuterProofLike>
void ValidateOuterProofShapeOrThrow(const FrobeniusPCSParams &params,
                                    const OuterProofLike &proof,
                                    const char *func_name) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const long expected_s_count = static_cast<long>(
      params.basis_data.normal_basis.beta.size());
  z2k_fixed_proof_serialize_detail::ValidateOuterProofShapeOrThrow(
      expected_s_count, proof.s_by_i.size(), "s_by_i",
      "Frobenius basis dimension", params.ell_prime, proof.h_by_level.size(),
      func_name);
}

template <typename Sink, typename OuterProofLike>
void SerializeOuterProofToSink(Sink &sink, const FrobeniusPCSParams &params,
                               const OuterProofLike &proof,
                               const FixedProofEncodingContext &ctx) {
  (void)params;
  z2k_fixed_proof_serialize_detail::SerializeOuterProofBodyToSink(
      sink, proof.s_by_i, "s_by_i", proof.h_by_level, proof.t_star, ctx);
}

}  // namespace

Bytes SerializeFrobeniusPCSOuterProofFixedBytes(
    const FrobeniusPCSParams &params, const FrobeniusPCSOuterEvalProof &proof,
    const FrobeniusProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "SerializeFrobeniusPCSOuterProofFixedBytes");
  const FixedProofEncodingContext ctx =
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildFrobeniusOuterEncodingContext");
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
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildFrobeniusOuterEncodingContext");
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
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildFrobeniusOuterEncodingContext");
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
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildFrobeniusOuterEncodingContext");
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
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildFrobeniusOuterEncodingContext");
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
  return z2k_fixed_proof_serialize_detail::AddU64OrThrow(
      z2k_fixed_proof_serialize_detail::AddU64OrThrow(
          outer_bytes, static_cast<std::uint64_t>(8),
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
