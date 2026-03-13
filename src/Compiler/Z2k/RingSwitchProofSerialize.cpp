#include "Compiler/Z2k/RingSwitchProofSerialize.hpp"

#include "Compiler/Z2k/ProofSerializeCommon.hpp"

using NTL::LogicError;

namespace basefold {
namespace {

template <typename OuterProofLike>
void ValidateOuterProofShapeOrThrow(const RingSwitchPCSParams &params,
                                    const OuterProofLike &proof,
                                    const char *func_name) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  const long expected_s_count =
      static_cast<long>(params.beta_basis.basis.size());
  z2k_fixed_proof_serialize_detail::ValidateOuterProofShapeOrThrow(
      expected_s_count, proof.s_by_u.size(), "s_by_u", "beta basis dimension",
      params.ell_prime, proof.h_by_level.size(), func_name);
}

template <typename Sink, typename OuterProofLike>
void SerializeOuterProofToSink(Sink &sink, const RingSwitchPCSParams &params,
                               const OuterProofLike &proof,
                               const FixedProofEncodingContext &ctx) {
  (void)params;
  z2k_fixed_proof_serialize_detail::SerializeOuterProofBodyToSink(
      sink, proof.s_by_u, "s_by_u", proof.h_by_level, proof.t_star, ctx);
}

}  // namespace

Bytes SerializeRingSwitchPCSOuterProofFixedBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSOuterEvalProof &proof,
    const RingSwitchProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "SerializeRingSwitchPCSOuterProofFixedBytes");
  const FixedProofEncodingContext ctx =
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildRingSwitchOuterEncodingContext");
  ByteBufferSink sink;
  SerializeOuterProofToSink(sink, params, proof, ctx);
  return sink.bytes();
}

Bytes SerializeRingSwitchPCSOuterProofFixedBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "SerializeRingSwitchPCSOuterProofFixedBytes");
  const FixedProofEncodingContext ctx =
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildRingSwitchOuterEncodingContext");
  ByteBufferSink sink;
  SerializeOuterProofToSink(sink, params, proof, ctx);
  return sink.bytes();
}

Bytes SerializeRingSwitchPCSEvalProofFixedBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "SerializeRingSwitchPCSEvalProofFixedBytes");
  const FixedProofEncodingContext ctx =
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildRingSwitchOuterEncodingContext");
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

std::uint64_t RingSwitchPCSOuterProofSizeBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSOuterEvalProof &proof,
    const RingSwitchProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "RingSwitchPCSOuterProofSizeBytes");
  const FixedProofEncodingContext ctx =
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildRingSwitchOuterEncodingContext");
  CountingSink sink;
  SerializeOuterProofToSink(sink, params, proof, ctx);
  return sink.bytes_written();
}

std::uint64_t RingSwitchPCSOuterProofSizeBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options) {
  ValidateOuterProofShapeOrThrow(params, proof,
                                 "RingSwitchPCSOuterProofSizeBytes");
  const FixedProofEncodingContext ctx =
      z2k_fixed_proof_serialize_detail::BuildOuterEncodingContextOrThrow(
          options.include_version_byte,
          "BuildRingSwitchOuterEncodingContext");
  CountingSink sink;
  SerializeOuterProofToSink(sink, params, proof, ctx);
  return sink.bytes_written();
}

double RingSwitchPCSOuterProofSizeKB(
    const RingSwitchPCSParams &params, const RingSwitchPCSOuterEvalProof &proof,
    const RingSwitchProofEncodingOptions &options) {
  return static_cast<double>(
             RingSwitchPCSOuterProofSizeBytes(params, proof, options)) /
         1024.0;
}

double RingSwitchPCSOuterProofSizeKB(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options) {
  return static_cast<double>(
             RingSwitchPCSOuterProofSizeBytes(params, proof, options)) /
         1024.0;
}

std::uint64_t RingSwitchPCSEvalProofSizeBytes(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options) {
  const std::uint64_t outer_bytes =
      RingSwitchPCSOuterProofSizeBytes(params, proof, options);
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
          "RingSwitchPCSEvalProofSizeBytes: byte count overflow"),
      backend_bytes,
      "RingSwitchPCSEvalProofSizeBytes: byte count overflow");
}

double RingSwitchPCSEvalProofSizeKB(
    const RingSwitchPCSParams &params, const RingSwitchPCSEvalProof &proof,
    const RingSwitchProofEncodingOptions &options) {
  return static_cast<double>(
             RingSwitchPCSEvalProofSizeBytes(params, proof, options)) /
         1024.0;
}

}  // namespace basefold
