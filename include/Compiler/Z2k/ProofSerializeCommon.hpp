#ifndef BASEFOLD_Z2K_PROOF_SERIALIZE_COMMON_HPP_
#define BASEFOLD_Z2K_PROOF_SERIALIZE_COMMON_HPP_

#include <NTL/ZZ_pE.h>

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "PCS/BaseFold/ProofSerialize.hpp"

namespace basefold {
namespace z2k_fixed_proof_serialize_detail {

inline FixedProofEncodingContext BuildOuterEncodingContextOrThrow(
    bool include_version_byte, const char *func_name) {
  FixedProofEncodingContext ctx;
  ctx.include_version_byte = include_version_byte;
  ctx.base_ext_degree = NTL::ZZ_pE::degree();
  if (ctx.base_ext_degree <= 0) {
    NTL::LogicError((std::string(func_name) +
                     ": invalid base extension degree")
                        .c_str());
  }
  ctx.coeff_bytes =
      fixed_proof_serialize_detail::ComputeFixedCoeffBytesOrThrow();
  ctx.field_elem_bytes = fixed_proof_serialize_detail::MulU64OrThrow(
      ctx.coeff_bytes,
      fixed_proof_serialize_detail::LongToU64OrThrow(
          ctx.base_ext_degree,
          (std::string(func_name) +
           ": base extension degree must be >= 0")
              .c_str()),
      (std::string(func_name) + ": field element byte width overflow")
          .c_str());
  return ctx;
}

inline void ValidateOuterProofShapeOrThrow(
    long expected_s_count, std::size_t actual_s_count, const char *s_label,
    const char *s_expected_desc, long expected_h_count,
    std::size_t actual_h_count, const char *func_name) {
  if (static_cast<long>(actual_s_count) != expected_s_count) {
    NTL::LogicError((std::string(func_name) + ": " + s_label +
                     " count must equal " + s_expected_desc)
                        .c_str());
  }
  if (static_cast<long>(actual_h_count) != expected_h_count) {
    NTL::LogicError((std::string(func_name) +
                     ": h_by_level count must equal ell_prime")
                        .c_str());
  }
}

template <typename Sink>
inline void SerializeOuterProofBodyToSink(
    Sink &sink, const std::vector<FieldElement> &s_values,
    const char *s_label, const std::vector<QuadraticPoly> &h_by_level,
    const FieldElement &t_star, const FixedProofEncodingContext &ctx) {
  fixed_proof_serialize_detail::SerializeVersion(sink, ctx);

  const std::uint64_t s_count = fixed_proof_serialize_detail::SizeToU64OrThrow(
      s_values.size(),
      (std::string("SerializeOuterProofBodyToSink: ") + s_label +
       " count overflow")
          .c_str());
  fixed_proof_serialize_detail::SerializeVecHeader(sink, s_count);
  for (const FieldElement &s_value : s_values) {
    fixed_proof_serialize_detail::SerializeFieldElementFixed(sink, s_value,
                                                             ctx);
  }

  const std::uint64_t h_count = fixed_proof_serialize_detail::SizeToU64OrThrow(
      h_by_level.size(),
      "SerializeOuterProofBodyToSink: h_by_level count overflow");
  fixed_proof_serialize_detail::SerializeVecHeader(sink, h_count);
  for (const QuadraticPoly &h : h_by_level) {
    fixed_proof_serialize_detail::SerializeQuadraticPolyFixed(sink, h, ctx);
  }

  fixed_proof_serialize_detail::SerializeFieldElementFixed(sink, t_star, ctx);
}

inline std::uint64_t AddU64OrThrow(std::uint64_t lhs, std::uint64_t rhs,
                                   const char *what) {
  if (lhs > std::numeric_limits<std::uint64_t>::max() - rhs) {
    NTL::LogicError(what);
  }
  return lhs + rhs;
}

}  // namespace z2k_fixed_proof_serialize_detail
}  // namespace basefold

#endif  // BASEFOLD_Z2K_PROOF_SERIALIZE_COMMON_HPP_
