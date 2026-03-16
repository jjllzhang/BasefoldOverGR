#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "PCS/Common/Hash.hpp"
#include "PCS/Common/Multilinear.hpp"
#include "PCS/Common/Sumcheck.hpp"
#include "Compiler/Z2k/BaseFoldBackendAdapter.hpp"
#include "Compiler/Z2k/RingSwitchPCS.hpp"
#include "Compiler/Z2k/RingSwitchProofSerialize.hpp"
#include "GaloisRing/Basis.hpp"
#include "GaloisRing/utils.hpp"
#include "tests/test_common.hpp"

using NTL::conv;
using NTL::clear;
using NTL::mat_ZZ_pE;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using std::cerr;
using std::cout;
using std::exception;
using std::function;
using std::string;

int g_test_failure_count = 0;

namespace {

basefold::FoldableCodeParams BuildParamsGR42(const ZZ &p, const ZZ_pE &alpha) {
  const long c = 2;
  const long k0 = 1;
  const long d = 2;
  const long n0 = c * k0;

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = p;
  params.zeta = alpha;

  mat_ZZ_pE G0;
  G0.SetDims(k0, n0);
  const ZZ_pE one = testutil::ConstZZpE(1);
  G0[0][0] = one;
  G0[0][1] = alpha;
  params.G0 = G0;

  params.diag_T.resize(static_cast<std::size_t>(d));
  params.diag_T[0].SetLength(2);
  params.diag_T[0][0] = one;
  params.diag_T[0][1] = alpha;

  params.diag_T[1].SetLength(4);
  params.diag_T[1][0] = alpha + one;
  params.diag_T[1][1] = one;
  params.diag_T[1][2] = alpha;
  params.diag_T[1][3] = alpha + one;
  return params;
}

basefold::FoldableCodeParams BuildParamsGR42Variant(const ZZ &p,
                                                    const ZZ_pE &alpha) {
  basefold::FoldableCodeParams params = BuildParamsGR42(p, alpha);
  const ZZ_pE one = testutil::ConstZZpE(1);
  params.zeta = alpha + one;
  params.G0[0][1] = one;
  params.diag_T[0][1] = one;
  params.diag_T[1][0] = one;
  params.diag_T[1][2] = alpha + one;
  return params;
}

basefold::FoldableCodeParams BuildParamsGR42D0(const ZZ &p, const ZZ_pE &alpha) {
  const long c = 2;
  const long k0 = 1;
  const long d = 0;
  const long n0 = c * k0;

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = p;
  params.zeta = alpha;

  mat_ZZ_pE G0;
  G0.SetDims(k0, n0);
  const ZZ_pE one = testutil::ConstZZpE(1);
  G0[0][0] = one;
  G0[0][1] = alpha;
  params.G0 = G0;
  return params;
}

enum class RingSwitchBackendVariant {
  kDefault,
  kVariant,
  kD0,
};

enum class BasisTransformFamily {
  kUnitriangular,
  kDense,
};

struct RingSwitchTestContextSpec {
  const char *name = nullptr;
  ZZ p;
  ZZ modulus;
  std::vector<long> extension_coeffs;
  std::vector<long> alpha_coeffs;
  RingSwitchBackendVariant backend_variant = RingSwitchBackendVariant::kDefault;
};

struct RingSwitchBasisCase {
  std::string label;
  std::vector<ZZ_pE> basis;
  std::vector<ZZ_pE> dual_basis;
};

struct RingSwitchEvalCaseData {
  vec_ZZ_pE t_table;
  std::vector<ZZ_pE> z;
  long num_queries = 2;
};

RingSwitchTestContextSpec MakeGR42ContextSpec() {
  RingSwitchTestContextSpec spec;
  spec.name = "GR(4,2)";
  spec.p = to_ZZ(2);
  spec.modulus = to_ZZ(4);
  spec.extension_coeffs = {1, 1, 1};
  spec.alpha_coeffs = {0, 1};
  spec.backend_variant = RingSwitchBackendVariant::kDefault;
  return spec;
}

RingSwitchTestContextSpec MakeGR42VariantContextSpec() {
  RingSwitchTestContextSpec spec = MakeGR42ContextSpec();
  spec.name = "GR(4,2)-variant";
  spec.backend_variant = RingSwitchBackendVariant::kVariant;
  return spec;
}

RingSwitchTestContextSpec MakeGR44ContextSpec() {
  RingSwitchTestContextSpec spec = MakeGR42ContextSpec();
  spec.name = "GR(4,4)";
  spec.extension_coeffs = {1, 1, 0, 0, 1};
  spec.alpha_coeffs = {0, 1, 0, 0};
  return spec;
}

ZZ_pX BuildExtensionPolynomial(const RingSwitchTestContextSpec &spec) {
  ZZ_pX F;
  for (long i = 0; i < static_cast<long>(spec.extension_coeffs.size()); ++i) {
    SetCoeff(F, i, spec.extension_coeffs[static_cast<std::size_t>(i)]);
  }
  return F;
}

ZZ_pE BuildElementFromPolynomialCoords(const std::vector<long> &coeffs) {
  ZZ_pX poly;
  for (long i = 0; i < static_cast<long>(coeffs.size()); ++i) {
    SetCoeff(poly, i, coeffs[static_cast<std::size_t>(i)]);
  }
  ZZ_pE out;
  conv(out, poly);
  return out;
}

basefold::FoldableCodeParams BuildBackendParamsFromSpec(
    const RingSwitchTestContextSpec &spec, const ZZ_pE &alpha) {
  switch (spec.backend_variant) {
    case RingSwitchBackendVariant::kDefault:
      return BuildParamsGR42(spec.p, alpha);
    case RingSwitchBackendVariant::kVariant:
      return BuildParamsGR42Variant(spec.p, alpha);
    case RingSwitchBackendVariant::kD0:
      return BuildParamsGR42D0(spec.p, alpha);
  }
  CHECK_MSG(false, "BuildBackendParamsFromSpec: unknown backend variant");
  return basefold::FoldableCodeParams{};
}

class ScopedRingSwitchTestContext {
 public:
  explicit ScopedRingSwitchTestContext(const RingSwitchTestContextSpec &spec)
      : spec_(spec),
        mod_push_(spec.modulus),
        F_(BuildExtensionPolynomial(spec_)),
        ext_push_(F_),
        alpha_(BuildElementFromPolynomialCoords(spec_.alpha_coeffs)) {}

  const RingSwitchTestContextSpec &spec() const { return spec_; }
  const ZZ_pX &F() const { return F_; }
  const ZZ_pE &alpha() const { return alpha_; }

 private:
  RingSwitchTestContextSpec spec_;
  ZZ_pPush mod_push_;
  ZZ_pX F_;
  ZZ_pEPush ext_push_;
  ZZ_pE alpha_;
};

std::vector<NTL::ZZ_p> LongVecToZZp(const std::vector<long> &coeffs) {
  std::vector<NTL::ZZ_p> out(coeffs.size());
  for (long i = 0; i < static_cast<long>(coeffs.size()); ++i) {
    conv(out[static_cast<std::size_t>(i)], coeffs[static_cast<std::size_t>(i)]);
  }
  return out;
}

long CurrentBaseModulusLong() { return NTL::conv<long>(ZZ_p::modulus()); }

long NormalizeBaseCoeff(long value) {
  const long modulus = CurrentBaseModulusLong();
  CHECK_MSG(modulus > 1, "NormalizeBaseCoeff: invalid base modulus");
  value %= modulus;
  if (value < 0) {
    value += modulus;
  }
  return value;
}

std::vector<std::vector<long>> IdentityBaseCoeffMatrix(long dimension) {
  std::vector<std::vector<long>> matrix(static_cast<std::size_t>(dimension),
                                        std::vector<long>(
                                            static_cast<std::size_t>(dimension),
                                            0));
  for (long i = 0; i < dimension; ++i) {
    matrix[static_cast<std::size_t>(i)][static_cast<std::size_t>(i)] = 1;
  }
  return matrix;
}

std::vector<std::vector<long>> MultiplyBaseCoeffMatrices(
    const std::vector<std::vector<long>> &lhs,
    const std::vector<std::vector<long>> &rhs) {
  const long dimension = static_cast<long>(lhs.size());
  CHECK_EQ(dimension, static_cast<long>(rhs.size()));
  std::vector<std::vector<long>> out(static_cast<std::size_t>(dimension),
                                     std::vector<long>(
                                         static_cast<std::size_t>(dimension),
                                         0));
  for (long i = 0; i < dimension; ++i) {
    CHECK_EQ(static_cast<long>(lhs[static_cast<std::size_t>(i)].size()),
             dimension);
    CHECK_EQ(static_cast<long>(rhs[static_cast<std::size_t>(i)].size()),
             dimension);
    for (long k = 0; k < dimension; ++k) {
      const long lhs_ik =
          lhs[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)];
      for (long j = 0; j < dimension; ++j) {
        out[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
            NormalizeBaseCoeff(out[static_cast<std::size_t>(i)]
                                   [static_cast<std::size_t>(j)] +
                               lhs_ik * rhs[static_cast<std::size_t>(k)]
                                           [static_cast<std::size_t>(j)]);
      }
    }
  }
  return out;
}

long NonZeroBaseCoeff(long seed_offset) {
  const long modulus = CurrentBaseModulusLong();
  CHECK_MSG(modulus > 1, "NonZeroBaseCoeff: invalid base modulus");
  return 1L + (seed_offset % (modulus - 1L));
}

std::vector<std::vector<long>> BuildUpperUnitriangularMatrix(long dimension,
                                                             long seed) {
  std::vector<std::vector<long>> matrix = IdentityBaseCoeffMatrix(dimension);
  for (long i = 0; i < dimension; ++i) {
    for (long j = i + 1; j < dimension; ++j) {
      matrix[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
          NonZeroBaseCoeff(seed + i + j);
    }
  }
  return matrix;
}

std::vector<std::vector<long>> BuildLowerUnitriangularMatrix(long dimension,
                                                             long seed) {
  std::vector<std::vector<long>> matrix = IdentityBaseCoeffMatrix(dimension);
  for (long i = 1; i < dimension; ++i) {
    for (long j = 0; j < i; ++j) {
      matrix[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
          NonZeroBaseCoeff(seed + i + j);
    }
  }
  return matrix;
}

std::vector<std::vector<long>> BuildDeterministicBasisMatrix(
    long dimension, long seed, BasisTransformFamily family) {
  switch (family) {
    case BasisTransformFamily::kUnitriangular:
      return BuildUpperUnitriangularMatrix(dimension, seed);
    case BasisTransformFamily::kDense:
      return MultiplyBaseCoeffMatrices(
          BuildLowerUnitriangularMatrix(dimension, seed),
          BuildUpperUnitriangularMatrix(dimension, seed + 1));
  }
  CHECK_MSG(false, "BuildDeterministicBasisMatrix: unknown family");
  return {};
}

std::vector<ZZ_pE> ApplyBasisChangeMatrixOrThrow(
    const std::vector<ZZ_pE> &source_basis,
    const std::vector<std::vector<long>> &matrix, const char *func_name) {
  const long dimension = static_cast<long>(source_basis.size());
  CHECK_EQ(static_cast<long>(matrix.size()), dimension);
  std::vector<ZZ_pE> out(static_cast<std::size_t>(dimension));
  for (long row = 0; row < dimension; ++row) {
    const std::vector<long> &coords = matrix[static_cast<std::size_t>(row)];
    CHECK_EQ(static_cast<long>(coords.size()), dimension);
    out[static_cast<std::size_t>(row)] = ComposeFromBasisCoordsOrThrow(
        source_basis, LongVecToZZp(coords), func_name);
  }
  return out;
}

RingSwitchBasisCase BuildPolynomialBasisCase() {
  RingSwitchBasisCase out;
  out.label = "polynomial";
  out.basis = BuildPolynomialBasis(ZZ_pE::degree());
  out.dual_basis = BuildDualBasisOrThrow(out.basis);
  return out;
}

RingSwitchBasisCase BuildTransformedBasisCase(long seed,
                                              BasisTransformFamily family) {
  const std::vector<ZZ_pE> polynomial_basis =
      BuildPolynomialBasis(ZZ_pE::degree());
  RingSwitchBasisCase out;
  out.label = (family == BasisTransformFamily::kUnitriangular)
                  ? "unitriangular"
                  : "dense";
  out.basis = ApplyBasisChangeMatrixOrThrow(
      polynomial_basis,
      BuildDeterministicBasisMatrix(ZZ_pE::degree(), seed, family),
      "BuildTransformedBasisCase");
  out.dual_basis = BuildDualBasisOrThrow(out.basis);
  return out;
}

RingSwitchBasisCase BuildSingularBasisCase(long seed,
                                           BasisTransformFamily family) {
  RingSwitchBasisCase out = BuildTransformedBasisCase(seed, family);
  if (!out.basis.empty()) {
    out.basis.back() = out.basis.front();
  }
  out.dual_basis.clear();
  out.label += "-singular";
  return out;
}

RingSwitchBasisCase BuildBrokenDualBasisCase(long seed,
                                             BasisTransformFamily family) {
  RingSwitchBasisCase out = BuildTransformedBasisCase(seed, family);
  if (!out.dual_basis.empty()) {
    out.dual_basis[0] += testutil::ConstZZpE(1);
  }
  out.label += "-broken-dual";
  return out;
}

RingSwitchBasisCase MakeBasisCase(const std::string &label,
                                  const GaloisRingBasisData &basis_data) {
  RingSwitchBasisCase out;
  out.label = label;
  out.basis = basis_data.basis;
  out.dual_basis = basis_data.dual_basis;
  return out;
}

vec_ZZ_pE BuildBaseRingCoeffVector(const std::vector<long> &coeffs);

RingSwitchEvalCaseData BuildEvalCaseData(
    const std::vector<long> &t_table_coeffs,
    const std::vector<std::vector<long>> &query_point_coeffs,
    long num_queries = 2) {
  RingSwitchEvalCaseData out;
  out.t_table = BuildBaseRingCoeffVector(t_table_coeffs);
  out.z.reserve(query_point_coeffs.size());
  for (const std::vector<long> &coeffs : query_point_coeffs) {
    out.z.push_back(BuildElementFromPolynomialCoords(coeffs));
  }
  out.num_queries = num_queries;
  return out;
}

long Pow2ForTest(long exponent) {
  CHECK_MSG(exponent >= 0, "Pow2ForTest: exponent must be non-negative");
  if (exponent < 0) {
    return 0;
  }
  return 1L << exponent;
}

vec_ZZ_pE BuildBaseRingCoeffVector(const std::vector<long> &coeffs) {
  vec_ZZ_pE out;
  out.SetLength(static_cast<long>(coeffs.size()));
  for (long i = 0; i < static_cast<long>(coeffs.size()); ++i) {
    out[i] = testutil::ConstZZpE(coeffs[static_cast<std::size_t>(i)]);
  }
  return out;
}

std::vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  std::vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] =
        testutil::ConstZZpE((index >> i) & 1L);
  }
  return point;
}

ZZ_pE PolynomialBasisElement(long degree) {
  ZZ_pX poly;
  SetCoeff(poly, degree, 1);
  ZZ_pE out;
  conv(out, poly);
  return out;
}

std::uint64_t FixedFieldElementBytesForCurrentContext() {
  const ZZ modulus_minus_one = NTL::ZZ_p::modulus() - 1;
  const long bits = NTL::NumBits(modulus_minus_one);
  CHECK_MSG(bits > 0,
            "FixedFieldElementBytesForCurrentContext: invalid modulus bit width");
  const std::uint64_t coeff_bytes = static_cast<std::uint64_t>((bits + 7) / 8);
  return coeff_bytes * static_cast<std::uint64_t>(ZZ_pE::degree());
}

void AppendU64ForTranscript(basefold::Bytes &out, std::uint64_t value) {
  for (int i = 0; i < 8; ++i) {
    out.push_back(static_cast<basefold::Byte>((value >> (8 * i)) & 0xff));
  }
}

void AppendSerializedFieldElementForTranscript(basefold::Bytes &out,
                                               const ZZ_pE &x) {
  const long r = ZZ_pE::degree();
  CHECK_MSG(r > 0, "AppendSerializedFieldElementForTranscript: invalid extension degree");
  const ZZ_pX poly = NTL::rep(x);
  AppendU64ForTranscript(out, static_cast<std::uint64_t>(r));
  for (long i = 0; i < r; ++i) {
    const ZZ coeff = NTL::rep(NTL::coeff(poly, i));
    const long n = NTL::NumBytes(coeff);
    AppendU64ForTranscript(out, static_cast<std::uint64_t>(n));
    if (n > 0) {
      const std::size_t old_size = out.size();
      out.resize(old_size + static_cast<std::size_t>(n));
      NTL::BytesFromZZ(
          reinterpret_cast<unsigned char *>(out.data() + old_size), coeff, n);
    }
  }
}

basefold::Bytes SerializeFieldElementForTranscript(const ZZ_pE &x) {
  basefold::Bytes out;
  AppendSerializedFieldElementForTranscript(out, x);
  return out;
}

basefold::Bytes TaggedHashForTranscript(basefold::Byte tag,
                                        const basefold::Bytes &state,
                                        const basefold::Bytes &payload) {
  basefold::Bytes in;
  in.push_back(tag);
  in.insert(in.end(), state.begin(), state.end());
  AppendU64ForTranscript(in, static_cast<std::uint64_t>(payload.size()));
  in.insert(in.end(), payload.begin(), payload.end());
  return basefold::HashBytes(in);
}

basefold::Bytes TaggedHashForTranscript(basefold::Byte tag,
                                        const basefold::Bytes &state,
                                        const std::string &payload) {
  basefold::Bytes bytes(payload.begin(), payload.end());
  return TaggedHashForTranscript(tag, state, bytes);
}

class TestRingSwitchTranscript {
 public:
  TestRingSwitchTranscript() {
    const std::string domain = "RingSwitchPCS/v1";
    state_ = TaggedHashForTranscript(static_cast<basefold::Byte>(0x42),
                                     basefold::Bytes{}, domain);
  }

  void AbsorbDigest(const basefold::Digest &digest) {
    const basefold::Bytes bytes(digest.begin(), digest.end());
    state_ = TaggedHashForTranscript(static_cast<basefold::Byte>(0x01), state_,
                                     bytes);
  }

  void AbsorbFieldElement(const ZZ_pE &x) {
    state_ = TaggedHashForTranscript(static_cast<basefold::Byte>(0x02), state_,
                                     SerializeFieldElementForTranscript(x));
  }

  void AbsorbQuadraticPoly(const basefold::QuadraticPoly &poly) {
    AbsorbFieldElement(poly.a0);
    AbsorbFieldElement(poly.a1);
    AbsorbFieldElement(poly.a2);
  }

  ZZ_pE ChallengeFieldElement(const std::string &label) const {
    const long r = ZZ_pE::degree();
    CHECK_MSG(r > 0, "TestRingSwitchTranscript::ChallengeFieldElement: invalid extension degree");
    const ZZ modulus = NTL::ZZ_p::modulus();
    CHECK_MSG(modulus > 1, "TestRingSwitchTranscript::ChallengeFieldElement: invalid base modulus");

    ChallengeStream stream(state_, "fe/" + label);
    ZZ_pX poly;
    clear(poly);
    for (long i = 0; i < r; ++i) {
      NTL::ZZ_p coeff;
      conv(coeff, stream.SampleZZLessThan(modulus));
      SetCoeff(poly, i, coeff);
    }
    ZZ_pE out;
    conv(out, poly);
    return out;
  }

 private:
  class ChallengeStream {
   public:
    ChallengeStream(const basefold::Bytes &state, const std::string &label)
        : state_(state), label_(label) {}

    ZZ SampleZZLessThan(const ZZ &upper_bound) {
      CHECK_MSG(upper_bound > 0, "ChallengeStream::SampleZZLessThan: upper_bound must be positive");
      if (upper_bound == 1) {
        return ZZ(0);
      }

      const ZZ upper_minus_one = upper_bound - 1;
      const long bits = NTL::NumBits(upper_minus_one);
      if (bits <= 0) {
        return ZZ(0);
      }
      const long byte_len = (bits + 7) / 8;
      const ZZ two_to_bits = ZZ(1) << bits;

      basefold::Bytes tmp(static_cast<std::size_t>(byte_len));
      while (true) {
        ReadBytes(reinterpret_cast<std::uint8_t *>(tmp.data()),
                  static_cast<std::size_t>(byte_len));
        ZZ x = NTL::ZZFromBytes(
            reinterpret_cast<const unsigned char *>(tmp.data()), byte_len);
        x %= two_to_bits;
        if (x < upper_bound) {
          return x;
        }
      }
    }

   private:
    void ReadBytes(std::uint8_t *out, std::size_t len) {
      std::size_t written = 0;
      while (written < len) {
        if (offset_ == buf_.size()) {
          buf_ = Digest(counter_++);
          offset_ = 0;
        }
        const std::size_t take =
            std::min(len - written, buf_.size() - offset_);
        std::memcpy(out + written, buf_.data() + offset_, take);
        written += take;
        offset_ += take;
      }
    }

    basefold::Bytes Digest(std::uint64_t counter) const {
      basefold::Bytes payload;
      AppendU64ForTranscript(payload, counter);
      const basefold::Bytes stream_state = TaggedHashForTranscript(
          static_cast<basefold::Byte>(0x20), state_, label_);
      return TaggedHashForTranscript(static_cast<basefold::Byte>(0x21),
                                     stream_state, payload);
    }

    basefold::Bytes state_;
    std::string label_;
    std::uint64_t counter_ = 0;
    basefold::Bytes buf_;
    std::size_t offset_ = 0;
  };

  basefold::Bytes state_;
};

struct RingSwitchChallengeTrace {
  std::vector<ZZ_pE> rprime_prefix;
  std::vector<ZZ_pE> rprime_suffix;
};

RingSwitchChallengeTrace ReplayRingSwitchChallenges(
    const basefold::MerkleRoot &commitment, const std::vector<ZZ_pE> &z,
    const ZZ_pE &claimed_s, const basefold::RingSwitchPCSEvalProof &proof,
    long kappa, long ell_prime) {
  TestRingSwitchTranscript transcript;
  transcript.AbsorbDigest(commitment);
  for (const ZZ_pE &zi : z) {
    transcript.AbsorbFieldElement(zi);
  }
  transcript.AbsorbFieldElement(claimed_s);
  for (const ZZ_pE &s_u : proof.s_by_u) {
    transcript.AbsorbFieldElement(s_u);
  }

  RingSwitchChallengeTrace trace;
  trace.rprime_prefix.resize(static_cast<std::size_t>(kappa));
  for (long i = 0; i < kappa; ++i) {
    trace.rprime_prefix[static_cast<std::size_t>(i)] =
        transcript.ChallengeFieldElement("rprime/prefix/" + std::to_string(i));
  }

  trace.rprime_suffix.resize(static_cast<std::size_t>(ell_prime));
  if (ell_prime == 0) {
    return trace;
  }

  transcript.AbsorbQuadraticPoly(proof.h_by_level[static_cast<std::size_t>(ell_prime - 1)]);
  for (long i = ell_prime; i-- > 0;) {
    trace.rprime_suffix[static_cast<std::size_t>(i)] =
        transcript.ChallengeFieldElement("rprime/suffix/" + std::to_string(i));
    if (i > 0) {
      transcript.AbsorbQuadraticPoly(proof.h_by_level[static_cast<std::size_t>(i - 1)]);
    }
  }
  return trace;
}

ZZ_pE SumOfPointwiseProducts(const vec_ZZ_pE &f_table,
                             const vec_ZZ_pE &g_table) {
  CHECK_EQ(f_table.length(), g_table.length());
  ZZ_pE acc;
  clear(acc);
  for (long i = 0; i < f_table.length(); ++i) {
    acc += f_table[i] * g_table[i];
  }
  return acc;
}

basefold::RingSwitchPCSParams BuildRingSwitchParamsFromSpec(
    long ell, long kappa, const ScopedRingSwitchTestContext &ctx) {
  basefold::RingSwitchPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = ctx.spec().modulus;
  input.extension_modulus = ctx.F();
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParamsFromSpec(ctx.spec(), ctx.alpha()));
  return basefold::RingSwitchPCSSetup(input);
}

basefold::RingSwitchPCSParams BuildRingSwitchParamsGR42(long ell, long kappa,
                                                        const ZZ &base_modulus,
                                                        const ZZ_pX &F,
                                                        const ZZ &p,
                                                        const ZZ_pE &alpha) {
  (void)base_modulus;
  (void)F;
  RingSwitchTestContextSpec spec = MakeGR42ContextSpec();
  spec.p = p;
  spec.modulus = base_modulus;
  spec.extension_coeffs = {1, 1, 1};
  ScopedRingSwitchTestContext ctx(spec);
  CHECK_EQ(ctx.F(), F);
  CHECK_EQ(ctx.alpha(), alpha);
  return BuildRingSwitchParamsFromSpec(ell, kappa, ctx);
}

basefold::RingSwitchPCSParams BuildRingSwitchParamsGR42D0(
    long ell, long kappa, const ZZ &base_modulus, const ZZ_pX &F, const ZZ &p,
    const ZZ_pE &alpha) {
  (void)base_modulus;
  (void)F;
  RingSwitchTestContextSpec spec = MakeGR42ContextSpec();
  spec.p = p;
  spec.modulus = base_modulus;
  spec.backend_variant = RingSwitchBackendVariant::kD0;
  ScopedRingSwitchTestContext ctx(spec);
  CHECK_EQ(ctx.F(), F);
  CHECK_EQ(ctx.alpha(), alpha);
  return BuildRingSwitchParamsFromSpec(ell, kappa, ctx);
}

basefold::RingSwitchPCSParams BuildProvidedRingSwitchParamsFromSpec(
    long ell, long kappa, const ScopedRingSwitchTestContext &ctx,
    const RingSwitchBasisCase &alpha_case, const RingSwitchBasisCase &beta_case,
    bool provide_alpha_dual = false, bool provide_beta_dual = false) {
  basefold::RingSwitchPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = ctx.spec().modulus;
  input.extension_modulus = ctx.F();
  input.use_provided_basis = true;
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.alpha_basis.basis = alpha_case.basis;
  if (provide_alpha_dual) {
    input.provided_basis.alpha_basis.dual_basis = alpha_case.dual_basis;
  }
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.beta_basis.basis = beta_case.basis;
  if (provide_beta_dual) {
    input.provided_basis.beta_basis.dual_basis = beta_case.dual_basis;
  }
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParamsFromSpec(ctx.spec(), ctx.alpha()));
  return basefold::RingSwitchPCSSetup(input);
}

basefold::RingSwitchPCSSetupInput BuildProvidedRingSwitchSetupInputFromSpec(
    long ell, long kappa, const ScopedRingSwitchTestContext &ctx) {
  basefold::RingSwitchPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = ctx.spec().modulus;
  input.extension_modulus = ctx.F();
  input.use_provided_basis = true;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParamsFromSpec(ctx.spec(), ctx.alpha()));
  return input;
}

struct ProvidedBasisProofFixture {
  basefold::RingSwitchPCSParams params;
  RingSwitchEvalCaseData eval_case;
  ZZ_pE claimed_s;
  basefold::RingSwitchPCSCommitArtifacts artifacts;
  basefold::RingSwitchPCSEvalProof proof;
};

ZZ_pE ComputeClaimedValue(const RingSwitchEvalCaseData &eval_case) {
  return basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(eval_case.t_table),
      eval_case.z);
}

ProvidedBasisProofFixture BuildProvidedBasisProofFixture(
    long ell, long kappa, const ScopedRingSwitchTestContext &ctx,
    const RingSwitchBasisCase &alpha_case, const RingSwitchBasisCase &beta_case,
    const RingSwitchEvalCaseData &eval_case, bool provide_alpha_dual = false,
    bool provide_beta_dual = false) {
  ProvidedBasisProofFixture fixture;
  fixture.params = BuildProvidedRingSwitchParamsFromSpec(
      ell, kappa, ctx, alpha_case, beta_case, provide_alpha_dual,
      provide_beta_dual);
  fixture.eval_case = eval_case;
  fixture.claimed_s = ComputeClaimedValue(fixture.eval_case);
  fixture.artifacts = basefold::RingSwitchPCSBuildCommitArtifacts(
      fixture.params, fixture.eval_case.t_table);
  fixture.proof = basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
      fixture.params, fixture.eval_case.t_table, fixture.eval_case.z,
      fixture.claimed_s, fixture.eval_case.num_queries, fixture.artifacts);
  return fixture;
}

void CheckHonestProvidedBasisProof(const ProvidedBasisProofFixture &fixture) {
  CHECK(basefold::RingSwitchPCSVerifyEval(
      fixture.params, fixture.artifacts.commitment, fixture.eval_case.z,
      fixture.claimed_s, fixture.eval_case.num_queries, fixture.proof));
}

std::vector<ZZ_pE> TakeSuffixPoint(const std::vector<ZZ_pE> &point,
                                   long suffix_dimension) {
  CHECK_GE(static_cast<long>(point.size()), suffix_dimension);
  return std::vector<ZZ_pE>(
      point.end() - static_cast<std::ptrdiff_t>(suffix_dimension), point.end());
}

void CheckPackingMatchesProvidedBetaBasis(const basefold::RingSwitchPCSParams &params,
                                          const vec_ZZ_pE &t_coeffs,
                                          const char *label) {
  const vec_ZZ_pE packed = basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
  CHECK_EQ(packed.length(), Pow2ForTest(params.ell_prime));

  for (long w = 0; w < packed.length(); ++w) {
    std::vector<NTL::ZZ_p> beta_coords(
        static_cast<std::size_t>(Pow2ForTest(params.kappa)));
    for (long v = 0; v < Pow2ForTest(params.kappa); ++v) {
      beta_coords[static_cast<std::size_t>(v)] =
          NTL::coeff(NTL::rep(t_coeffs[v + (w << params.kappa)]), 0);
    }
    CHECK_EQ(packed[w],
             ComposeFromBasisCoordsOrThrow(params.beta_basis.basis, beta_coords,
                                           label));

    const vec_ZZ_pE recovered = basefold::DecomposeGRElementToBaseCoeffs(
        params, packed[w], params.beta_basis);
    CHECK_EQ(recovered.length(), Pow2ForTest(params.kappa));
    for (long v = 0; v < recovered.length(); ++v) {
      CHECK_EQ(recovered[v], t_coeffs[v + (w << params.kappa)]);
    }
  }
}

void CheckSuffixEqualityRecoveryInProvidedAlphaBasis(
    const basefold::RingSwitchPCSParams &params,
    const std::vector<ZZ_pE> &r_suffix, const char *label) {
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long w = 0; w < num_w; ++w) {
    std::vector<NTL::ZZ_p> alpha_coords(
        static_cast<std::size_t>(tensor.basis_dimension));
    for (long u = 0; u < tensor.basis_dimension; ++u) {
      alpha_coords[static_cast<std::size_t>(u)] =
          NTL::coeff(NTL::rep(tensor.a_by_u_then_w[u * num_w + w]), 0);
    }
    const ZZ_pE reconstructed = ComposeFromBasisCoordsOrThrow(
        params.alpha_basis.basis, alpha_coords, label);
    const ZZ_pE expected =
        basefold::EqPolynomial(r_suffix,
                               BooleanPointFromIndex(w, params.ell_prime));
    CHECK_EQ(reconstructed, expected);
  }
}

void CheckPartialEvaluationRecoveryWithProvidedBases(
    const basefold::RingSwitchPCSParams &params, const vec_ZZ_pE &t_coeffs,
    const std::vector<ZZ_pE> &r_suffix, const char *label) {
  const vec_ZZ_pE packed = basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);

  vec_ZZ_pE s_by_u;
  s_by_u.SetLength(tensor.basis_dimension);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long u = 0; u < tensor.basis_dimension; ++u) {
    ZZ_pE acc;
    clear(acc);
    for (long w = 0; w < num_w; ++w) {
      acc += tensor.a_by_u_then_w[u * num_w + w] * packed[w];
    }
    s_by_u[u] = acc;
  }

  std::vector<vec_ZZ_pE> beta_coeff_rows(
      static_cast<std::size_t>(tensor.basis_dimension));
  for (long u = 0; u < tensor.basis_dimension; ++u) {
    beta_coeff_rows[static_cast<std::size_t>(u)] =
        basefold::DecomposeGRElementToBaseCoeffs(params, s_by_u[u],
                                                 params.beta_basis);
  }

  for (long v = 0; v < tensor.basis_dimension; ++v) {
    std::vector<NTL::ZZ_p> alpha_coords(
        static_cast<std::size_t>(tensor.basis_dimension));
    for (long u = 0; u < tensor.basis_dimension; ++u) {
      alpha_coords[static_cast<std::size_t>(u)] = NTL::coeff(
          NTL::rep(beta_coeff_rows[static_cast<std::size_t>(u)][v]), 0);
    }
    const ZZ_pE recovered_partial = ComposeFromBasisCoordsOrThrow(
        params.alpha_basis.basis, alpha_coords, label);

    vec_ZZ_pE slice;
    slice.SetLength(num_w);
    for (long w = 0; w < num_w; ++w) {
      slice[w] = t_coeffs[v + (w << params.kappa)];
    }
    const vec_ZZ_pE slice_monomial =
        basefold::BooleanHypercubeTableToMonomialCoeffs(slice);
    const ZZ_pE direct_partial =
        basefold::EvalMultilinearMonomialCoeffs(slice_monomial, r_suffix);
    CHECK_EQ(recovered_partial, direct_partial);
  }
}

void CheckProvidedBasisProofSizeMatchesBytes(
    const ProvidedBasisProofFixture &fixture) {
  const basefold::RingSwitchPCSOuterEvalProof outer_proof =
      basefold::RingSwitchPCSProveOuterEvalFromCommitArtifacts(
          fixture.params, fixture.eval_case.t_table, fixture.artifacts.commitment,
          fixture.eval_case.z, fixture.claimed_s, fixture.eval_case.num_queries,
          basefold::RingSwitchPCSBuildOuterCommitArtifacts(
              fixture.params, fixture.eval_case.t_table));

  CHECK_EQ(static_cast<long>(outer_proof.s_by_u.size()),
           static_cast<long>(fixture.params.beta_basis.basis.size()));

  const std::uint64_t field_elem_bytes = FixedFieldElementBytesForCurrentContext();
  const std::uint64_t expected_outer_bytes =
      1U + 8U +
      static_cast<std::uint64_t>(outer_proof.s_by_u.size()) * field_elem_bytes +
      8U +
      static_cast<std::uint64_t>(outer_proof.h_by_level.size()) * 3U *
          field_elem_bytes +
      field_elem_bytes;

  const basefold::Bytes outer_bytes =
      basefold::SerializeRingSwitchPCSOuterProofFixedBytes(fixture.params,
                                                           outer_proof);
  const std::uint64_t outer_size =
      basefold::RingSwitchPCSOuterProofSizeBytes(fixture.params, outer_proof);
  CHECK_EQ(outer_bytes.size(), static_cast<std::size_t>(outer_size));
  CHECK_EQ(outer_size, expected_outer_bytes);

  const basefold::Bytes backend_bytes = basefold::Z2kPCSBackendSerializeEvalProof(
      fixture.params.backend, fixture.proof.backend_proof);
  const std::uint64_t backend_size = basefold::Z2kPCSBackendEvalProofSizeBytes(
      fixture.params.backend, fixture.proof.backend_proof);
  CHECK_EQ(backend_bytes.size(), static_cast<std::size_t>(backend_size));

  const basefold::Bytes composed_bytes =
      basefold::SerializeRingSwitchPCSEvalProofFixedBytes(fixture.params,
                                                          fixture.proof);
  const std::uint64_t composed_size =
      basefold::RingSwitchPCSEvalProofSizeBytes(fixture.params, fixture.proof);
  CHECK_EQ(composed_bytes.size(), static_cast<std::size_t>(composed_size));
  CHECK_EQ(composed_size, outer_size + 8U + backend_size);
  CHECK(std::equal(outer_bytes.begin(), outer_bytes.end(),
                   composed_bytes.begin()));
}

basefold::Z2kPCSBackendEvalProof MutateBaseFoldBackendSubproof(
    const basefold::Z2kPCSBackendEvalProof &opaque_proof) {
  const auto basefold_ptr =
      std::static_pointer_cast<const basefold::BaseFoldPCSEvalProof>(
          opaque_proof.payload);
  CHECK(static_cast<bool>(basefold_ptr));

  basefold::BaseFoldPCSEvalProof mutated = *basefold_ptr;
  if (!mutated.h_by_level.empty()) {
    mutated.h_by_level[0].a0 += testutil::ConstZZpE(1);
  } else if (mutated.pi0_codeword.length() > 0) {
    mutated.pi0_codeword[0] += testutil::ConstZZpE(1);
  } else {
    CHECK_MSG(false, "MutateBaseFoldBackendSubproof: no mutable payload found");
  }

  basefold::Z2kPCSBackendEvalProof out = opaque_proof;
  out.payload = std::static_pointer_cast<const void>(
      std::make_shared<basefold::BaseFoldPCSEvalProof>(std::move(mutated)));
  return out;
}

void ExpectChildFailureContains(const function<void()> &fn, const string &needle,
                                const string &label) {
#if defined(__unix__) || defined(__APPLE__)
  int pipe_fds[2];
  CHECK_MSG(pipe(pipe_fds) == 0, label + ": pipe() failed");
  if (g_test_failure_count != 0) {
    return;
  }

  cout.flush();
  cerr.flush();
  const pid_t pid = fork();
  CHECK_MSG(pid >= 0, label + ": fork() failed");
  if (g_test_failure_count != 0) {
    close(pipe_fds[0]);
    close(pipe_fds[1]);
    return;
  }

  if (pid == 0) {
    close(pipe_fds[0]);
    dup2(pipe_fds[1], STDOUT_FILENO);
    dup2(pipe_fds[1], STDERR_FILENO);
    close(pipe_fds[1]);
    fn();
    _exit(0);
  }

  close(pipe_fds[1]);
  string child_output;
  char buffer[256];
  while (true) {
    const ssize_t bytes_read = read(pipe_fds[0], buffer, sizeof(buffer));
    if (bytes_read <= 0) {
      break;
    }
    child_output.append(buffer, static_cast<std::size_t>(bytes_read));
  }
  close(pipe_fds[0]);

  int status = 0;
  CHECK_MSG(waitpid(pid, &status, 0) == pid, label + ": waitpid() failed");
  if (g_test_failure_count != 0) {
    return;
  }

  const bool child_failed =
      (WIFEXITED(status) && WEXITSTATUS(status) != 0) || WIFSIGNALED(status);
  CHECK_MSG(child_failed, label + ": child unexpectedly succeeded");
  CHECK_MSG(child_output.find(needle) != string::npos,
            label + ": child output did not contain expected text: " +
                child_output);
#else
  (void)fn;
  (void)needle;
  testutil::PrintInfo(label + ": skipped negative child-process assertion on this platform");
#endif
}

void TestBaseFoldBackendAdapter_Smoke() {
  testutil::PrintInfo("Z2k backend: BaseFold adapter commit/prove/verify works");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FoldableCodeParams params = BuildParamsGR42(p, alpha);
  const basefold::Z2kPCSBackendHandle backend =
      basefold::MakeBaseFoldZ2kPCSBackend(params);

  CHECK_EQ(string(basefold::Z2kPCSBackendName(backend)), string("basefold"));
  CHECK_EQ(basefold::Z2kPCSBackendMessageLength(backend), 4L);
  CHECK_EQ(basefold::Z2kPCSBackendPointDimension(backend), 2L);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(4);
  f_coeffs[0] = testutil::ConstZZpE(0);
  f_coeffs[1] = testutil::ConstZZpE(1);
  f_coeffs[2] = alpha;
  f_coeffs[3] = testutil::ConstZZpE(2);

  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                testutil::ConstZZpE(3)};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
  const basefold::MerkleRoot commitment =
      basefold::Z2kPCSBackendCommit(backend, f_coeffs);
  const basefold::Z2kPCSBackendCommitArtifacts commit_artifacts =
      basefold::Z2kPCSBackendBuildCommitArtifacts(backend, f_coeffs);

  CHECK_EQ(commit_artifacts.commitment, commitment);
  const basefold::Z2kPCSBackendEvalProof proof = basefold::Z2kPCSBackendProveEval(
      backend, f_coeffs, z, y, /*num_queries=*/3, &commit_artifacts);

  CHECK(basefold::Z2kPCSBackendVerifyEval(backend, commitment, z, y,
                                          /*num_queries=*/3, proof));
  CHECK_GT(basefold::Z2kPCSBackendEvalProofSizeBytes(backend, proof), 0U);
}

void TestBaseFoldBackendAdapter_RejectsSameDegreeContextSwitch() {
  testutil::PrintInfo("Z2k backend: handle rejects same-degree context switches");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F1;
        SetCoeff(F1, 2, 1);
        SetCoeff(F1, 1, 1);
        SetCoeff(F1, 0, 1);
        ZZ_pEPush ext_push1(F1);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        const basefold::Z2kPCSBackendHandle backend =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

        ZZ_pX F2;
        SetCoeff(F2, 2, 1);
        SetCoeff(F2, 1, 3);
        SetCoeff(F2, 0, 1);
        ZZ_pEPush ext_push2(F2);

        basefold::RingSwitchPCSSetupInput input;
        input.ell = 3;
        input.kappa = 1;
        input.base_modulus = modulus;
        input.extension_modulus = F2;
        input.backend = backend;
        (void)basefold::RingSwitchPCSSetup(input);
      },
      "current ZZ_pE modulus does not match backend context",
      "TestBaseFoldBackendAdapter_RejectsSameDegreeContextSwitch");
}

void TestBaseFoldBackendAdapter_RejectsArtifactsFromDifferentHandle() {
  testutil::PrintInfo("Z2k backend: commit artifacts are bound to one backend instance");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F;
        SetCoeff(F, 2, 1);
        SetCoeff(F, 1, 1);
        SetCoeff(F, 0, 1);
        ZZ_pEPush ext_push(F);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        const basefold::Z2kPCSBackendHandle backend_a =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
        const basefold::Z2kPCSBackendHandle backend_b =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42Variant(p, alpha));

        vec_ZZ_pE f_coeffs;
        f_coeffs.SetLength(4);
        f_coeffs[0] = testutil::ConstZZpE(0);
        f_coeffs[1] = testutil::ConstZZpE(1);
        f_coeffs[2] = alpha;
        f_coeffs[3] = testutil::ConstZZpE(2);
        const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                      testutil::ConstZZpE(3)};
        const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

        const basefold::Z2kPCSBackendCommitArtifacts artifacts_a =
            basefold::Z2kPCSBackendBuildCommitArtifacts(backend_a, f_coeffs);
        (void)basefold::Z2kPCSBackendProveEval(backend_b, f_coeffs, z, y,
                                               /*num_queries=*/3, &artifacts_a);
      },
      "commit artifacts must belong to the same backend params instance",
      "TestBaseFoldBackendAdapter_RejectsArtifactsFromDifferentHandle");
}

void TestBaseFoldBackendAdapter_RejectsProofReuseAcrossHandles() {
  testutil::PrintInfo("Z2k backend: proofs are bound to one backend instance");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F;
        SetCoeff(F, 2, 1);
        SetCoeff(F, 1, 1);
        SetCoeff(F, 0, 1);
        ZZ_pEPush ext_push(F);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        const basefold::Z2kPCSBackendHandle backend_a =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
        const basefold::Z2kPCSBackendHandle backend_b =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42Variant(p, alpha));

        vec_ZZ_pE f_coeffs;
        f_coeffs.SetLength(4);
        f_coeffs[0] = testutil::ConstZZpE(0);
        f_coeffs[1] = testutil::ConstZZpE(1);
        f_coeffs[2] = alpha;
        f_coeffs[3] = testutil::ConstZZpE(2);
        const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                      testutil::ConstZZpE(3)};
        const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
        const basefold::MerkleRoot commitment =
            basefold::Z2kPCSBackendCommit(backend_a, f_coeffs);
        const basefold::Z2kPCSBackendEvalProof proof_a =
            basefold::Z2kPCSBackendProveEval(backend_a, f_coeffs, z, y,
                                             /*num_queries=*/3);
        (void)basefold::Z2kPCSBackendVerifyEval(backend_b, commitment, z, y,
                                                /*num_queries=*/3, proof_a);
      },
      "proof must belong to the same backend params instance",
      "TestBaseFoldBackendAdapter_RejectsProofReuseAcrossHandles");
}

void TestBaseFoldBackendAdapter_RejectsProofSizeAfterContextSwitch() {
  testutil::PrintInfo("Z2k backend: proof-size counting rejects context drift");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F1;
        SetCoeff(F1, 2, 1);
        SetCoeff(F1, 1, 1);
        SetCoeff(F1, 0, 1);
        ZZ_pEPush ext_push1(F1);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        const basefold::Z2kPCSBackendHandle backend =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
        vec_ZZ_pE f_coeffs;
        f_coeffs.SetLength(4);
        f_coeffs[0] = testutil::ConstZZpE(0);
        f_coeffs[1] = testutil::ConstZZpE(1);
        f_coeffs[2] = alpha;
        f_coeffs[3] = testutil::ConstZZpE(2);
        const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                      testutil::ConstZZpE(3)};
        const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
        const basefold::Z2kPCSBackendEvalProof proof =
            basefold::Z2kPCSBackendProveEval(backend, f_coeffs, z, y,
                                             /*num_queries=*/3);

        ZZ_pX F2;
        SetCoeff(F2, 2, 1);
        SetCoeff(F2, 1, 3);
        SetCoeff(F2, 0, 1);
        ZZ_pEPush ext_push2(F2);

        (void)basefold::Z2kPCSBackendEvalProofSizeBytes(backend, proof);
      },
      "current ZZ_pE modulus does not match backend context",
      "TestBaseFoldBackendAdapter_RejectsProofSizeAfterContextSwitch");
}

void TestRingSwitchSetup_DefaultBasesUseCurrentDegree() {
  testutil::PrintInfo("Ring-switch setup: default alpha/beta bases reflect current ZZ_pE degree");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  const basefold::RingSwitchPCSParams params = basefold::RingSwitchPCSSetup(input);
  CHECK_EQ(static_cast<long>(params.alpha_basis.basis.size()), 2L);
  CHECK_EQ(static_cast<long>(params.beta_basis.basis.size()), 2L);
  CHECK_EQ(static_cast<long>(params.alpha_basis.dual_basis.size()), 2L);
  CHECK_EQ(static_cast<long>(params.beta_basis.dual_basis.size()), 2L);
  CHECK_EQ(params.alpha_basis.basis[0], testutil::ConstZZpE(1));
  CHECK_EQ(params.beta_basis.basis[0], testutil::ConstZZpE(1));
  CHECK_EQ(params.alpha_basis.basis[1], alpha);
  CHECK_EQ(params.beta_basis.basis[1], alpha);
}

void TestValidateCurrentZ2kRingContext_SucceedsForGR42() {
  testutil::PrintInfo("Ring-switch setup: current context validation accepts GR(4,2)");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  basefold::ValidateCurrentZ2kRingContextOrThrow(modulus, F, /*kappa=*/1);
}

void TestValidateCurrentZ2kRingContext_RejectsReducibleMod2Polynomial() {
  testutil::PrintInfo("Ring-switch setup: current context validation rejects non-basic extension");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F_bad;
  SetCoeff(F_bad, 2, 1);
  SetCoeff(F_bad, 0, 1);
  ZZ_pEPush ext_push(F_bad);

  ExpectChildFailureContains(
      [&]() {
        basefold::ValidateCurrentZ2kRingContextOrThrow(modulus, F_bad,
                                                       /*kappa=*/1);
      },
      "basic irreducible modulo 2",
      "TestValidateCurrentZ2kRingContext_RejectsReducibleMod2Polynomial");
}

void TestRingSwitchSetup_Succeeds() {
  testutil::PrintInfo("Ring-switch setup: lightweight Setup validates and returns params");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  const basefold::RingSwitchPCSParams params = basefold::RingSwitchPCSSetup(input);
  CHECK_EQ(params.ell, 3L);
  CHECK_EQ(params.kappa, 1L);
  CHECK_EQ(params.ell_prime, 2L);
  CHECK_EQ(static_cast<long>(params.alpha_basis.basis.size()), 2L);
  CHECK_EQ(static_cast<long>(params.beta_basis.basis.size()), 2L);
  CHECK_EQ(static_cast<long>(params.alpha_basis.dual_basis.size()), 2L);
  CHECK_EQ(static_cast<long>(params.beta_basis.dual_basis.size()), 2L);
  CHECK_EQ(basefold::Z2kPCSBackendPointDimension(params.backend), 2L);
}

void TestRingSwitchSetup_ProvidedBasisAcceptsValidNonPolynomialBases() {
  testutil::PrintInfo("Ring-switch setup: provided alpha/beta bases can be valid non-polynomial bases");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/1, BasisTransformFamily::kDense);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = ctx.spec().modulus;
  input.extension_modulus = ctx.F();
  input.use_provided_basis = true;
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.alpha_basis.basis = alpha_case.basis;
  input.provided_basis.beta_basis.basis = beta_case.basis;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParamsFromSpec(ctx.spec(), ctx.alpha()));

  const basefold::RingSwitchPCSParams params = basefold::RingSwitchPCSSetup(input);
  CHECK_EQ(params.alpha_basis.basis, input.provided_basis.alpha_basis.basis);
  CHECK_EQ(params.beta_basis.basis, input.provided_basis.beta_basis.basis);
  CHECK_EQ(static_cast<long>(params.alpha_basis.dual_basis.size()), 2L);
  CHECK_EQ(static_cast<long>(params.beta_basis.dual_basis.size()), 2L);
}

void TestRingSwitchSetup_ProvidedBasisDerivesMissingDualBases() {
  testutil::PrintInfo("Ring-switch setup: provided alpha/beta bases derive missing dual bases");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const basefold::RingSwitchPCSParams auto_params =
      BuildRingSwitchParamsFromSpec(/*ell=*/3, /*kappa=*/1, ctx);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = ctx.spec().modulus;
  input.extension_modulus = ctx.F();
  input.use_provided_basis = true;
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.alpha_basis.basis = auto_params.alpha_basis.basis;
  input.provided_basis.beta_basis.basis = auto_params.beta_basis.basis;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParamsFromSpec(ctx.spec(), ctx.alpha()));

  const basefold::RingSwitchPCSParams params = basefold::RingSwitchPCSSetup(input);
  CHECK_EQ(params.alpha_basis.basis, auto_params.alpha_basis.basis);
  CHECK_EQ(params.beta_basis.basis, auto_params.beta_basis.basis);
  CHECK_EQ(params.alpha_basis.dual_basis, auto_params.alpha_basis.dual_basis);
  CHECK_EQ(params.beta_basis.dual_basis, auto_params.beta_basis.dual_basis);
}

void TestRingSwitchWP3Helpers_BuildReusableContextAndBasisCases() {
  testutil::PrintInfo("Ring-switch WP3: reusable context and basis helpers span variant GR(4,2) and GR(4,4)");

  {
    const ScopedRingSwitchTestContext ctx(MakeGR42VariantContextSpec());
    const RingSwitchBasisCase alpha_case = BuildTransformedBasisCase(
        /*seed=*/2, BasisTransformFamily::kUnitriangular);
    const RingSwitchBasisCase beta_case =
        BuildTransformedBasisCase(/*seed=*/3, BasisTransformFamily::kDense);
    const basefold::RingSwitchPCSParams params =
        BuildProvidedRingSwitchParamsFromSpec(
            /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case,
            /*provide_alpha_dual=*/true, /*provide_beta_dual=*/true);
    CHECK_EQ(params.alpha_basis.basis, alpha_case.basis);
    CHECK_EQ(params.beta_basis.basis, beta_case.basis);
    CHECK_EQ(params.alpha_basis.dual_basis, alpha_case.dual_basis);
    CHECK_EQ(params.beta_basis.dual_basis, beta_case.dual_basis);
  }

  {
    const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
    const basefold::RingSwitchPCSParams params =
        BuildRingSwitchParamsFromSpec(/*ell=*/4, /*kappa=*/2, ctx);
    const RingSwitchBasisCase dense_case =
        BuildTransformedBasisCase(/*seed=*/4, BasisTransformFamily::kDense);
    CHECK_EQ(static_cast<long>(params.alpha_basis.basis.size()), 4L);
    CHECK_EQ(static_cast<long>(dense_case.basis.size()), 4L);
    CHECK(IsBasisOverBaseRing(dense_case.basis));
  }
}

void TestRingSwitchSetup_RejectsMismatchedBaseModulus() {
  testutil::PrintInfo("Ring-switch setup: Setup rejects mismatched base modulus");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = to_ZZ(8);
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "current ZZ_p modulus does not match base_modulus",
      "TestRingSwitchSetup_RejectsMismatchedBaseModulus");
}

void TestRingSwitchSetup_RejectsMismatchedExtensionModulus() {
  testutil::PrintInfo("Ring-switch setup: Setup rejects mismatched extension modulus");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX F_bad;
  SetCoeff(F_bad, 2, 1);
  SetCoeff(F_bad, 0, 1);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F_bad;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "current ZZ_pE modulus does not match extension_modulus",
      "TestRingSwitchSetup_RejectsMismatchedExtensionModulus");
}

void TestRingSwitchSetup_ProvidedBasisRejectsMissingAlphaOrBeta() {
  testutil::PrintInfo("Ring-switch setup: provided-basis mode requires both alpha and beta");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.use_provided_basis = true;
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.has_beta_basis = false;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "provided alpha_basis and beta_basis are both required",
      "TestRingSwitchSetup_ProvidedBasisRejectsMissingAlphaOrBeta");
}

void TestRingSwitchSetup_ProvidedBasisRejectsMissingAlpha_GR44() {
  testutil::PrintInfo("Ring-switch WP5: GR(4,4) provided-basis mode rejects missing alpha");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/5, BasisTransformFamily::kDense);

  basefold::RingSwitchPCSSetupInput input =
      BuildProvidedRingSwitchSetupInputFromSpec(/*ell=*/4, /*kappa=*/2, ctx);
  input.provided_basis.has_alpha_basis = false;
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.beta_basis.basis = beta_case.basis;

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "provided alpha_basis and beta_basis are both required",
      "TestRingSwitchSetup_ProvidedBasisRejectsMissingAlpha_GR44");
}

void TestRingSwitchSetup_ProvidedBasisRejectsWrongDimension() {
  testutil::PrintInfo("Ring-switch setup: provided bases reject wrong dimension");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams auto_params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.use_provided_basis = true;
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.alpha_basis = auto_params.alpha_basis;
  input.provided_basis.alpha_basis.basis.pop_back();
  input.provided_basis.beta_basis = auto_params.beta_basis;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "provided_basis.alpha_basis.basis.size() must equal ZZ_pE::degree()",
      "TestRingSwitchSetup_ProvidedBasisRejectsWrongDimension");
}

void TestRingSwitchSetup_ProvidedBasisRejectsWrongBetaDimension_GR44() {
  testutil::PrintInfo("Ring-switch WP5: GR(4,4) provided bases reject wrong beta dimension");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const basefold::RingSwitchPCSParams auto_params =
      BuildRingSwitchParamsFromSpec(/*ell=*/4, /*kappa=*/2, ctx);

  basefold::RingSwitchPCSSetupInput input =
      BuildProvidedRingSwitchSetupInputFromSpec(/*ell=*/4, /*kappa=*/2, ctx);
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.alpha_basis = auto_params.alpha_basis;
  input.provided_basis.beta_basis = auto_params.beta_basis;
  input.provided_basis.beta_basis.basis.pop_back();

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "provided_basis.beta_basis.basis.size() must equal ZZ_pE::degree()",
      "TestRingSwitchSetup_ProvidedBasisRejectsWrongBetaDimension_GR44");
}

void TestRingSwitchSetup_ProvidedBasisRejectsBrokenDualBasis() {
  testutil::PrintInfo("Ring-switch setup: provided bases reject malformed dual bases");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const basefold::RingSwitchPCSParams auto_params =
      BuildRingSwitchParamsFromSpec(/*ell=*/3, /*kappa=*/1, ctx);
  const RingSwitchBasisCase broken_alpha =
      BuildBrokenDualBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = ctx.spec().modulus;
  input.extension_modulus = ctx.F();
  input.use_provided_basis = true;
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.alpha_basis.basis = broken_alpha.basis;
  input.provided_basis.alpha_basis.dual_basis = broken_alpha.dual_basis;
  input.provided_basis.beta_basis = auto_params.beta_basis;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParamsFromSpec(ctx.spec(), ctx.alpha()));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "dual-basis trace identity",
      "TestRingSwitchSetup_ProvidedBasisRejectsBrokenDualBasis");
}

void TestRingSwitchSetup_ProvidedBasisRejectsBrokenBetaDualBasis_GR44() {
  testutil::PrintInfo("Ring-switch WP5: GR(4,4) provided bases reject malformed beta dual bases");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/4, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase broken_beta =
      BuildBrokenDualBasisCase(/*seed=*/5, BasisTransformFamily::kDense);

  basefold::RingSwitchPCSSetupInput input =
      BuildProvidedRingSwitchSetupInputFromSpec(/*ell=*/4, /*kappa=*/2, ctx);
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.alpha_basis.basis = alpha_case.basis;
  input.provided_basis.alpha_basis.dual_basis = alpha_case.dual_basis;
  input.provided_basis.beta_basis.basis = broken_beta.basis;
  input.provided_basis.beta_basis.dual_basis = broken_beta.dual_basis;

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "provided_basis.beta_basis.dual_basis",
      "TestRingSwitchSetup_ProvidedBasisRejectsBrokenBetaDualBasis_GR44");
}

void TestRingSwitchSetup_ProvidedBasisRejectsNonBasisAlpha() {
  testutil::PrintInfo("Ring-switch WP5: provided setup rejects a non-basis alpha");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const RingSwitchBasisCase bad_alpha =
      BuildSingularBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/1, BasisTransformFamily::kDense);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = ctx.spec().modulus;
  input.extension_modulus = ctx.F();
  input.use_provided_basis = true;
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.alpha_basis.basis = bad_alpha.basis;
  input.provided_basis.beta_basis.basis = beta_case.basis;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParamsFromSpec(ctx.spec(), ctx.alpha()));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "provided_basis.alpha_basis.basis is not a basis over the base ring",
      "TestRingSwitchSetup_ProvidedBasisRejectsNonBasisAlpha");
}

void TestRingSwitchSetup_ProvidedBasisRejectsNonBasisBeta() {
  testutil::PrintInfo("Ring-switch WP5: provided setup rejects a non-basis beta");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase bad_beta =
      BuildSingularBasisCase(/*seed=*/1, BasisTransformFamily::kDense);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = ctx.spec().modulus;
  input.extension_modulus = ctx.F();
  input.use_provided_basis = true;
  input.provided_basis.has_alpha_basis = true;
  input.provided_basis.has_beta_basis = true;
  input.provided_basis.alpha_basis.basis = alpha_case.basis;
  input.provided_basis.beta_basis.basis = bad_beta.basis;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParamsFromSpec(ctx.spec(), ctx.alpha()));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "provided_basis.beta_basis.basis is not a basis over the base ring",
      "TestRingSwitchSetup_ProvidedBasisRejectsNonBasisBeta");
}

void TestRingSwitchSetup_RejectsBackendDimensionMismatch() {
  testutil::PrintInfo("Ring-switch setup: Setup rejects backend dimension mismatch");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 4;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "backend message length must equal 2^(ell-kappa)",
      "TestRingSwitchSetup_RejectsBackendDimensionMismatch");
}

void TestPackZ2kCoeffsToGREvals_RoundTripsSmallExample() {
  testutil::PrintInfo("Ring-switch WP3: default polynomial packing still round-trips");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_coeffs =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});

  const vec_ZZ_pE packed = basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
  CHECK_EQ(packed.length(), 4L);
  CHECK_EQ(packed[0], LongVecToZZpE({0, 1}));
  CHECK_EQ(packed[1], LongVecToZZpE({2, 3}));
  CHECK_EQ(packed[2], LongVecToZZpE({1, 0}));
  CHECK_EQ(packed[3], LongVecToZZpE({3, 2}));

  for (long w = 0; w < packed.length(); ++w) {
    const vec_ZZ_pE decomposed =
        basefold::DecomposeGRElementToBaseCoeffsPolynomialBasis(params,
                                                                packed[w]);
    CHECK_EQ(decomposed.length(), Pow2ForTest(params.kappa));
    for (long v = 0; v < decomposed.length(); ++v) {
      CHECK_EQ(decomposed[v], t_coeffs[v + (w << params.kappa)]);
    }
  }
}

void TestPackZ2kCoeffsToGREvals_ComposesAgainstProvidedBetaBasis() {
  testutil::PrintInfo("Ring-switch WP3: packing composes t'(w) against the provided beta basis");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/1, BasisTransformFamily::kDense);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 0, 3, 2},
      {{0, 1}, {1, 0}, {1, 1}});
  const vec_ZZ_pE &t_coeffs = eval_case.t_table;

  const vec_ZZ_pE packed = basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
  CHECK_EQ(packed.length(), 4L);

  for (long w = 0; w < packed.length(); ++w) {
    std::vector<NTL::ZZ_p> beta_coords(static_cast<std::size_t>(
        Pow2ForTest(params.kappa)));
    for (long v = 0; v < Pow2ForTest(params.kappa); ++v) {
      beta_coords[static_cast<std::size_t>(v)] =
          NTL::coeff(NTL::rep(t_coeffs[v + (w << params.kappa)]), 0);
    }
    CHECK_EQ(packed[w], ComposeFromBasisCoordsOrThrow(
                           params.beta_basis.basis, beta_coords,
                           "TestPackZ2kCoeffsToGREvals_ComposesAgainstProvidedBetaBasis"));

    const vec_ZZ_pE recovered = basefold::DecomposeGRElementToBaseCoeffs(
        params, packed[w], params.beta_basis);
    CHECK_EQ(recovered.length(), Pow2ForTest(params.kappa));
    for (long v = 0; v < recovered.length(); ++v) {
      CHECK_EQ(recovered[v], t_coeffs[v + (w << params.kappa)]);
    }
  }
}

void TestBooleanHypercubeTableToMonomialCoeffs_EvaluatesTableCorrectly() {
  testutil::PrintInfo("Ring-switch WP1: Boolean-table to monomial conversion matches table semantics");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  const vec_ZZ_pE table = BuildBaseRingCoeffVector({1, 2, 3, 0});
  const vec_ZZ_pE monomial =
      basefold::BooleanHypercubeTableToMonomialCoeffs(table);

  CHECK_EQ(basefold::EvalMultilinearMonomialCoeffs(
               monomial, BooleanPointFromIndex(/*index=*/0, /*dimension=*/2)),
           table[0]);
  CHECK_EQ(basefold::EvalMultilinearMonomialCoeffs(
               monomial, BooleanPointFromIndex(/*index=*/1, /*dimension=*/2)),
           table[1]);
  CHECK_EQ(basefold::EvalMultilinearMonomialCoeffs(
               monomial, BooleanPointFromIndex(/*index=*/2, /*dimension=*/2)),
           table[2]);
  CHECK_EQ(basefold::EvalMultilinearMonomialCoeffs(
               monomial, BooleanPointFromIndex(/*index=*/3, /*dimension=*/2)),
           table[3]);
}

void TestPackZ2kCoeffsToGREvals_RejectsNonBaseRingInputs() {
  testutil::PrintInfo("Ring-switch WP3: packing rejects non-constant Z2k inputs");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F;
        SetCoeff(F, 2, 1);
        SetCoeff(F, 1, 1);
        SetCoeff(F, 0, 1);
        ZZ_pEPush ext_push(F);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        const basefold::RingSwitchPCSParams params = BuildRingSwitchParamsGR42(
            /*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
        vec_ZZ_pE t_coeffs = BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});
        t_coeffs[3] = alpha;
        (void)basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
      },
      "must be a base-ring constant",
      "TestPackZ2kCoeffsToGREvals_RejectsNonBaseRingInputs");
}

void TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValues() {
  testutil::PrintInfo("Ring-switch WP3: default alpha basis still reconstructs suffix equalities");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const std::vector<ZZ_pE> r_suffix = {alpha + testutil::ConstZZpE(1),
                                       testutil::ConstZZpE(3)};
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);

  CHECK_EQ(tensor.basis_dimension, Pow2ForTest(params.kappa));
  CHECK_EQ(tensor.ell_prime, params.ell_prime);
  CHECK_EQ(tensor.a_by_u_then_w.length(), Pow2ForTest(params.ell));
  CHECK_EQ(tensor.r_table.length(), Pow2ForTest(params.ell));
  CHECK_EQ(tensor.r_monomial_coeffs.length(), Pow2ForTest(params.ell));
  CHECK_EQ(tensor.r_monomial_coeffs,
           basefold::BooleanHypercubeTableToMonomialCoeffs(tensor.r_table));

  const ZZ_pE x = PolynomialBasisElement(1);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long w = 0; w < num_w; ++w) {
    ZZ_pE reconstructed = tensor.a_by_u_then_w[w];
    reconstructed += tensor.a_by_u_then_w[num_w + w] * x;
    const ZZ_pE expected =
        basefold::EqPolynomial(r_suffix, BooleanPointFromIndex(w, params.ell_prime));
    CHECK_EQ(reconstructed, expected);
  }
}

void TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValuesWithProvidedAlphaBasis() {
  testutil::PrintInfo("Ring-switch WP3: A_{u||w} reconstructs suffix equalities in the provided alpha basis");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/1, BasisTransformFamily::kDense);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case);
  const std::vector<ZZ_pE> r_suffix = {
      ctx.alpha() + testutil::ConstZZpE(1), testutil::ConstZZpE(3)};
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);

  const long num_w = Pow2ForTest(params.ell_prime);
  for (long w = 0; w < num_w; ++w) {
    std::vector<NTL::ZZ_p> alpha_coords(static_cast<std::size_t>(
        tensor.basis_dimension));
    for (long u = 0; u < tensor.basis_dimension; ++u) {
      alpha_coords[static_cast<std::size_t>(u)] = NTL::coeff(
          NTL::rep(tensor.a_by_u_then_w[u * num_w + w]), 0);
    }
    const ZZ_pE reconstructed = ComposeFromBasisCoordsOrThrow(
        params.alpha_basis.basis, alpha_coords,
        "TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValuesWithProvidedAlphaBasis");
    const ZZ_pE expected =
        basefold::EqPolynomial(r_suffix, BooleanPointFromIndex(w, params.ell_prime));
    CHECK_EQ(reconstructed, expected);
  }
}

void TestBuildRingSwitchComponentTensor_RecoversPartialEvaluations() {
  testutil::PrintInfo("Ring-switch WP3: default alpha=beta=polynomial still recovers Appendix C.1 partial evaluations");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_coeffs =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});
  const vec_ZZ_pE packed = basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
  const std::vector<ZZ_pE> r_suffix = {alpha + testutil::ConstZZpE(1),
                                       testutil::ConstZZpE(3)};
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);

  vec_ZZ_pE s_by_u;
  s_by_u.SetLength(tensor.basis_dimension);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long u = 0; u < tensor.basis_dimension; ++u) {
    ZZ_pE acc;
    clear(acc);
    for (long w = 0; w < num_w; ++w) {
      acc += tensor.a_by_u_then_w[u * num_w + w] * packed[w];
    }
    s_by_u[u] = acc;
  }

  const ZZ_pE x = PolynomialBasisElement(1);
  for (long v = 0; v < tensor.basis_dimension; ++v) {
    ZZ_pE recovered_partial;
    clear(recovered_partial);
    for (long u = 0; u < tensor.basis_dimension; ++u) {
      const vec_ZZ_pE s_u_coeffs =
          basefold::DecomposeGRElementToBaseCoeffsPolynomialBasis(params,
                                                                  s_by_u[u]);
      const ZZ_pE alpha_u = (u == 0) ? testutil::ConstZZpE(1) : x;
      recovered_partial += s_u_coeffs[v] * alpha_u;
    }

    vec_ZZ_pE slice;
    slice.SetLength(num_w);
    for (long w = 0; w < num_w; ++w) {
      slice[w] = t_coeffs[v + (w << params.kappa)];
    }
    const vec_ZZ_pE slice_monomial =
        basefold::BooleanHypercubeTableToMonomialCoeffs(slice);
    const ZZ_pE direct_partial =
        basefold::EvalMultilinearMonomialCoeffs(slice_monomial, r_suffix);
    CHECK_EQ(recovered_partial, direct_partial);
  }
}

void TestBuildRingSwitchComponentTensor_RecoversPartialEvaluationsWithIndependentAlphaBeta() {
  testutil::PrintInfo("Ring-switch WP3: s_u recovers partial evaluations with independent provided alpha/beta bases");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/1, BasisTransformFamily::kDense);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 0, 3, 2},
      {{0, 1}, {1, 0}, {1, 1}});
  const vec_ZZ_pE &t_coeffs = eval_case.t_table;
  const vec_ZZ_pE packed = basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
  const std::vector<ZZ_pE> r_suffix = {
      ctx.alpha() + testutil::ConstZZpE(1), testutil::ConstZZpE(3)};
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);

  vec_ZZ_pE s_by_u;
  s_by_u.SetLength(tensor.basis_dimension);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long u = 0; u < tensor.basis_dimension; ++u) {
    ZZ_pE acc;
    clear(acc);
    for (long w = 0; w < num_w; ++w) {
      acc += tensor.a_by_u_then_w[u * num_w + w] * packed[w];
    }
    s_by_u[u] = acc;
  }

  std::vector<vec_ZZ_pE> beta_coeff_rows(
      static_cast<std::size_t>(tensor.basis_dimension));
  for (long u = 0; u < tensor.basis_dimension; ++u) {
    beta_coeff_rows[static_cast<std::size_t>(u)] =
        basefold::DecomposeGRElementToBaseCoeffs(params, s_by_u[u],
                                                 params.beta_basis);
  }

  for (long v = 0; v < tensor.basis_dimension; ++v) {
    std::vector<NTL::ZZ_p> alpha_coords(static_cast<std::size_t>(
        tensor.basis_dimension));
    for (long u = 0; u < tensor.basis_dimension; ++u) {
      alpha_coords[static_cast<std::size_t>(u)] = NTL::coeff(
          NTL::rep(beta_coeff_rows[static_cast<std::size_t>(u)][v]), 0);
    }
    const ZZ_pE recovered_partial = ComposeFromBasisCoordsOrThrow(
        params.alpha_basis.basis, alpha_coords,
        "TestBuildRingSwitchComponentTensor_RecoversPartialEvaluationsWithIndependentAlphaBeta");

    vec_ZZ_pE slice;
    slice.SetLength(num_w);
    for (long w = 0; w < num_w; ++w) {
      slice[w] = t_coeffs[v + (w << params.kappa)];
    }
    const vec_ZZ_pE slice_monomial =
        basefold::BooleanHypercubeTableToMonomialCoeffs(slice);
    const ZZ_pE direct_partial =
        basefold::EvalMultilinearMonomialCoeffs(slice_monomial, r_suffix);
    CHECK_EQ(recovered_partial, direct_partial);
  }
}

void TestPackZ2kCoeffsToGREvals_ComposesAgainstProvidedBetaBasis_GR44() {
  testutil::PrintInfo("Ring-switch WP4: GR(4,4) packing composes t'(w) against the provided beta basis");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/4, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/5, BasisTransformFamily::kDense);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/4, /*kappa=*/2, ctx, alpha_case, beta_case);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 3, 0, 2, 2, 0, 3, 1, 3, 2, 1, 0},
      {{0, 1, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {3, 0, 1, 0}});
  CheckPackingMatchesProvidedBetaBasis(
      params, eval_case.t_table,
      "TestPackZ2kCoeffsToGREvals_ComposesAgainstProvidedBetaBasis_GR44");
}

void TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValuesWithProvidedAlphaBasis_GR44() {
  testutil::PrintInfo("Ring-switch WP4: GR(4,4) A_{u||w} reconstructs suffix equalities in the provided alpha basis");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/4, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/5, BasisTransformFamily::kDense);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/4, /*kappa=*/2, ctx, alpha_case, beta_case);
  const std::vector<ZZ_pE> r_suffix = {
      BuildElementFromPolynomialCoords({1, 1, 0, 0}),
      BuildElementFromPolynomialCoords({3, 0, 1, 0})};
  CheckSuffixEqualityRecoveryInProvidedAlphaBasis(
      params, r_suffix,
      "TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValuesWithProvidedAlphaBasis_GR44");
}

void TestBuildRingSwitchComponentTensor_RecoversPartialEvaluationsWithIndependentAlphaBeta_GR44() {
  testutil::PrintInfo("Ring-switch WP4: GR(4,4) s_u recovers partial evaluations with independent provided alpha/beta bases");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/4, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/5, BasisTransformFamily::kDense);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/4, /*kappa=*/2, ctx, alpha_case, beta_case);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 3, 0, 2, 2, 0, 3, 1, 3, 2, 1, 0},
      {{0, 1, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {3, 0, 1, 0}});
  CheckPartialEvaluationRecoveryWithProvidedBases(
      params, eval_case.t_table, TakeSuffixPoint(eval_case.z, params.ell_prime),
      "TestBuildRingSwitchComponentTensor_RecoversPartialEvaluationsWithIndependentAlphaBeta_GR44");
}

void TestBuildRingSwitchComponentTensor_RCoeffsEvaluateAsExpected() {
  testutil::PrintInfo("Ring-switch WP3: verifier polynomial r evaluates with v||w flattening");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const std::vector<ZZ_pE> r_suffix = {alpha + testutil::ConstZZpE(1),
                                       testutil::ConstZZpE(3)};
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);

  const std::vector<ZZ_pE> r_prime = {alpha, testutil::ConstZZpE(1),
                                      alpha + testutil::ConstZZpE(1)};
  const ZZ_pE got =
      basefold::EvalMultilinearMonomialCoeffs(tensor.r_monomial_coeffs, r_prime);

  const std::vector<ZZ_pE> r_prefix = {r_prime[0]};
  const std::vector<ZZ_pE> r_suffix_prime = {r_prime[1], r_prime[2]};
  ZZ_pE expected;
  clear(expected);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long w = 0; w < num_w; ++w) {
    const ZZ_pE eq_w =
        basefold::EqPolynomial(r_suffix_prime, BooleanPointFromIndex(w, params.ell_prime));
    ZZ_pE inner;
    clear(inner);
    for (long u = 0; u < tensor.basis_dimension; ++u) {
      const ZZ_pE eq_u =
          basefold::EqPolynomial(r_prefix, BooleanPointFromIndex(u, params.kappa));
      inner += tensor.a_by_u_then_w[u * num_w + w] * eq_u;
    }
    expected += inner * eq_w;
  }
  CHECK_EQ(got, expected);
}

void TestBuildRingSwitchComponentTensor_RejectsWrongSuffixDimension() {
  testutil::PrintInfo("Ring-switch WP3: component tensor rejects wrong suffix dimension");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F;
        SetCoeff(F, 2, 1);
        SetCoeff(F, 1, 1);
        SetCoeff(F, 0, 1);
        ZZ_pEPush ext_push(F);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        const basefold::RingSwitchPCSParams params = BuildRingSwitchParamsGR42(
            /*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
        const std::vector<ZZ_pE> wrong_suffix = {alpha};
        (void)basefold::BuildRingSwitchComponentTensor(params, wrong_suffix);
      },
      "r_suffix dimension must equal ell_prime",
      "TestBuildRingSwitchComponentTensor_RejectsWrongSuffixDimension");
}

void TestRingSwitchCommit_MatchesDirectBackendCommit() {
  testutil::PrintInfo("Ring-switch WP3: compiler commitment matches direct backend commitment on packed polynomial");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});

  const basefold::MerkleRoot compiler_commitment =
      basefold::RingSwitchPCSCommit(params, t_table);

  const vec_ZZ_pE packed_table = basefold::PackZ2kCoeffsToGREvals(params, t_table);
  const vec_ZZ_pE packed_monomial =
      basefold::BooleanHypercubeTableToMonomialCoeffs(packed_table);
  const basefold::MerkleRoot direct_backend_commitment =
      basefold::Z2kPCSBackendCommit(params.backend, packed_monomial);

  CHECK_EQ(compiler_commitment, direct_backend_commitment);
}

void TestRingSwitchBuildCommitArtifacts_CachesPackedRepresentations() {
  testutil::PrintInfo("Ring-switch WP3: commit artifacts cache both packed representations and backend artifacts");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 2, 1, 0, 1, 2, 3, 0});

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);

  const vec_ZZ_pE expected_packed_table =
      basefold::PackZ2kCoeffsToGREvals(params, t_table);
  const vec_ZZ_pE expected_packed_monomial =
      basefold::BooleanHypercubeTableToMonomialCoeffs(expected_packed_table);
  const basefold::MerkleRoot direct_backend_commitment =
      basefold::Z2kPCSBackendCommit(params.backend, expected_packed_monomial);

  CHECK_EQ(artifacts.t_packed_table, expected_packed_table);
  CHECK_EQ(artifacts.t_packed_monomial_coeffs, expected_packed_monomial);
  CHECK_EQ(artifacts.commitment, direct_backend_commitment);
  CHECK_EQ(artifacts.backend_commit_artifacts.commitment,
           direct_backend_commitment);
}

void TestRingSwitchProveEvalFromCommitArtifacts_HonestProofIsSelfConsistent() {
  testutil::PrintInfo("Ring-switch WP4: artifact-backed prove path yields a self-consistent honest proof");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});
  const std::vector<ZZ_pE> z = {alpha, testutil::ConstZZpE(1),
                                alpha + testutil::ConstZZpE(1)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK_EQ(static_cast<long>(proof.s_by_u.size()), Pow2ForTest(params.kappa));
  CHECK_EQ(static_cast<long>(proof.h_by_level.size()), params.ell_prime);

  const std::vector<ZZ_pE> z_suffix(z.begin() + params.kappa, z.end());
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, z_suffix);

  const ZZ_pE basis_x = PolynomialBasisElement(1);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long v = 0; v < tensor.basis_dimension; ++v) {
    ZZ_pE recovered_partial;
    clear(recovered_partial);
    for (long u = 0; u < tensor.basis_dimension; ++u) {
      const vec_ZZ_pE s_u_coeffs =
          basefold::DecomposeGRElementToBaseCoeffsPolynomialBasis(
              params, proof.s_by_u[static_cast<std::size_t>(u)]);
      const ZZ_pE alpha_u = (u == 0) ? testutil::ConstZZpE(1) : basis_x;
      recovered_partial += s_u_coeffs[v] * alpha_u;
    }

    vec_ZZ_pE slice;
    slice.SetLength(num_w);
    for (long w = 0; w < num_w; ++w) {
      slice[w] = t_table[v + (w << params.kappa)];
    }
    const ZZ_pE direct_partial = basefold::EvalMultilinearMonomialCoeffs(
        basefold::BooleanHypercubeTableToMonomialCoeffs(slice), z_suffix);
    CHECK_EQ(recovered_partial, direct_partial);
  }

  const RingSwitchChallengeTrace trace = ReplayRingSwitchChallenges(
      artifacts.commitment, z, claimed_s, proof, params.kappa, params.ell_prime);
  ZZ_pE initial_claim;
  clear(initial_claim);
  for (long u = 0; u < tensor.basis_dimension; ++u) {
    initial_claim +=
        proof.s_by_u[static_cast<std::size_t>(u)] *
        basefold::EqPolynomial(
            trace.rprime_prefix,
            BooleanPointFromIndex(u, params.kappa));
  }
  CHECK(basefold::CheckProductSumcheckChain(initial_claim, proof.h_by_level,
                                            trace.rprime_suffix));

  CHECK_EQ(proof.t_star, basefold::EvalMultilinearMonomialCoeffs(
                             artifacts.t_packed_monomial_coeffs,
                             trace.rprime_suffix));

  std::vector<ZZ_pE> rprime_full = trace.rprime_prefix;
  rprime_full.insert(rprime_full.end(), trace.rprime_suffix.begin(),
                     trace.rprime_suffix.end());
  const ZZ_pE g_star =
      basefold::EvalMultilinearMonomialCoeffs(tensor.r_monomial_coeffs,
                                              rprime_full);
  const ZZ_pE final_claim = proof.h_by_level.empty()
                                ? initial_claim
                                : proof.h_by_level[0].Eval(trace.rprime_suffix[0]);
  CHECK_EQ(final_claim, proof.t_star * g_star);

  CHECK(basefold::Z2kPCSBackendVerifyEval(
      params.backend, artifacts.commitment, trace.rprime_suffix, proof.t_star,
      /*num_queries=*/2, proof.backend_proof));
}

void TestRingSwitchProveEval_DirectAndArtifactPathsAgreeOnOuterMessages() {
  testutil::PrintInfo("Ring-switch WP4: direct and artifact-backed prove paths agree on outer proof messages");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), alpha,
                                testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSEvalProof direct =
      basefold::RingSwitchPCSProveEval(params, t_table, z, claimed_s,
                                      /*num_queries=*/2);
  const basefold::RingSwitchPCSEvalProof cached =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK_EQ(direct.s_by_u, cached.s_by_u);
  CHECK_EQ(static_cast<long>(direct.h_by_level.size()),
           static_cast<long>(cached.h_by_level.size()));
  for (long i = 0; i < static_cast<long>(direct.h_by_level.size()); ++i) {
    CHECK_EQ(direct.h_by_level[static_cast<std::size_t>(i)].a0,
             cached.h_by_level[static_cast<std::size_t>(i)].a0);
    CHECK_EQ(direct.h_by_level[static_cast<std::size_t>(i)].a1,
             cached.h_by_level[static_cast<std::size_t>(i)].a1);
    CHECK_EQ(direct.h_by_level[static_cast<std::size_t>(i)].a2,
             cached.h_by_level[static_cast<std::size_t>(i)].a2);
  }
  CHECK_EQ(direct.t_star, cached.t_star);
  CHECK_EQ(basefold::Z2kPCSBackendEvalProofSizeBytes(params.backend,
                                                     direct.backend_proof),
           basefold::Z2kPCSBackendEvalProofSizeBytes(params.backend,
                                                     cached.backend_proof));
}

void TestRingSwitchVerifyEval_AcceptsHonestProof() {
  testutil::PrintInfo("Ring-switch WP5: verifier accepts an honest proof");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});
  const std::vector<ZZ_pE> z = {alpha, testutil::ConstZZpE(1),
                                alpha + testutil::ConstZZpE(1)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK(basefold::RingSwitchPCSVerifyEval(params, artifacts.commitment, z,
                                          claimed_s, /*num_queries=*/2, proof));
}

void TestRingSwitchVerifyEval_AcceptsHonestProofFromDirectProvePath() {
  testutil::PrintInfo("Ring-switch WP6: direct prove path verifies end-to-end");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), alpha,
                                testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::MerkleRoot commitment =
      basefold::RingSwitchPCSCommit(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEval(params, t_table, z, claimed_s,
                                       /*num_queries=*/2);

  CHECK(basefold::RingSwitchPCSVerifyEval(params, commitment, z, claimed_s,
                                          /*num_queries=*/2, proof));
}

void TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedAlphaBetaBases() {
  testutil::PrintInfo("Ring-switch WP3: verifier accepts an honest proof with provided non-polynomial alpha/beta bases");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/1, BasisTransformFamily::kDense);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 0, 3, 2},
      {{0, 1}, {1, 0}, {1, 1}});
  const vec_ZZ_pE &t_table = eval_case.t_table;
  const std::vector<ZZ_pE> &z = eval_case.z;
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK(basefold::RingSwitchPCSVerifyEval(params, artifacts.commitment, z,
                                          claimed_s, /*num_queries=*/2, proof));
}

void TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedNonPolynomialAlphaOnly() {
  testutil::PrintInfo("Ring-switch WP5: verifier accepts an honest proof when only alpha is non-polynomial");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const basefold::RingSwitchPCSParams default_params =
      BuildRingSwitchParamsFromSpec(/*ell=*/3, /*kappa=*/1, ctx);
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      MakeBasisCase("default-beta", default_params.beta_basis);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {3, 1, 0, 2, 1, 2, 3, 0},
      {{1, 1}, {0, 1}, {3, 0}});
  const vec_ZZ_pE &t_table = eval_case.t_table;
  const std::vector<ZZ_pE> &z = eval_case.z;
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::MerkleRoot commitment =
      basefold::RingSwitchPCSCommit(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEval(params, t_table, z, claimed_s,
                                       /*num_queries=*/2);

  CHECK(basefold::RingSwitchPCSVerifyEval(params, commitment, z, claimed_s,
                                          /*num_queries=*/2, proof));
}

void TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedNonPolynomialBetaOnly() {
  testutil::PrintInfo("Ring-switch WP5: verifier accepts an honest proof when only beta is non-polynomial");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const basefold::RingSwitchPCSParams default_params =
      BuildRingSwitchParamsFromSpec(/*ell=*/3, /*kappa=*/1, ctx);
  const RingSwitchBasisCase alpha_case =
      MakeBasisCase("default-alpha", default_params.alpha_basis);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/1, BasisTransformFamily::kDense);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {3, 1, 0, 2, 1, 2, 3, 0},
      {{1, 1}, {0, 1}, {3, 0}});
  const vec_ZZ_pE &t_table = eval_case.t_table;
  const std::vector<ZZ_pE> &z = eval_case.z;
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::MerkleRoot commitment =
      basefold::RingSwitchPCSCommit(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEval(params, t_table, z, claimed_s,
                                       /*num_queries=*/2);

  CHECK(basefold::RingSwitchPCSVerifyEval(params, commitment, z, claimed_s,
                                          /*num_queries=*/2, proof));
}

void TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedExplicitDualBases() {
  testutil::PrintInfo("Ring-switch WP4: verifier accepts an honest proof with provided explicit alpha/beta dual bases");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/1, BasisTransformFamily::kDense);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 0, 3, 2},
      {{0, 1}, {1, 0}, {1, 1}});
  const ProvidedBasisProofFixture fixture = BuildProvidedBasisProofFixture(
      /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case, eval_case,
      /*provide_alpha_dual=*/true, /*provide_beta_dual=*/true);
  CheckHonestProvidedBasisProof(fixture);
}

void TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedAlphaBetaBasesOnVariantBackend() {
  testutil::PrintInfo("Ring-switch WP4: verifier accepts an honest provided alpha!=beta proof on the GR(4,2) variant backend");

  const ScopedRingSwitchTestContext ctx(MakeGR42VariantContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/2, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/3, BasisTransformFamily::kDense);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 3, 1, 2, 2, 1, 3, 0},
      {{1, 0}, {1, 1}, {3, 1}});
  const ProvidedBasisProofFixture fixture = BuildProvidedBasisProofFixture(
      /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case, eval_case);
  CheckHonestProvidedBasisProof(fixture);
}

void TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedNonPolynomialAlphaOnly_GR44() {
  testutil::PrintInfo("Ring-switch WP4: verifier accepts an honest GR(4,4) proof when only alpha is non-polynomial");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const basefold::RingSwitchPCSParams default_params =
      BuildRingSwitchParamsFromSpec(/*ell=*/4, /*kappa=*/2, ctx);
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/4, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      MakeBasisCase("default-beta", default_params.beta_basis);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 3, 0, 2, 2, 0, 3, 1, 3, 2, 1, 0},
      {{0, 1, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {3, 0, 1, 0}});
  const ProvidedBasisProofFixture fixture = BuildProvidedBasisProofFixture(
      /*ell=*/4, /*kappa=*/2, ctx, alpha_case, beta_case, eval_case);
  CheckHonestProvidedBasisProof(fixture);
}

void TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedNonPolynomialBetaOnly_GR44() {
  testutil::PrintInfo("Ring-switch WP4: verifier accepts an honest GR(4,4) proof when only beta is non-polynomial");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const basefold::RingSwitchPCSParams default_params =
      BuildRingSwitchParamsFromSpec(/*ell=*/4, /*kappa=*/2, ctx);
  const RingSwitchBasisCase alpha_case =
      MakeBasisCase("default-alpha", default_params.alpha_basis);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/5, BasisTransformFamily::kDense);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 3, 0, 2, 2, 0, 3, 1, 3, 2, 1, 0},
      {{0, 1, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {3, 0, 1, 0}});
  const ProvidedBasisProofFixture fixture = BuildProvidedBasisProofFixture(
      /*ell=*/4, /*kappa=*/2, ctx, alpha_case, beta_case, eval_case);
  CheckHonestProvidedBasisProof(fixture);
}

void TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedAlphaBetaBases_GR44() {
  testutil::PrintInfo("Ring-switch WP4: verifier accepts an honest GR(4,4) proof with provided non-polynomial alpha!=beta bases");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/4, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/5, BasisTransformFamily::kDense);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 3, 0, 2, 2, 0, 3, 1, 3, 2, 1, 0},
      {{0, 1, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {3, 0, 1, 0}});
  const ProvidedBasisProofFixture fixture = BuildProvidedBasisProofFixture(
      /*ell=*/4, /*kappa=*/2, ctx, alpha_case, beta_case, eval_case);
  CheckHonestProvidedBasisProof(fixture);
}

void TestRingSwitchPaperAPI_AcceptsHonestProof() {
  testutil::PrintInfo("Ring-switch WP6: staged setup/commit/prove/verify API matches the legacy flow");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
  const basefold::RingSwitchPCSSetupOutput api =
      basefold::RingSwitchPCSSetupProtocol(input);

  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), alpha,
                                testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommittedWitness committed =
      api.prover.Commit(t_table);
  CHECK_EQ(committed.commitment,
           basefold::RingSwitchPCSCommit(api.params, t_table));
  CHECK_EQ(committed.commitment, committed.commit_artifacts.commitment);

  const basefold::RingSwitchPCSEvalProof proof =
      api.prover.Prove(committed, z, claimed_s, /*num_queries=*/2);
  CHECK(api.verifier.Verify(committed.commitment, z, claimed_s,
                            /*num_queries=*/2, proof));
}

void TestRingSwitchVerifyEval_RejectsTampering() {
  testutil::PrintInfo("Ring-switch WP5: verifier rejects wrong claim, outer tampering, and backend tampering");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), alpha,
                                testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK(basefold::RingSwitchPCSVerifyEval(params, artifacts.commitment, z,
                                          claimed_s, /*num_queries=*/2, proof));

  CHECK(!basefold::RingSwitchPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s + testutil::ConstZZpE(1),
      /*num_queries=*/2, proof));

  basefold::RingSwitchPCSEvalProof tampered_s_by_u = proof;
  tampered_s_by_u.s_by_u[0] += testutil::ConstZZpE(1);
  CHECK(!basefold::RingSwitchPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      tampered_s_by_u));

  basefold::RingSwitchPCSEvalProof tampered_h = proof;
  tampered_h.h_by_level[0].a0 += testutil::ConstZZpE(1);
  CHECK(!basefold::RingSwitchPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      tampered_h));

  basefold::RingSwitchPCSEvalProof tampered_t_star = proof;
  tampered_t_star.t_star += testutil::ConstZZpE(1);
  CHECK(!basefold::RingSwitchPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      tampered_t_star));

  basefold::RingSwitchPCSEvalProof tampered_backend = proof;
  tampered_backend.backend_proof =
      MutateBaseFoldBackendSubproof(proof.backend_proof);
  CHECK(!basefold::RingSwitchPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      tampered_backend));
}

void TestRingSwitchVerifyEval_DimensionZeroUsesNoSumcheckRounds() {
  testutil::PrintInfo("Ring-switch WP5: verifier handles ell_prime=0 as a zero-round sumcheck");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params = BuildRingSwitchParamsGR42D0(
      /*ell=*/1, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table = BuildBaseRingCoeffVector({1, 3});
  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK(proof.h_by_level.empty());
  CHECK(basefold::RingSwitchPCSVerifyEval(params, artifacts.commitment, z,
                                          claimed_s, /*num_queries=*/2, proof));
}

void TestRingSwitchOuterProveVerify_AcceptsHonestProof() {
  testutil::PrintInfo("Ring-switch outer: prover/verifier accepts an honest outer-only proof");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), alpha,
                                testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts composed_artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSOuterCommitArtifacts outer_artifacts =
      basefold::RingSwitchPCSBuildOuterCommitArtifacts(params, t_table);
  CHECK_EQ(outer_artifacts.t_packed_table, composed_artifacts.t_packed_table);
  CHECK_EQ(outer_artifacts.t_packed_monomial_coeffs,
           composed_artifacts.t_packed_monomial_coeffs);

  const basefold::RingSwitchPCSOuterEvalProof outer_proof =
      basefold::RingSwitchPCSProveOuterEvalFromCommitArtifacts(
          params, t_table, composed_artifacts.commitment, z, claimed_s,
          /*num_queries=*/2, outer_artifacts);
  CHECK(basefold::RingSwitchPCSVerifyOuterEval(
      params, composed_artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      outer_proof));

  const basefold::RingSwitchPCSEvalProof composed_proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, composed_artifacts);
  CHECK_EQ(outer_proof.s_by_u, composed_proof.s_by_u);
  CHECK_EQ(outer_proof.h_by_level.size(), composed_proof.h_by_level.size());
  for (long i = 0; i < static_cast<long>(outer_proof.h_by_level.size()); ++i) {
    CHECK_EQ(outer_proof.h_by_level[static_cast<std::size_t>(i)].a0,
             composed_proof.h_by_level[static_cast<std::size_t>(i)].a0);
    CHECK_EQ(outer_proof.h_by_level[static_cast<std::size_t>(i)].a1,
             composed_proof.h_by_level[static_cast<std::size_t>(i)].a1);
    CHECK_EQ(outer_proof.h_by_level[static_cast<std::size_t>(i)].a2,
             composed_proof.h_by_level[static_cast<std::size_t>(i)].a2);
  }
  CHECK_EQ(outer_proof.t_star, composed_proof.t_star);
}

void TestRingSwitchOuterProveFromCommitArtifacts_CheckedAndUncheckedAgree() {
  testutil::PrintInfo(
      "Ring-switch outer: checked and unchecked outer prove paths agree on honest inputs");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), alpha,
                                testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts composed_artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSOuterCommitArtifacts outer_artifacts =
      basefold::RingSwitchPCSBuildOuterCommitArtifacts(params, t_table);

  const basefold::RingSwitchPCSOuterEvalProof checked_proof =
      basefold::RingSwitchPCSProveOuterEvalFromCommitArtifacts(
          params, t_table, composed_artifacts.commitment, z, claimed_s,
          /*num_queries=*/2, outer_artifacts);
  const basefold::RingSwitchPCSOuterEvalProof unchecked_proof =
      basefold::RingSwitchPCSProveOuterEvalFromCommitArtifactsUnchecked(
          params, t_table, composed_artifacts.commitment, z, claimed_s,
          /*num_queries=*/2, outer_artifacts);

  CHECK(basefold::RingSwitchPCSVerifyOuterEval(
      params, composed_artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      checked_proof));
  CHECK(basefold::RingSwitchPCSVerifyOuterEval(
      params, composed_artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      unchecked_proof));
  CHECK_EQ(checked_proof.s_by_u, unchecked_proof.s_by_u);
  CHECK_EQ(checked_proof.h_by_level.size(), unchecked_proof.h_by_level.size());
  for (long i = 0; i < static_cast<long>(checked_proof.h_by_level.size()); ++i) {
    CHECK_EQ(checked_proof.h_by_level[static_cast<std::size_t>(i)].a0,
             unchecked_proof.h_by_level[static_cast<std::size_t>(i)].a0);
    CHECK_EQ(checked_proof.h_by_level[static_cast<std::size_t>(i)].a1,
             unchecked_proof.h_by_level[static_cast<std::size_t>(i)].a1);
    CHECK_EQ(checked_proof.h_by_level[static_cast<std::size_t>(i)].a2,
             unchecked_proof.h_by_level[static_cast<std::size_t>(i)].a2);
  }
  CHECK_EQ(checked_proof.t_star, unchecked_proof.t_star);
}

void TestRingSwitchProveEvalFromCommitArtifacts_CheckedAndUncheckedAgree() {
  testutil::PrintInfo(
      "Ring-switch outer: checked and unchecked composed prove paths agree on honest inputs");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), alpha,
                                testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSEvalProof checked_proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);
  const basefold::RingSwitchPCSEvalProof unchecked_proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifactsUnchecked(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK(basefold::RingSwitchPCSVerifyEval(params, artifacts.commitment, z,
                                          claimed_s, /*num_queries=*/2,
                                          checked_proof));
  CHECK(basefold::RingSwitchPCSVerifyEval(params, artifacts.commitment, z,
                                          claimed_s, /*num_queries=*/2,
                                          unchecked_proof));
  CHECK_EQ(checked_proof.s_by_u, unchecked_proof.s_by_u);
  CHECK_EQ(checked_proof.h_by_level.size(),
           unchecked_proof.h_by_level.size());
  for (long i = 0; i < static_cast<long>(checked_proof.h_by_level.size()); ++i) {
    CHECK_EQ(checked_proof.h_by_level[static_cast<std::size_t>(i)].a0,
             unchecked_proof.h_by_level[static_cast<std::size_t>(i)].a0);
    CHECK_EQ(checked_proof.h_by_level[static_cast<std::size_t>(i)].a1,
             unchecked_proof.h_by_level[static_cast<std::size_t>(i)].a1);
    CHECK_EQ(checked_proof.h_by_level[static_cast<std::size_t>(i)].a2,
             unchecked_proof.h_by_level[static_cast<std::size_t>(i)].a2);
  }
  CHECK_EQ(checked_proof.t_star, unchecked_proof.t_star);

  const basefold::Bytes checked_bytes =
      basefold::SerializeRingSwitchPCSEvalProofFixedBytes(params,
                                                          checked_proof);
  const basefold::Bytes unchecked_bytes =
      basefold::SerializeRingSwitchPCSEvalProofFixedBytes(params,
                                                          unchecked_proof);
  CHECK_EQ(checked_bytes.size(), unchecked_bytes.size());
  CHECK(std::equal(checked_bytes.begin(), checked_bytes.end(),
                   unchecked_bytes.begin()));
}

void TestRingSwitchProofSerialize_ComposedSizeMatchesBytes() {
  testutil::PrintInfo("Ring-switch WP7: serializer bytes match outer and composed proof-size accounting");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), alpha,
                                testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);
  const basefold::RingSwitchPCSOuterEvalProof outer_proof =
      basefold::RingSwitchPCSProveOuterEvalFromCommitArtifacts(
          params, t_table, artifacts.commitment, z, claimed_s,
          /*num_queries=*/2,
          basefold::RingSwitchPCSBuildOuterCommitArtifacts(params, t_table));

  const std::uint64_t field_elem_bytes = FixedFieldElementBytesForCurrentContext();
  const std::uint64_t expected_outer_bytes =
      1U + 8U +
      static_cast<std::uint64_t>(outer_proof.s_by_u.size()) * field_elem_bytes + 8U +
      static_cast<std::uint64_t>(outer_proof.h_by_level.size()) * 3U *
          field_elem_bytes +
      field_elem_bytes;

  const basefold::Bytes outer_bytes =
      basefold::SerializeRingSwitchPCSOuterProofFixedBytes(params, outer_proof);
  const std::uint64_t outer_size =
      basefold::RingSwitchPCSOuterProofSizeBytes(params, outer_proof);
  CHECK_EQ(outer_bytes.size(), static_cast<std::size_t>(outer_size));
  CHECK_EQ(outer_size, expected_outer_bytes);

  const basefold::Bytes backend_bytes =
      basefold::Z2kPCSBackendSerializeEvalProof(params.backend,
                                                proof.backend_proof);
  const std::uint64_t backend_size =
      basefold::Z2kPCSBackendEvalProofSizeBytes(params.backend,
                                                proof.backend_proof);
  CHECK_EQ(backend_bytes.size(), static_cast<std::size_t>(backend_size));

  const basefold::Bytes composed_bytes =
      basefold::SerializeRingSwitchPCSEvalProofFixedBytes(params, proof);
  const std::uint64_t composed_size =
      basefold::RingSwitchPCSEvalProofSizeBytes(params, proof);
  CHECK_EQ(composed_bytes.size(), static_cast<std::size_t>(composed_size));
  CHECK_EQ(composed_size, outer_size + 8U + backend_size);
  CHECK(std::equal(outer_bytes.begin(), outer_bytes.end(),
                   composed_bytes.begin()));
}

void TestRingSwitchProofSerialize_ComposedSizeMatchesBytesWithProvidedAlphaBeta() {
  testutil::PrintInfo("Ring-switch WP4: serializer and proof-size helpers stay exact under provided non-polynomial alpha/beta bases");

  const ScopedRingSwitchTestContext ctx(MakeGR42ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/0, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/1, BasisTransformFamily::kDense);
  const basefold::RingSwitchPCSParams params = BuildProvidedRingSwitchParamsFromSpec(
      /*ell=*/3, /*kappa=*/1, ctx, alpha_case, beta_case);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {3, 1, 0, 2, 1, 2, 3, 0},
      {{1, 1}, {0, 1}, {3, 0}});
  const vec_ZZ_pE &t_table = eval_case.t_table;
  const std::vector<ZZ_pE> &z = eval_case.z;
  const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
      basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

  const basefold::RingSwitchPCSCommitArtifacts artifacts =
      basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);
  const basefold::RingSwitchPCSEvalProof proof =
      basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);
  const basefold::RingSwitchPCSOuterEvalProof outer_proof =
      basefold::RingSwitchPCSProveOuterEvalFromCommitArtifacts(
          params, t_table, artifacts.commitment, z, claimed_s,
          /*num_queries=*/2,
          basefold::RingSwitchPCSBuildOuterCommitArtifacts(params, t_table));

  CHECK_EQ(static_cast<long>(outer_proof.s_by_u.size()),
           static_cast<long>(params.beta_basis.basis.size()));

  const std::uint64_t field_elem_bytes = FixedFieldElementBytesForCurrentContext();
  const std::uint64_t expected_outer_bytes =
      1U + 8U +
      static_cast<std::uint64_t>(outer_proof.s_by_u.size()) * field_elem_bytes + 8U +
      static_cast<std::uint64_t>(outer_proof.h_by_level.size()) * 3U *
          field_elem_bytes +
      field_elem_bytes;

  const basefold::Bytes outer_bytes =
      basefold::SerializeRingSwitchPCSOuterProofFixedBytes(params, outer_proof);
  const std::uint64_t outer_size =
      basefold::RingSwitchPCSOuterProofSizeBytes(params, outer_proof);
  CHECK_EQ(outer_bytes.size(), static_cast<std::size_t>(outer_size));
  CHECK_EQ(outer_size, expected_outer_bytes);

  const basefold::Bytes backend_bytes =
      basefold::Z2kPCSBackendSerializeEvalProof(params.backend,
                                                proof.backend_proof);
  const std::uint64_t backend_size =
      basefold::Z2kPCSBackendEvalProofSizeBytes(params.backend,
                                                proof.backend_proof);
  CHECK_EQ(backend_bytes.size(), static_cast<std::size_t>(backend_size));

  const basefold::Bytes composed_bytes =
      basefold::SerializeRingSwitchPCSEvalProofFixedBytes(params, proof);
  const std::uint64_t composed_size =
      basefold::RingSwitchPCSEvalProofSizeBytes(params, proof);
  CHECK_EQ(composed_bytes.size(), static_cast<std::size_t>(composed_size));
  CHECK_EQ(composed_size, outer_size + 8U + backend_size);
  CHECK(std::equal(outer_bytes.begin(), outer_bytes.end(),
                   composed_bytes.begin()));
}

void TestRingSwitchProofSerialize_ComposedSizeMatchesBytesWithProvidedAlphaBeta_GR44() {
  testutil::PrintInfo("Ring-switch WP4: GR(4,4) serializer and proof-size helpers stay exact under provided non-polynomial alpha/beta bases");

  const ScopedRingSwitchTestContext ctx(MakeGR44ContextSpec());
  const RingSwitchBasisCase alpha_case =
      BuildTransformedBasisCase(/*seed=*/4, BasisTransformFamily::kUnitriangular);
  const RingSwitchBasisCase beta_case =
      BuildTransformedBasisCase(/*seed=*/5, BasisTransformFamily::kDense);
  const RingSwitchEvalCaseData eval_case = BuildEvalCaseData(
      {0, 1, 2, 3, 1, 3, 0, 2, 2, 0, 3, 1, 3, 2, 1, 0},
      {{0, 1, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {3, 0, 1, 0}});
  const ProvidedBasisProofFixture fixture = BuildProvidedBasisProofFixture(
      /*ell=*/4, /*kappa=*/2, ctx, alpha_case, beta_case, eval_case);
  CheckProvidedBasisProofSizeMatchesBytes(fixture);
}

void TestProductSumcheckProver_BooleanTablesPasses() {
  testutil::PrintInfo("Ring-switch WP2: product sumcheck passes on honest Boolean tables");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  vec_ZZ_pE f_table;
  f_table.SetLength(4);
  f_table[0] = testutil::ConstZZpE(1);
  f_table[1] = alpha;
  f_table[2] = testutil::ConstZZpE(3);
  f_table[3] = alpha + testutil::ConstZZpE(1);

  vec_ZZ_pE g_table;
  g_table.SetLength(4);
  g_table[0] = alpha + testutil::ConstZZpE(2);
  g_table[1] = testutil::ConstZZpE(1);
  g_table[2] = alpha;
  g_table[3] = testutil::ConstZZpE(0);

  const ZZ_pE initial_claim = SumOfPointwiseProducts(f_table, g_table);
  basefold::ProductSumcheckProver prover(f_table, g_table);

  std::vector<basefold::QuadraticPoly> h_by_level(2);
  h_by_level[1] = prover.CurrentPolynomial();

  std::vector<ZZ_pE> r_by_level = {alpha + testutil::ConstZZpE(1),
                                   testutil::ConstZZpE(3)};
  prover.ReceiveChallenge(r_by_level[1]);
  h_by_level[0] = prover.CurrentPolynomial();
  prover.ReceiveChallenge(r_by_level[0]);

  CHECK(basefold::CheckProductSumcheckChain(initial_claim, h_by_level,
                                            r_by_level));

  const vec_ZZ_pE f_monomial =
      basefold::BooleanHypercubeTableToMonomialCoeffs(f_table);
  const vec_ZZ_pE g_monomial =
      basefold::BooleanHypercubeTableToMonomialCoeffs(g_table);
  const ZZ_pE expected_final =
      basefold::EvalMultilinearMonomialCoeffs(f_monomial, r_by_level) *
      basefold::EvalMultilinearMonomialCoeffs(g_monomial, r_by_level);
  CHECK_EQ(h_by_level[0].Eval(r_by_level[0]), expected_final);
}

void TestProductSumcheckChain_RejectsTamperedPolynomial() {
  testutil::PrintInfo("Ring-switch WP2: product sumcheck chain rejects tampering");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  vec_ZZ_pE f_table;
  f_table.SetLength(4);
  f_table[0] = testutil::ConstZZpE(0);
  f_table[1] = testutil::ConstZZpE(1);
  f_table[2] = alpha;
  f_table[3] = testutil::ConstZZpE(2);

  vec_ZZ_pE g_table;
  g_table.SetLength(4);
  g_table[0] = alpha + testutil::ConstZZpE(1);
  g_table[1] = testutil::ConstZZpE(3);
  g_table[2] = testutil::ConstZZpE(1);
  g_table[3] = alpha;

  const ZZ_pE initial_claim = SumOfPointwiseProducts(f_table, g_table);
  basefold::ProductSumcheckProver prover(f_table, g_table);

  std::vector<basefold::QuadraticPoly> h_by_level(2);
  h_by_level[1] = prover.CurrentPolynomial();

  const std::vector<ZZ_pE> r_by_level = {alpha, testutil::ConstZZpE(1)};
  prover.ReceiveChallenge(r_by_level[1]);
  h_by_level[0] = prover.CurrentPolynomial();

  CHECK(basefold::CheckProductSumcheckChain(initial_claim, h_by_level,
                                            r_by_level));

  std::vector<basefold::QuadraticPoly> tampered = h_by_level;
  tampered[0].a0 += testutil::ConstZZpE(1);
  CHECK(!basefold::CheckProductSumcheckChain(initial_claim, tampered,
                                             r_by_level));
}

void TestProductSumcheckProver_FromMonomialCoeffsMatchesTables() {
  testutil::PrintInfo("Ring-switch WP2: monomial helper matches table-native product sumcheck");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(4);
  f_coeffs[0] = testutil::ConstZZpE(1);
  f_coeffs[1] = alpha;
  f_coeffs[2] = testutil::ConstZZpE(2);
  f_coeffs[3] = alpha + testutil::ConstZZpE(1);

  vec_ZZ_pE g_coeffs;
  g_coeffs.SetLength(4);
  g_coeffs[0] = alpha + testutil::ConstZZpE(1);
  g_coeffs[1] = testutil::ConstZZpE(3);
  g_coeffs[2] = testutil::ConstZZpE(0);
  g_coeffs[3] = testutil::ConstZZpE(1);

  vec_ZZ_pE f_table;
  vec_ZZ_pE g_table;
  f_table.SetLength(4);
  g_table.SetLength(4);
  for (long i = 0; i < 4; ++i) {
    const std::vector<ZZ_pE> point = BooleanPointFromIndex(i, 2);
    f_table[i] = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, point);
    g_table[i] = basefold::EvalMultilinearMonomialCoeffs(g_coeffs, point);
  }

  basefold::ProductSumcheckProver from_tables(f_table, g_table);
  basefold::ProductSumcheckProver from_coeffs =
      basefold::ProductSumcheckProver::FromMonomialCoeffs(f_coeffs, g_coeffs);

  CHECK_EQ(from_tables.CurrentPolynomial().a0, from_coeffs.CurrentPolynomial().a0);
  CHECK_EQ(from_tables.CurrentPolynomial().a1, from_coeffs.CurrentPolynomial().a1);
  CHECK_EQ(from_tables.CurrentPolynomial().a2, from_coeffs.CurrentPolynomial().a2);

  const std::vector<ZZ_pE> r_by_level = {alpha + testutil::ConstZZpE(1),
                                         testutil::ConstZZpE(2)};
  from_tables.ReceiveChallenge(r_by_level[1]);
  from_coeffs.ReceiveChallenge(r_by_level[1]);

  CHECK_EQ(from_tables.CurrentPolynomial().a0, from_coeffs.CurrentPolynomial().a0);
  CHECK_EQ(from_tables.CurrentPolynomial().a1, from_coeffs.CurrentPolynomial().a1);
  CHECK_EQ(from_tables.CurrentPolynomial().a2, from_coeffs.CurrentPolynomial().a2);

  from_tables.ReceiveChallenge(r_by_level[0]);
  from_coeffs.ReceiveChallenge(r_by_level[0]);
}

void TestProductSumcheckProver_DimensionZeroUsesNoRounds() {
  testutil::PrintInfo("Ring-switch WP2: product sumcheck has an explicit d=0 no-round contract");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  vec_ZZ_pE f_table;
  f_table.SetLength(1);
  f_table[0] = alpha + testutil::ConstZZpE(1);

  vec_ZZ_pE g_table;
  g_table.SetLength(1);
  g_table[0] = alpha;

  basefold::ProductSumcheckProver prover(f_table, g_table);
  CHECK_EQ(prover.Dimension(), 0L);
  CHECK_EQ(prover.RemainingVars(), 0L);

  std::vector<basefold::QuadraticPoly> h_by_level;
  std::vector<ZZ_pE> r_by_level;
  CHECK(basefold::CheckProductSumcheckChain(
      f_table[0] * g_table[0], h_by_level, r_by_level));
}

}  // namespace

int main() {
  try {
    RUN_TEST(TestBaseFoldBackendAdapter_Smoke);
    RUN_TEST(TestBaseFoldBackendAdapter_RejectsSameDegreeContextSwitch);
    RUN_TEST(TestBaseFoldBackendAdapter_RejectsArtifactsFromDifferentHandle);
    RUN_TEST(TestBaseFoldBackendAdapter_RejectsProofReuseAcrossHandles);
    RUN_TEST(TestBaseFoldBackendAdapter_RejectsProofSizeAfterContextSwitch);
    RUN_TEST(TestRingSwitchSetup_DefaultBasesUseCurrentDegree);
    RUN_TEST(TestValidateCurrentZ2kRingContext_SucceedsForGR42);
    RUN_TEST(TestValidateCurrentZ2kRingContext_RejectsReducibleMod2Polynomial);
    RUN_TEST(TestRingSwitchSetup_Succeeds);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisAcceptsValidNonPolynomialBases);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisDerivesMissingDualBases);
    RUN_TEST(TestRingSwitchWP3Helpers_BuildReusableContextAndBasisCases);
    RUN_TEST(TestRingSwitchSetup_RejectsMismatchedBaseModulus);
    RUN_TEST(TestRingSwitchSetup_RejectsMismatchedExtensionModulus);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisRejectsMissingAlphaOrBeta);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisRejectsMissingAlpha_GR44);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisRejectsWrongDimension);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisRejectsWrongBetaDimension_GR44);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisRejectsBrokenDualBasis);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisRejectsBrokenBetaDualBasis_GR44);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisRejectsNonBasisAlpha);
    RUN_TEST(TestRingSwitchSetup_ProvidedBasisRejectsNonBasisBeta);
    RUN_TEST(TestRingSwitchSetup_RejectsBackendDimensionMismatch);
    RUN_TEST(TestPackZ2kCoeffsToGREvals_RoundTripsSmallExample);
    RUN_TEST(TestPackZ2kCoeffsToGREvals_ComposesAgainstProvidedBetaBasis);
    RUN_TEST(TestPackZ2kCoeffsToGREvals_ComposesAgainstProvidedBetaBasis_GR44);
    RUN_TEST(TestBooleanHypercubeTableToMonomialCoeffs_EvaluatesTableCorrectly);
    RUN_TEST(TestPackZ2kCoeffsToGREvals_RejectsNonBaseRingInputs);
    RUN_TEST(TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValues);
    RUN_TEST(
        TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValuesWithProvidedAlphaBasis);
    RUN_TEST(
        TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValuesWithProvidedAlphaBasis_GR44);
    RUN_TEST(TestBuildRingSwitchComponentTensor_RecoversPartialEvaluations);
    RUN_TEST(
        TestBuildRingSwitchComponentTensor_RecoversPartialEvaluationsWithIndependentAlphaBeta);
    RUN_TEST(
        TestBuildRingSwitchComponentTensor_RecoversPartialEvaluationsWithIndependentAlphaBeta_GR44);
    RUN_TEST(TestBuildRingSwitchComponentTensor_RCoeffsEvaluateAsExpected);
    RUN_TEST(TestBuildRingSwitchComponentTensor_RejectsWrongSuffixDimension);
    RUN_TEST(TestRingSwitchCommit_MatchesDirectBackendCommit);
    RUN_TEST(TestRingSwitchBuildCommitArtifacts_CachesPackedRepresentations);
    RUN_TEST(TestRingSwitchProveEvalFromCommitArtifacts_HonestProofIsSelfConsistent);
    RUN_TEST(TestRingSwitchProveEval_DirectAndArtifactPathsAgreeOnOuterMessages);
    RUN_TEST(TestRingSwitchVerifyEval_AcceptsHonestProof);
    RUN_TEST(TestRingSwitchVerifyEval_AcceptsHonestProofFromDirectProvePath);
    RUN_TEST(TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedAlphaBetaBases);
    RUN_TEST(
        TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedNonPolynomialAlphaOnly);
    RUN_TEST(
        TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedNonPolynomialBetaOnly);
    RUN_TEST(
        TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedExplicitDualBases);
    RUN_TEST(
        TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedAlphaBetaBasesOnVariantBackend);
    RUN_TEST(
        TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedNonPolynomialAlphaOnly_GR44);
    RUN_TEST(
        TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedNonPolynomialBetaOnly_GR44);
    RUN_TEST(
        TestRingSwitchVerifyEval_AcceptsHonestProofWithProvidedAlphaBetaBases_GR44);
    RUN_TEST(TestRingSwitchPaperAPI_AcceptsHonestProof);
    RUN_TEST(TestRingSwitchVerifyEval_RejectsTampering);
    RUN_TEST(TestRingSwitchVerifyEval_DimensionZeroUsesNoSumcheckRounds);
    RUN_TEST(TestRingSwitchOuterProveVerify_AcceptsHonestProof);
    RUN_TEST(TestRingSwitchOuterProveFromCommitArtifacts_CheckedAndUncheckedAgree);
    RUN_TEST(TestRingSwitchProveEvalFromCommitArtifacts_CheckedAndUncheckedAgree);
    RUN_TEST(TestRingSwitchProofSerialize_ComposedSizeMatchesBytes);
    RUN_TEST(
        TestRingSwitchProofSerialize_ComposedSizeMatchesBytesWithProvidedAlphaBeta);
    RUN_TEST(
        TestRingSwitchProofSerialize_ComposedSizeMatchesBytesWithProvidedAlphaBeta_GR44);
    RUN_TEST(TestProductSumcheckProver_BooleanTablesPasses);
    RUN_TEST(TestProductSumcheckChain_RejectsTamperedPolynomial);
    RUN_TEST(TestProductSumcheckProver_FromMonomialCoeffsMatchesTables);
    RUN_TEST(TestProductSumcheckProver_DimensionZeroUsesNoRounds);
  } catch (const exception &e) {
    cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    cerr << "Unhandled non-std exception\n";
    return 2;
  }

  if (g_test_failure_count == 0) {
    cout << "\nAll tests passed.\n";
    return 0;
  }

  cerr << "\n" << g_test_failure_count << " test(s) failed.\n";
  return 1;
}
