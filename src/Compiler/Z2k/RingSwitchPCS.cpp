#include "Compiler/Z2k/RingSwitchPCS.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pXFactoring.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <string>

#include "GaloisRing/FrobeniusBasis.hpp"
#include "PCS/Common/Hash.hpp"
#include "PCS/Common/Multilinear.hpp"
#include "PCS/Common/Profile.hpp"
#include "PCS/Common/Transcript.hpp"

using NTL::LogicError;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEBak;
using NTL::ZZ_pX;
using NTL::ZZ_pBak;
using NTL::vec_ZZ_pE;

namespace basefold {
namespace {

long Pow2LongOrThrow(long exponent, const char *what) {
  if (exponent < 0) {
    LogicError(what);
  }
  if (exponent >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError(what);
  }
  return 1L << exponent;
}

bool IsPowerOfTwoZZ(const ZZ &value) {
  if (value <= 0) {
    return false;
  }
  ZZ t = value;
  while ((t % 2) == 0) {
    t /= 2;
  }
  return t == 1;
}

ZZ NormalizeModNonNegative(const ZZ &value, const ZZ &modulus) {
  ZZ out = value % modulus;
  if (out < 0) {
    out += modulus;
  }
  return out;
}

ZZ_pX ReduceZZpXModPrime(const ZZ_pX &poly_over_pk, const ZZ &prime) {
  ZZ_pX out;
  NTL::clear(out);
  const long degree = NTL::deg(poly_over_pk);
  for (long i = 0; i <= degree; ++i) {
    ZZ_p coeff_mod_prime;
    NTL::conv(coeff_mod_prime,
              NormalizeModNonNegative(NTL::rep(NTL::coeff(poly_over_pk, i)),
                                      prime));
    if (coeff_mod_prime != 0) {
      NTL::SetCoeff(out, i, coeff_mod_prime);
    }
  }
  out.normalize();
  return out;
}

bool IsPowerOfTwoLong(long value) {
  return value > 0 && (value & (value - 1)) == 0;
}

long Log2ExactPowerOfTwoLongOrThrow(long value, const char *what) {
  if (!IsPowerOfTwoLong(value)) {
    LogicError(what);
  }
  long out = 0;
  while (value > 1) {
    value >>= 1;
    ++out;
  }
  return out;
}

ZZ_pE BaseRingConstant(const ZZ_p &value) {
  ZZ_pX poly;
  if (value != 0) {
    NTL::SetCoeff(poly, 0, value);
  }
  ZZ_pE out;
  NTL::conv(out, poly);
  return out;
}

ZZ_pE BaseRingConstant(long value) {
  return BaseRingConstant(NTL::to_ZZ_p(value));
}

std::vector<ZZ_pE> BuildActivePolynomialBasisOrThrow(const char *func_name) {
  if (!ZZ_pE::initialized()) {
    LogicError((std::string(func_name) +
                ": ZZ_pE context must be initialized")
                   .c_str());
  }
  const long degree = ZZ_pE::degree();
  if (degree <= 0) {
    LogicError((std::string(func_name) +
                ": current ZZ_pE degree must be positive")
                   .c_str());
  }

  std::vector<ZZ_pE> basis(static_cast<std::size_t>(degree));
  for (long i = 0; i < degree; ++i) {
    ZZ_pX poly;
    NTL::SetCoeff(poly, i, NTL::to_ZZ_p(1));
    NTL::conv(basis[static_cast<std::size_t>(i)], poly);
  }
  return basis;
}

GaloisRingBasisData BuildActivePolynomialBasisDataOrThrow(
    const char *func_name) {
  GaloisRingBasisData basis_data;
  basis_data.basis = BuildActivePolynomialBasisOrThrow(func_name);
  basis_data.dual_basis = ::BuildDualBasisOrThrow(basis_data.basis);
  return basis_data;
}

long BasisDimensionOrThrow(const GaloisRingBasisData &basis, const char *label,
                           const char *func_name) {
  const long dimension = static_cast<long>(basis.basis.size());
  if (dimension <= 0) {
    LogicError((std::string(func_name) + ": " + label +
                ".basis must be non-empty")
                   .c_str());
  }
  return dimension;
}

bool IsBaseRingConstant(const ZZ_pE &value) {
  const ZZ_pX poly = NTL::rep(value);
  const long degree = NTL::deg(poly);
  for (long i = 1; i <= degree; ++i) {
    if (NTL::coeff(poly, i) != 0) {
      return false;
    }
  }
  return true;
}

void ValidateBaseRingConstantOrThrow(const ZZ_pE &value, const char *label,
                                     long index, const char *func_name) {
  if (!IsBaseRingConstant(value)) {
    LogicError((std::string(func_name) + ": " + label + "[" +
                std::to_string(index) + "] must be a base-ring constant")
                   .c_str());
  }
}

void ValidateBaseRingVectorOrThrow(const vec_ZZ_pE &values, const char *label,
                                   const char *func_name) {
  for (long i = 0; i < values.length(); ++i) {
    ValidateBaseRingConstantOrThrow(values[i], label, i, func_name);
  }
}

std::vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  std::vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] = BaseRingConstant((index >> i) & 1L);
  }
  return point;
}

long CheckedMultiplyLong(long a, long b, const char *what) {
  if (a < 0 || b < 0) {
    LogicError(what);
  }
  if (a == 0 || b == 0) {
    return 0;
  }
  if (a > std::numeric_limits<long>::max() / b) {
    LogicError(what);
  }
  return a * b;
}

vec_ZZ_pE DecomposeGRElementToBaseCoeffsPolynomialBasisUnchecked(
    long basis_dimension, const ZZ_pE &element) {
  vec_ZZ_pE coeffs;
  coeffs.SetLength(basis_dimension);
  const ZZ_pX poly = NTL::rep(element);
  for (long i = 0; i < basis_dimension; ++i) {
    coeffs[i] = BaseRingConstant(NTL::coeff(poly, i));
  }
  return coeffs;
}

struct PackedCommitInputs {
  vec_ZZ_pE t_packed_table;
  vec_ZZ_pE t_packed_monomial_coeffs;
};

void ValidateActivePolynomialBasisDataOrThrow(const GaloisRingBasisData &basis,
                                              long expected_dimension,
                                              const char *label,
                                              const char *func_name) {
  if (static_cast<long>(basis.basis.size()) != expected_dimension) {
    LogicError((std::string(func_name) + ": " + label +
                ".basis size must match current ZZ_pE degree")
                   .c_str());
  }
  if (!basis.dual_basis.empty() &&
      static_cast<long>(basis.dual_basis.size()) != expected_dimension) {
    LogicError((std::string(func_name) + ": " + label +
                ".dual_basis size must match current ZZ_pE degree")
                   .c_str());
  }

  const std::vector<ZZ_pE> expected_basis =
      BuildActivePolynomialBasisOrThrow(func_name);
  for (long i = 0; i < expected_dimension; ++i) {
    if (basis.basis[static_cast<std::size_t>(i)] !=
        expected_basis[static_cast<std::size_t>(i)]) {
      LogicError((std::string(func_name) + ": " + label +
                  " must equal the active polynomial basis while general-basis support is pending")
                     .c_str());
    }
  }

  if (basis.dual_basis.empty()) {
    return;
  }
  const std::vector<ZZ_pE> expected_dual = ::BuildDualBasisOrThrow(expected_basis);
  for (long i = 0; i < expected_dimension; ++i) {
    if (basis.dual_basis[static_cast<std::size_t>(i)] !=
        expected_dual[static_cast<std::size_t>(i)]) {
      LogicError((std::string(func_name) + ": " + label +
                  ".dual_basis must equal the dual of the active polynomial basis while general-basis support is pending")
                     .c_str());
    }
  }
}

void ValidateBasicIrreducibilityModTwoOrThrow(const ZZ_pX &extension_modulus,
                                              long expected_degree,
                                              const char *func_name) {
  ZZ_pBak modulus_bak;
  modulus_bak.save();
  ZZ_pEBak extension_bak;
  extension_bak.save();

  const ZZ two(2);
  ZZ_p::init(two);
  const ZZ_pX reduced = ReduceZZpXModPrime(extension_modulus, two);
  if (NTL::deg(reduced) != expected_degree) {
    LogicError((std::string(func_name) +
                ": extension_modulus must stay full-degree after mod-2 reduction")
                   .c_str());
  }
  if (NTL::IterIrredTest(reduced) != 1) {
    LogicError((std::string(func_name) +
                ": extension_modulus must be basic irreducible modulo 2")
                   .c_str());
  }
}

PackedCommitInputs BuildPackedCommitInputs(const RingSwitchPCSParams &params,
                                           const vec_ZZ_pE &t_table) {
  PackedCommitInputs out;
  out.t_packed_table = PackZ2kCoeffsToGREvals(params, t_table);
  out.t_packed_monomial_coeffs =
      BooleanHypercubeTableToMonomialCoeffs(out.t_packed_table);
  return out;
}

HashTranscript MakeRingSwitchTranscript() {
  HashTranscriptConfig config;
  config.domain_separator = "RingSwitchPCS/v1";
  config.byte_order = TranscriptByteOrder::kLittleEndian;
  config.error_prefix = "RingSwitchHashTranscript";
  return HashTranscript(config);
}

void AbsorbPublicInput(HashTranscript &transcript,
                       const MerkleRoot &commitment,
                       const std::vector<FieldElement> &z,
                       const FieldElement &claimed_s) {
  transcript.AbsorbDigest(commitment);
  for (const FieldElement &zi : z) {
    transcript.AbsorbFieldElement(zi);
  }
  transcript.AbsorbFieldElement(claimed_s);
}

void ValidateEvalInputsOrThrow(const RingSwitchPCSParams &params,
                               const vec_ZZ_pE &t_table,
                               const std::vector<FieldElement> &z,
                               long num_queries, const char *func_name) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  const long expected_t_length = Pow2LongOrThrow(
      params.ell, (std::string(func_name) + ": ell is too large for long").c_str());
  if (t_table.length() != expected_t_length) {
    LogicError((std::string(func_name) + ": t_table length must equal 2^ell")
                   .c_str());
  }
  ValidateBaseRingVectorOrThrow(t_table, "t_table", func_name);
  if (static_cast<long>(z.size()) != params.ell) {
    LogicError((std::string(func_name) + ": z dimension must equal ell").c_str());
  }
  if (num_queries < 0) {
    LogicError((std::string(func_name) + ": num_queries must be non-negative")
                   .c_str());
  }
}

void ValidateCommitArtifactsOrThrow(const RingSwitchPCSParams &params,
                                    const RingSwitchPCSCommitArtifacts &artifacts,
                                    const char *func_name) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  const long expected_packed_length = Pow2LongOrThrow(
      params.ell_prime,
      (std::string(func_name) + ": ell_prime is too large for long").c_str());
  if (artifacts.t_packed_table.length() != expected_packed_length) {
    LogicError((std::string(func_name) +
                ": t_packed_table length must equal 2^(ell-kappa)")
                   .c_str());
  }
  if (artifacts.t_packed_monomial_coeffs.length() != expected_packed_length) {
    LogicError((std::string(func_name) +
                ": t_packed_monomial_coeffs length must equal 2^(ell-kappa)")
                   .c_str());
  }
  if (artifacts.commitment != artifacts.backend_commit_artifacts.commitment) {
    LogicError((std::string(func_name) +
                ": commitment must match backend_commit_artifacts.commitment")
                   .c_str());
  }
}

void ValidateOuterCommitArtifactsOrThrow(
    const RingSwitchPCSParams &params,
    const RingSwitchPCSOuterCommitArtifacts &artifacts,
    const char *func_name) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  const long expected_packed_length = Pow2LongOrThrow(
      params.ell_prime,
      (std::string(func_name) + ": ell_prime is too large for long").c_str());
  if (artifacts.t_packed_table.length() != expected_packed_length) {
    LogicError((std::string(func_name) +
                ": t_packed_table length must equal 2^(ell-kappa)")
                   .c_str());
  }
  if (artifacts.t_packed_monomial_coeffs.length() != expected_packed_length) {
    LogicError((std::string(func_name) +
                ": t_packed_monomial_coeffs length must equal 2^(ell-kappa)")
                   .c_str());
  }
}

bool HasExpectedEvalProofShape(const RingSwitchPCSParams &params,
                               const RingSwitchPCSEvalProof &proof) {
  return static_cast<long>(proof.s_by_u.size()) ==
             static_cast<long>(params.beta_basis.basis.size()) &&
         static_cast<long>(proof.h_by_level.size()) == params.ell_prime;
}

bool HasExpectedOuterEvalProofShape(const RingSwitchPCSParams &params,
                                    const RingSwitchPCSOuterEvalProof &proof) {
  return static_cast<long>(proof.s_by_u.size()) ==
             static_cast<long>(params.beta_basis.basis.size()) &&
         static_cast<long>(proof.h_by_level.size()) == params.ell_prime;
}

bool HasCompatibleBackendEvalSubproof(const RingSwitchPCSParams &params,
                                      const Z2kPCSBackendEvalProof &proof) {
  return proof.vtable == params.backend.vtable && proof.payload &&
         proof.params_owner &&
         proof.params_owner.get() == params.backend.params.get();
}

std::vector<FieldElement> SlicePoint(const std::vector<FieldElement> &z,
                                     long begin, long count) {
  return std::vector<FieldElement>(z.begin() + begin, z.begin() + begin + count);
}

vec_ZZ_pE BuildEqualityTable(const std::vector<FieldElement> &point) {
  const long dimension = static_cast<long>(point.size());
  const long length = Pow2LongOrThrow(
      dimension, "BuildEqualityTable: dimension is too large for long");
  vec_ZZ_pE table;
  table.SetLength(length);
  for (long idx = 0; idx < length; ++idx) {
    table[idx] = EqPolynomial(point, BooleanPointFromIndex(idx, dimension));
  }
  return table;
}

std::vector<FieldElement> RecoverPartialEvaluationsFromSByU(
    const RingSwitchPCSParams &params,
    const std::vector<FieldElement> &s_by_u) {
  const long basis_dimension =
      BasisDimensionOrThrow(params.alpha_basis, "alpha_basis",
                            "RecoverPartialEvaluationsFromSByU");
  if (static_cast<long>(s_by_u.size()) != basis_dimension) {
    LogicError(
        "RecoverPartialEvaluationsFromSByU: s_by_u size must equal basis dimension");
  }

  std::vector<FieldVec> coeff_rows(static_cast<std::size_t>(basis_dimension));
  for (long u = 0; u < basis_dimension; ++u) {
    coeff_rows[static_cast<std::size_t>(u)] =
        DecomposeGRElementToBaseCoeffsPolynomialBasis(params,
                                                      s_by_u[static_cast<std::size_t>(u)]);
  }

  std::vector<FieldElement> partials(static_cast<std::size_t>(basis_dimension),
                                     FieldElement(0));
  for (long v = 0; v < basis_dimension; ++v) {
    FieldElement acc = FieldElement(0);
    for (long u = 0; u < basis_dimension; ++u) {
      acc += coeff_rows[static_cast<std::size_t>(u)][v] *
             params.alpha_basis.basis[static_cast<std::size_t>(u)];
    }
    partials[static_cast<std::size_t>(v)] = acc;
  }
  return partials;
}

std::vector<FieldElement> ComputeDirectPartialEvaluations(
    const RingSwitchPCSParams &params, const vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z_suffix) {
  const long basis_dimension =
      BasisDimensionOrThrow(params.alpha_basis, "alpha_basis",
                            "ComputeDirectPartialEvaluations");
  const long num_w = Pow2LongOrThrow(
      params.ell_prime,
      "ComputeDirectPartialEvaluations: ell_prime is too large for long");
  std::vector<FieldElement> partials(static_cast<std::size_t>(basis_dimension),
                                     FieldElement(0));
  for (long v = 0; v < basis_dimension; ++v) {
    vec_ZZ_pE slice;
    slice.SetLength(num_w);
    for (long w = 0; w < num_w; ++w) {
      slice[w] = t_table[v + (w << params.kappa)];
    }
    const vec_ZZ_pE slice_monomial = BooleanHypercubeTableToMonomialCoeffs(slice);
    partials[static_cast<std::size_t>(v)] =
        EvalMultilinearMonomialCoeffs(slice_monomial, z_suffix);
  }
  return partials;
}

FieldElement RecombineClaimFromPartials(const std::vector<FieldElement> &partials,
                                        const std::vector<FieldElement> &z_prefix) {
  FieldElement acc = FieldElement(0);
  for (long v = 0; v < static_cast<long>(partials.size()); ++v) {
    acc += partials[static_cast<std::size_t>(v)] *
           EqPolynomial(z_prefix, BooleanPointFromIndex(v,
                                                        static_cast<long>(z_prefix.size())));
  }
  return acc;
}

std::vector<FieldElement> ComputeSByU(const RingSwitchComponentTensor &tensor,
                                      const vec_ZZ_pE &t_packed_table) {
  const long num_u = tensor.basis_dimension;
  const long num_w = Pow2LongOrThrow(
      tensor.ell_prime, "ComputeSByU: ell_prime is too large for long");
  if (t_packed_table.length() != num_w) {
    LogicError("ComputeSByU: t_packed_table length mismatch");
  }

  std::vector<FieldElement> s_by_u(static_cast<std::size_t>(num_u),
                                   FieldElement(0));
  for (long u = 0; u < num_u; ++u) {
    FieldElement acc = FieldElement(0);
    for (long w = 0; w < num_w; ++w) {
      acc += tensor.a_by_u_then_w[u * num_w + w] * t_packed_table[w];
    }
    s_by_u[static_cast<std::size_t>(u)] = acc;
  }
  return s_by_u;
}

vec_ZZ_pE BuildBatchedGTable(const RingSwitchComponentTensor &tensor,
                             const std::vector<FieldElement> &rprime_prefix) {
  const long num_u = tensor.basis_dimension;
  const long num_w = Pow2LongOrThrow(
      tensor.ell_prime, "BuildBatchedGTable: ell_prime is too large for long");
  vec_ZZ_pE out;
  out.SetLength(num_w);

  const vec_ZZ_pE eq_prefix = BuildEqualityTable(rprime_prefix);
  if (eq_prefix.length() != num_u) {
    LogicError("BuildBatchedGTable: prefix equality table length mismatch");
  }

  for (long w = 0; w < num_w; ++w) {
    FieldElement acc = FieldElement(0);
    for (long u = 0; u < num_u; ++u) {
      acc += tensor.a_by_u_then_w[u * num_w + w] * eq_prefix[u];
    }
    out[w] = acc;
  }
  return out;
}

FieldElement ComputeInitialBatchedClaim(const std::vector<FieldElement> &s_by_u,
                                        const std::vector<FieldElement> &rprime_prefix) {
  const vec_ZZ_pE eq_prefix = BuildEqualityTable(rprime_prefix);
  if (eq_prefix.length() != static_cast<long>(s_by_u.size())) {
    LogicError(
        "ComputeInitialBatchedClaim: prefix equality table length mismatch");
  }

  FieldElement acc = FieldElement(0);
  for (long u = 0; u < static_cast<long>(s_by_u.size()); ++u) {
    acc += s_by_u[static_cast<std::size_t>(u)] * eq_prefix[u];
  }
  return acc;
}

FieldElement ComputeOriginalEvaluation(const vec_ZZ_pE &t_table,
                                       const std::vector<FieldElement> &z) {
  const vec_ZZ_pE t_monomial = BooleanHypercubeTableToMonomialCoeffs(t_table);
  return EvalMultilinearMonomialCoeffs(t_monomial, z);
}

struct OuterProveEvalResult {
  RingSwitchPCSOuterEvalProof proof;
  std::vector<FieldElement> rprime_suffix;
};

OuterProveEvalResult ProveOuterEvalFromCommitArtifactsInternal(
    const RingSwitchPCSParams &params, const vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries,
    const RingSwitchPCSOuterCommitArtifacts &commit_artifacts,
    const char *func_name) {
  ValidateEvalInputsOrThrow(params, t_table, z, num_queries, func_name);
  ValidateOuterCommitArtifactsOrThrow(params, commit_artifacts, func_name);

  const FieldElement direct_eval = ComputeOriginalEvaluation(t_table, z);
  if (direct_eval != claimed_s) {
    LogicError((std::string(func_name) + ": claimed_s must equal t(z)").c_str());
  }

  const std::vector<FieldElement> z_prefix = SlicePoint(z, 0, params.kappa);
  const std::vector<FieldElement> z_suffix =
      SlicePoint(z, params.kappa, params.ell_prime);
  const RingSwitchComponentTensor tensor =
      BuildRingSwitchComponentTensor(params, z_suffix);

  OuterProveEvalResult out;
  out.proof.s_by_u = ComputeSByU(tensor, commit_artifacts.t_packed_table);
  out.proof.h_by_level.resize(static_cast<std::size_t>(params.ell_prime));

  const std::vector<FieldElement> recovered_partials =
      RecoverPartialEvaluationsFromSByU(params, out.proof.s_by_u);
  const std::vector<FieldElement> direct_partials =
      ComputeDirectPartialEvaluations(params, t_table, z_suffix);
  if (recovered_partials != direct_partials) {
    LogicError((std::string(func_name) +
                ": recovered partial evaluations do not match Appendix C.1 reconstruction")
                   .c_str());
  }
  if (RecombineClaimFromPartials(recovered_partials, z_prefix) != claimed_s) {
    LogicError((std::string(func_name) +
                ": Equality Check 1 failed on honest witness")
                   .c_str());
  }

  HashTranscript transcript = MakeRingSwitchTranscript();
  AbsorbPublicInput(transcript, commitment, z, claimed_s);
  for (const FieldElement &s_u : out.proof.s_by_u) {
    transcript.AbsorbFieldElement(s_u);
  }

  std::vector<FieldElement> rprime_prefix(
      static_cast<std::size_t>(params.kappa));
  for (long i = 0; i < params.kappa; ++i) {
    rprime_prefix[static_cast<std::size_t>(i)] =
        transcript.ChallengeFieldElement("rprime/prefix/" + std::to_string(i));
  }

  const FieldElement initial_claim =
      ComputeInitialBatchedClaim(out.proof.s_by_u, rprime_prefix);
  const vec_ZZ_pE g_table = BuildBatchedGTable(tensor, rprime_prefix);

  out.rprime_suffix.resize(static_cast<std::size_t>(params.ell_prime));
  if (params.ell_prime > 0) {
    ProductSumcheckProver sumcheck(commit_artifacts.t_packed_table, g_table);
    out.proof.h_by_level[static_cast<std::size_t>(params.ell_prime - 1)] =
        sumcheck.CurrentPolynomial();
    AbsorbQuadraticPoly(
        transcript,
        out.proof.h_by_level[static_cast<std::size_t>(params.ell_prime - 1)]);

    for (long i = params.ell_prime; i-- > 0;) {
      const FieldElement r_i = transcript.ChallengeFieldElement(
          "rprime/suffix/" + std::to_string(i));
      out.rprime_suffix[static_cast<std::size_t>(i)] = r_i;
      sumcheck.ReceiveChallenge(r_i);
      if (i > 0) {
        out.proof.h_by_level[static_cast<std::size_t>(i - 1)] =
            sumcheck.CurrentPolynomial();
        AbsorbQuadraticPoly(
            transcript,
            out.proof.h_by_level[static_cast<std::size_t>(i - 1)]);
      }
    }

    if (!CheckProductSumcheckChain(initial_claim, out.proof.h_by_level,
                                   out.rprime_suffix)) {
      LogicError((std::string(func_name) +
                  ": honest product sumcheck chain is inconsistent")
                     .c_str());
    }
  }

  out.proof.t_star = EvalMultilinearMonomialCoeffs(
      commit_artifacts.t_packed_monomial_coeffs, out.rprime_suffix);

  std::vector<FieldElement> rprime_full = rprime_prefix;
  rprime_full.insert(rprime_full.end(), out.rprime_suffix.begin(),
                     out.rprime_suffix.end());
  const FieldElement g_star =
      EvalMultilinearMonomialCoeffs(tensor.r_monomial_coeffs, rprime_full);
  const FieldElement final_sumcheck_claim =
      (params.ell_prime == 0)
          ? initial_claim
          : out.proof.h_by_level[0].Eval(out.rprime_suffix[0]);
  if (final_sumcheck_claim != out.proof.t_star * g_star) {
    LogicError((std::string(func_name) +
                ": honest Equality Check 3 failed")
                   .c_str());
  }

  return out;
}

bool VerifyOuterEvalAndMaybeRecoverSuffix(
    const RingSwitchPCSParams &params, const MerkleRoot &commitment,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const RingSwitchPCSOuterEvalProof &proof,
    std::vector<FieldElement> *rprime_suffix_out) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  if (static_cast<long>(z.size()) != params.ell) {
    return false;
  }
  if (num_queries < 0) {
    return false;
  }
  if (!HasExpectedOuterEvalProofShape(params, proof)) {
    return false;
  }

  const std::vector<FieldElement> z_prefix = SlicePoint(z, 0, params.kappa);
  const std::vector<FieldElement> z_suffix =
      SlicePoint(z, params.kappa, params.ell_prime);
  const RingSwitchComponentTensor tensor =
      BuildRingSwitchComponentTensor(params, z_suffix);

  const std::vector<FieldElement> recovered_partials =
      RecoverPartialEvaluationsFromSByU(params, proof.s_by_u);
  if (RecombineClaimFromPartials(recovered_partials, z_prefix) != claimed_s) {
    return false;
  }

  HashTranscript transcript = MakeRingSwitchTranscript();
  AbsorbPublicInput(transcript, commitment, z, claimed_s);
  for (const FieldElement &s_u : proof.s_by_u) {
    transcript.AbsorbFieldElement(s_u);
  }

  std::vector<FieldElement> rprime_prefix(
      static_cast<std::size_t>(params.kappa));
  for (long i = 0; i < params.kappa; ++i) {
    rprime_prefix[static_cast<std::size_t>(i)] =
        transcript.ChallengeFieldElement("rprime/prefix/" + std::to_string(i));
  }

  const FieldElement initial_claim =
      ComputeInitialBatchedClaim(proof.s_by_u, rprime_prefix);

  std::vector<FieldElement> rprime_suffix(
      static_cast<std::size_t>(params.ell_prime));
  if (params.ell_prime > 0) {
    AbsorbQuadraticPoly(
        transcript,
        proof.h_by_level[static_cast<std::size_t>(params.ell_prime - 1)]);
    for (long i = params.ell_prime; i-- > 0;) {
      rprime_suffix[static_cast<std::size_t>(i)] = transcript.ChallengeFieldElement(
          "rprime/suffix/" + std::to_string(i));
      if (i > 0) {
        AbsorbQuadraticPoly(
            transcript,
            proof.h_by_level[static_cast<std::size_t>(i - 1)]);
      }
    }
  }

  if (!CheckProductSumcheckChain(initial_claim, proof.h_by_level,
                                 rprime_suffix)) {
    return false;
  }

  std::vector<FieldElement> rprime_full = rprime_prefix;
  rprime_full.insert(rprime_full.end(), rprime_suffix.begin(),
                     rprime_suffix.end());
  const FieldElement g_star =
      EvalMultilinearMonomialCoeffs(tensor.r_monomial_coeffs, rprime_full);
  const FieldElement final_sumcheck_claim =
      (params.ell_prime == 0)
          ? initial_claim
          : proof.h_by_level[0].Eval(rprime_suffix[0]);
  if (final_sumcheck_claim != proof.t_star * g_star) {
    return false;
  }

  if (rprime_suffix_out != nullptr) {
    *rprime_suffix_out = std::move(rprime_suffix);
  }
  return true;
}

}  // namespace

void ValidateCurrentZ2kRingContextOrThrow(const ZZ &base_modulus,
                                          const ZZ_pX &extension_modulus,
                                          long kappa) {
  if (kappa < 1) {
    LogicError("ValidateCurrentZ2kRingContextOrThrow: kappa must be >= 1");
  }
  if (base_modulus <= 1) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: base_modulus must be > 1");
  }
  if (!IsPowerOfTwoZZ(base_modulus)) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: base_modulus must be a power of two");
  }
  if (ZZ_p::modulus() != base_modulus) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: current ZZ_p modulus does not match base_modulus");
  }
  if (!ZZ_pE::initialized()) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: ZZ_pE context must be initialized");
  }

  const long expected_degree = Pow2LongOrThrow(
      kappa,
      "ValidateCurrentZ2kRingContextOrThrow: kappa is too large for long");
  if (ZZ_pE::degree() != expected_degree) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: current ZZ_pE degree must equal 2^kappa");
  }
  if (NTL::deg(extension_modulus) != expected_degree) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: extension_modulus degree must equal 2^kappa");
  }
  if (ZZ_pE::modulus().val() != extension_modulus) {
    LogicError(
        "ValidateCurrentZ2kRingContextOrThrow: current ZZ_pE modulus does not match extension_modulus");
  }
  ValidateBasicIrreducibilityModTwoOrThrow(
      extension_modulus, expected_degree,
      "ValidateCurrentZ2kRingContextOrThrow");
}

void ValidateRingSwitchPCSParamsOrThrow(const RingSwitchPCSParams &params) {
  if (params.kappa < 1) {
    LogicError("ValidateRingSwitchPCSParamsOrThrow: kappa must be >= 1");
  }
  if (params.ell < params.kappa) {
    LogicError("ValidateRingSwitchPCSParamsOrThrow: ell must be >= kappa");
  }

  ValidateCurrentZ2kRingContextOrThrow(params.base_modulus,
                                       params.extension_modulus, params.kappa);
  Z2kPCSBackendValidateParamsOrThrow(params.backend);

  const long expected_degree = ZZ_pE::degree();
  ValidateActivePolynomialBasisDataOrThrow(
      params.alpha_basis, expected_degree, "alpha_basis",
      "ValidateRingSwitchPCSParamsOrThrow");
  ValidateActivePolynomialBasisDataOrThrow(
      params.beta_basis, expected_degree, "beta_basis",
      "ValidateRingSwitchPCSParamsOrThrow");

  const long ell_prime = params.ell - params.kappa;
  if (params.ell_prime != ell_prime) {
    LogicError(
        "ValidateRingSwitchPCSParamsOrThrow: ell_prime must equal ell-kappa");
  }

  const long expected_backend_message_length = Pow2LongOrThrow(
      ell_prime,
      "ValidateRingSwitchPCSParamsOrThrow: ell-kappa is too large for long");
  if (Z2kPCSBackendMessageLength(params.backend) !=
      expected_backend_message_length) {
    LogicError(
        "ValidateRingSwitchPCSParamsOrThrow: backend message length must equal 2^(ell-kappa)");
  }
  if (Z2kPCSBackendPointDimension(params.backend) != ell_prime) {
    LogicError(
        "ValidateRingSwitchPCSParamsOrThrow: backend point dimension must equal ell-kappa");
  }
}

RingSwitchPCSParams RingSwitchPCSSetup(const RingSwitchPCSSetupInput &input) {
  RingSwitchPCSParams params;
  params.ell = input.ell;
  params.kappa = input.kappa;
  params.ell_prime = input.ell - input.kappa;
  params.base_modulus = input.base_modulus;
  params.extension_modulus = input.extension_modulus;
  if (input.use_provided_basis) {
    LogicError(
        "RingSwitchPCSSetup: caller-provided alpha/beta bases are not implemented yet; continue after WP2");
  }
  params.alpha_basis =
      BuildActivePolynomialBasisDataOrThrow("RingSwitchPCSSetup");
  params.beta_basis =
      BuildActivePolynomialBasisDataOrThrow("RingSwitchPCSSetup");
  params.backend = input.backend;
  ValidateRingSwitchPCSParamsOrThrow(params);
  return params;
}

vec_ZZ_pE BooleanHypercubeTableToMonomialCoeffs(const vec_ZZ_pE &table_values) {
  const long n = table_values.length();
  if (!IsPowerOfTwoLong(n)) {
    LogicError(
        "BooleanHypercubeTableToMonomialCoeffs: table length must be a power of two");
  }
  const long dimension = Log2ExactPowerOfTwoLongOrThrow(
      n,
      "BooleanHypercubeTableToMonomialCoeffs: table length must be a power of two");
  vec_ZZ_pE coeffs = table_values;
  for (long bit = 0; bit < dimension; ++bit) {
    const long step = 1L << bit;
    for (long mask = 0; mask < n; ++mask) {
      if (mask & step) {
        coeffs[mask] -= coeffs[mask ^ step];
      }
    }
  }
  return coeffs;
}

vec_ZZ_pE PackZ2kCoeffsToGREvals(const RingSwitchPCSParams &params,
                                 const vec_ZZ_pE &t_table) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  const long expected_length = Pow2LongOrThrow(
      params.ell, "PackZ2kCoeffsToGREvals: ell is too large for long");
  if (t_table.length() != expected_length) {
    LogicError("PackZ2kCoeffsToGREvals: t_table length must equal 2^ell");
  }
  ValidateBaseRingVectorOrThrow(t_table, "t_table",
                                "PackZ2kCoeffsToGREvals");

  const long basis_dimension =
      BasisDimensionOrThrow(params.beta_basis, "beta_basis",
                            "PackZ2kCoeffsToGREvals");
  const long packed_length = Pow2LongOrThrow(
      params.ell_prime,
      "PackZ2kCoeffsToGREvals: ell_prime is too large for long");
  vec_ZZ_pE packed;
  packed.SetLength(packed_length);

  for (long w = 0; w < packed_length; ++w) {
    ZZ_pX poly;
    for (long v = 0; v < basis_dimension; ++v) {
      const long coeff_index = v + w * basis_dimension;
      const ZZ_p constant_coeff = NTL::coeff(NTL::rep(t_table[coeff_index]), 0);
      if (constant_coeff != 0) {
        NTL::SetCoeff(poly, v, constant_coeff);
      }
    }
    NTL::conv(packed[w], poly);
  }

  return packed;
}

MerkleRoot RingSwitchPCSCommit(const RingSwitchPCSParams &params,
                               const vec_ZZ_pE &t_table) {
  return RingSwitchPCSBuildCommitArtifacts(params, t_table).commitment;
}

RingSwitchPCSOuterCommitArtifacts RingSwitchPCSBuildOuterCommitArtifacts(
    const RingSwitchPCSParams &params, const vec_ZZ_pE &t_table) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  const PackedCommitInputs packed = BuildPackedCommitInputs(params, t_table);

  RingSwitchPCSOuterCommitArtifacts out;
  out.t_packed_table = packed.t_packed_table;
  out.t_packed_monomial_coeffs = packed.t_packed_monomial_coeffs;
  return out;
}

RingSwitchPCSCommitArtifacts RingSwitchPCSBuildCommitArtifacts(
    const RingSwitchPCSParams &params, const vec_ZZ_pE &t_table) {
  const RingSwitchPCSOuterCommitArtifacts outer =
      RingSwitchPCSBuildOuterCommitArtifacts(params, t_table);
  RingSwitchPCSCommitArtifacts out;
  out.t_packed_table = outer.t_packed_table;
  out.t_packed_monomial_coeffs = outer.t_packed_monomial_coeffs;
  out.backend_commit_artifacts =
      Z2kPCSBackendBuildCommitArtifacts(params.backend,
                                        out.t_packed_monomial_coeffs);
  out.commitment = out.backend_commit_artifacts.commitment;
  return out;
}

RingSwitchPCSEvalProof RingSwitchPCSProveEval(
    const RingSwitchPCSParams &params, const vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries) {
  const RingSwitchPCSCommitArtifacts commit_artifacts =
      RingSwitchPCSBuildCommitArtifacts(params, t_table);
  return RingSwitchPCSProveEvalFromCommitArtifacts(
      params, t_table, z, claimed_s, num_queries, commit_artifacts);
}

RingSwitchPCSOuterEvalProof RingSwitchPCSProveOuterEval(
    const RingSwitchPCSParams &params, const vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries) {
  const RingSwitchPCSOuterCommitArtifacts commit_artifacts =
      RingSwitchPCSBuildOuterCommitArtifacts(params, t_table);
  return RingSwitchPCSProveOuterEvalFromCommitArtifacts(
      params, t_table, commitment, z, claimed_s, num_queries, commit_artifacts);
}

RingSwitchPCSOuterEvalProof RingSwitchPCSProveOuterEvalFromCommitArtifacts(
    const RingSwitchPCSParams &params, const vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries,
    const RingSwitchPCSOuterCommitArtifacts &commit_artifacts) {
  return ProveOuterEvalFromCommitArtifactsInternal(
             params, t_table, commitment, z, claimed_s, num_queries,
             commit_artifacts,
             "RingSwitchPCSProveOuterEvalFromCommitArtifacts")
      .proof;
}

RingSwitchPCSEvalProof RingSwitchPCSProveEvalFromCommitArtifacts(
    const RingSwitchPCSParams &params, const vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const RingSwitchPCSCommitArtifacts &commit_artifacts) {
  ValidateCommitArtifactsOrThrow(params, commit_artifacts,
                                 "RingSwitchPCSProveEvalFromCommitArtifacts");
  RingSwitchPCSOuterCommitArtifacts outer_commit_artifacts;
  outer_commit_artifacts.t_packed_table = commit_artifacts.t_packed_table;
  outer_commit_artifacts.t_packed_monomial_coeffs =
      commit_artifacts.t_packed_monomial_coeffs;
  OuterProveEvalResult outer = ProveOuterEvalFromCommitArtifactsInternal(
      params, t_table, commit_artifacts.commitment, z, claimed_s, num_queries,
      outer_commit_artifacts, "RingSwitchPCSProveEvalFromCommitArtifacts");

  RingSwitchPCSEvalProof proof;
  proof.s_by_u = outer.proof.s_by_u;
  proof.h_by_level = outer.proof.h_by_level;
  proof.t_star = outer.proof.t_star;

  {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->z2k_backend_prove_ns : nullptr,
                      prof ? &prof->z2k_backend_prove_calls : nullptr);
    proof.backend_proof = Z2kPCSBackendProveEval(
        params.backend, commit_artifacts.t_packed_monomial_coeffs,
        outer.rprime_suffix, proof.t_star, num_queries,
        &commit_artifacts.backend_commit_artifacts);
  }
  return proof;
}

bool RingSwitchPCSVerifyEval(const RingSwitchPCSParams &params,
                             const MerkleRoot &commitment,
                             const std::vector<FieldElement> &z,
                             const FieldElement &claimed_s, long num_queries,
                             const RingSwitchPCSEvalProof &proof) {
  if (!HasExpectedEvalProofShape(params, proof) ||
      !HasCompatibleBackendEvalSubproof(params, proof.backend_proof)) {
    return false;
  }

  RingSwitchPCSOuterEvalProof outer_proof;
  outer_proof.s_by_u = proof.s_by_u;
  outer_proof.h_by_level = proof.h_by_level;
  outer_proof.t_star = proof.t_star;

  std::vector<FieldElement> rprime_suffix;
  if (!VerifyOuterEvalAndMaybeRecoverSuffix(params, commitment, z, claimed_s,
                                            num_queries, outer_proof,
                                            &rprime_suffix)) {
    return false;
  }

  {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->z2k_backend_verify_ns : nullptr,
                      prof ? &prof->z2k_backend_verify_calls : nullptr);
    return Z2kPCSBackendVerifyEval(params.backend, commitment, rprime_suffix,
                                   proof.t_star, num_queries,
                                   proof.backend_proof);
  }
}

bool RingSwitchPCSVerifyOuterEval(const RingSwitchPCSParams &params,
                                  const MerkleRoot &commitment,
                                  const std::vector<FieldElement> &z,
                                  const FieldElement &claimed_s,
                                  long num_queries,
                                  const RingSwitchPCSOuterEvalProof &proof) {
  return VerifyOuterEvalAndMaybeRecoverSuffix(params, commitment, z, claimed_s,
                                              num_queries, proof, nullptr);
}

vec_ZZ_pE DecomposeGRElementToBaseCoeffsPolynomialBasis(
    const RingSwitchPCSParams &params, const ZZ_pE &element) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  return DecomposeGRElementToBaseCoeffsPolynomialBasisUnchecked(
      BasisDimensionOrThrow(params.alpha_basis, "alpha_basis",
                            "DecomposeGRElementToBaseCoeffsPolynomialBasis"),
      element);
}

vec_ZZ_pE DecomposeGRElementToBaseCoeffs(
    const RingSwitchPCSParams &params, const ZZ_pE &element,
    const GaloisRingBasisData &basis) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  ValidateActivePolynomialBasisDataOrThrow(
      basis, ZZ_pE::degree(), "basis", "DecomposeGRElementToBaseCoeffs");
  return DecomposeGRElementToBaseCoeffsPolynomialBasisUnchecked(
      BasisDimensionOrThrow(basis, "basis", "DecomposeGRElementToBaseCoeffs"),
      element);
}

RingSwitchComponentTensor BuildRingSwitchComponentTensor(
    const RingSwitchPCSParams &params, const std::vector<ZZ_pE> &r_suffix) {
  ValidateRingSwitchPCSParamsOrThrow(params);
  if (static_cast<long>(r_suffix.size()) != params.ell_prime) {
    LogicError(
        "BuildRingSwitchComponentTensor: r_suffix dimension must equal ell_prime");
  }

  const long basis_dimension =
      BasisDimensionOrThrow(params.alpha_basis, "alpha_basis",
                            "BuildRingSwitchComponentTensor");
  const long num_w = Pow2LongOrThrow(
      params.ell_prime,
      "BuildRingSwitchComponentTensor: ell_prime is too large for long");
  const long total_coeffs = CheckedMultiplyLong(
      basis_dimension, num_w,
      "BuildRingSwitchComponentTensor: coefficient table is too large for long");

  RingSwitchComponentTensor tensor;
  tensor.basis_dimension = basis_dimension;
  tensor.ell_prime = params.ell_prime;
  tensor.a_by_u_then_w.SetLength(total_coeffs);
  tensor.r_table.SetLength(total_coeffs);

  for (long w = 0; w < num_w; ++w) {
    const std::vector<ZZ_pE> bool_point =
        BooleanPointFromIndex(w, params.ell_prime);
    const ZZ_pE eq_at_w = EqPolynomial(r_suffix, bool_point);
    const vec_ZZ_pE coeffs =
        DecomposeGRElementToBaseCoeffsPolynomialBasisUnchecked(basis_dimension,
                                                               eq_at_w);
    for (long u = 0; u < basis_dimension; ++u) {
      const long row_major_index = u * num_w + w;
      const long flat_index = u + w * basis_dimension;
      tensor.a_by_u_then_w[row_major_index] = coeffs[u];
      tensor.r_table[flat_index] = coeffs[u];
    }
  }
  tensor.r_monomial_coeffs =
      BooleanHypercubeTableToMonomialCoeffs(tensor.r_table);

  return tensor;
}

}  // namespace basefold
