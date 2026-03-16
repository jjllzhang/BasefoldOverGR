#include "BaseFoldPCSInternal.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pEXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include "PCS/Common/MerkleMultiproofPlanner.hpp"
#include "PCS/Common/MerkleMultiproofReplay.hpp"
#include "PCS/Common/Multilinear.hpp"

using NTL::LogicError;
using NTL::NumBits;
using NTL::NumBytes;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEX;
using NTL::ZZ_pX;
using NTL::vec_ZZ_pE;

namespace basefold {
namespace {

using LocalMultiproofPlan = multiproof_planner::Plan;

void AppendU64(Bytes &out, std::uint64_t value) {
  for (int i = 7; i >= 0; --i) {
    out.push_back(static_cast<Byte>((value >> (8 * i)) & 0xff));
  }
}

void StoreU64BigEndian(Byte *dst, std::uint64_t value) {
  for (int i = 0; i < 8; ++i) {
    dst[i] = static_cast<Byte>((value >> (8 * (7 - i))) & 0xff);
  }
}

void AppendSerializedFieldElement(Bytes &out, const FieldElement &x,
                                  const char *func_name) {
  const long r = ZZ_pE::degree();
  if (r <= 0) {
    const std::string msg =
        std::string(func_name) + ": invalid extension degree";
    LogicError(msg.c_str());
  }

  const ZZ_pX &poly = NTL::rep(x);
  AppendU64(out, static_cast<std::uint64_t>(r));
  for (long i = 0; i < r; ++i) {
    const ZZ coeff = NTL::rep(NTL::coeff(poly, i));
    const long n = NumBytes(coeff);
    AppendU64(out, static_cast<std::uint64_t>(n));
    if (n > 0) {
      const std::size_t old_size = out.size();
      out.resize(old_size + static_cast<std::size_t>(n));
      NTL::BytesFromZZ(
          reinterpret_cast<unsigned char *>(out.data() + old_size), coeff, n);
    }
  }
}

Bytes SerializeFieldElement(const FieldElement &x) {
  Bytes out;
  AppendSerializedFieldElement(out, x, "SerializeFieldElement");
  return out;
}

ZZ NormalizeModNonNegative(const ZZ &a, const ZZ &m) {
  ZZ r = a % m;
  if (r < 0) {
    r += m;
  }
  return r;
}

bool IsPrimePowerOf(const ZZ &n, const ZZ &p) {
  if (n <= 1 || p <= 1) {
    return false;
  }
  ZZ t = n;
  long exp = 0;
  while (true) {
    const ZZ r = t % p;
    if (r != 0) {
      break;
    }
    t /= p;
    ++exp;
  }
  return exp > 0 && t == 1;
}

ZZ_pX ReduceZZpXModPrime(const ZZ_pX &poly_over_pk, const ZZ &p) {
  ZZ_pX out;
  NTL::clear(out);
  const long d = NTL::deg(poly_over_pk);
  for (long i = 0; i <= d; ++i) {
    ZZ_p c_mod_p;
    NTL::conv(c_mod_p,
              NormalizeModNonNegative(NTL::rep(NTL::coeff(poly_over_pk, i)), p));
    if (c_mod_p != 0) {
      NTL::SetCoeff(out, i, c_mod_p);
    }
  }
  out.normalize();
  return out;
}

ZZ_pEX ReduceExtensionPolynomialToResidueField(const ZZ_pEX &poly_over_pk_ext,
                                               const ZZ &p,
                                               long base_degree) {
  ZZ_pEX out;
  NTL::clear(out);
  const long d = NTL::deg(poly_over_pk_ext);
  for (long i = 0; i <= d; ++i) {
    const ZZ_pX coeff_poly_over_pk = NTL::rep(NTL::coeff(poly_over_pk_ext, i));
    ZZ_pX coeff_poly_mod_p;
    NTL::clear(coeff_poly_mod_p);
    for (long j = 0; j < base_degree; ++j) {
      ZZ_p c_mod_p;
      NTL::conv(
          c_mod_p,
          NormalizeModNonNegative(
              NTL::rep(NTL::coeff(coeff_poly_over_pk, j)), p));
      if (c_mod_p != 0) {
        NTL::SetCoeff(coeff_poly_mod_p, j, c_mod_p);
      }
    }
    coeff_poly_mod_p.normalize();
    ZZ_pE coeff_ext;
    NTL::conv(coeff_ext, coeff_poly_mod_p);
    NTL::SetCoeff(out, i, coeff_ext);
  }
  out.normalize();
  return out;
}

void ValidateChallengeConfigOrThrow(
    const BaseFoldPCSChallengeConfig &challenge_cfg,
    const FoldableCodeParams &params) {
  if (!challenge_cfg.use_extension_challenges) {
    return;
  }

  const long ext_degree = NTL::deg(challenge_cfg.challenge_extension_modulus);
  if (ext_degree <= 0) {
    LogicError(
        "ValidateChallengeConfigOrThrow: extension modulus must have positive degree");
  }

  FieldElement one;
  NTL::set(one);
  if (NTL::LeadCoeff(challenge_cfg.challenge_extension_modulus) != one) {
    LogicError("ValidateChallengeConfigOrThrow: extension modulus must be monic");
  }

  const ZZ modulus = ZZ_p::modulus();
  if (modulus <= 1) {
    LogicError("ValidateChallengeConfigOrThrow: invalid ZZ_p modulus");
  }

  const long base_degree = ZZ_pE::degree();
  if (base_degree <= 0) {
    LogicError("ValidateChallengeConfigOrThrow: invalid ZZ_pE degree");
  }

  ZZ base_prime = params.p;
  if (base_prime <= 1) {
    if (NTL::ProbPrime(modulus)) {
      base_prime = modulus;
    } else {
      LogicError(
          "ValidateChallengeConfigOrThrow: params.p must be set to the base prime in ring mode");
    }
  }
  if (!NTL::ProbPrime(base_prime)) {
    LogicError("ValidateChallengeConfigOrThrow: params.p must be prime");
  }
  if (!IsPrimePowerOf(modulus, base_prime)) {
    LogicError(
        "ValidateChallengeConfigOrThrow: current ZZ_p modulus must be a power of params.p");
  }

  if (NTL::ProbPrime(modulus)) {
    if (NTL::IterIrredTest(challenge_cfg.challenge_extension_modulus) != 1) {
      LogicError(
          "ValidateChallengeConfigOrThrow: extension modulus must be irreducible over the current base field");
    }
    return;
  }

  const ZZ_pX base_modulus_over_pk = ZZ_pE::modulus().val();

  NTL::ZZ_pBak modulus_bak;
  modulus_bak.save();
  NTL::ZZ_pEBak extension_bak;
  extension_bak.save();

  ZZ_p::init(base_prime);

  const ZZ_pX base_modulus_over_p =
      ReduceZZpXModPrime(base_modulus_over_pk, base_prime);
  if (NTL::deg(base_modulus_over_p) != base_degree) {
    LogicError(
        "ValidateChallengeConfigOrThrow: current ZZ_pE modulus must stay full-degree after mod-p reduction");
  }
  if (NTL::IterIrredTest(base_modulus_over_p) != 1) {
    LogicError(
        "ValidateChallengeConfigOrThrow: current ZZ_pE modulus must reduce to an irreducible polynomial modulo p");
  }
  ZZ_pE::init(base_modulus_over_p);

  const ZZ_pEX ext_modulus_over_p = ReduceExtensionPolynomialToResidueField(
      challenge_cfg.challenge_extension_modulus, base_prime, base_degree);
  if (NTL::deg(ext_modulus_over_p) != ext_degree) {
    LogicError(
        "ValidateChallengeConfigOrThrow: extension modulus must stay full-degree after mod-p reduction");
  }
  if (NTL::IterIrredTest(ext_modulus_over_p) != 1) {
    LogicError(
        "ValidateChallengeConfigOrThrow: extension modulus must be basic irreducible (irreducible after mod-p reduction)");
  }
}

Bytes SerializeExtensionPolynomial(const ZZ_pEX &poly) {
  const long d = NTL::deg(poly);
  Bytes out;
  const std::uint64_t coeff_count =
      (d < 0) ? 0ULL : static_cast<std::uint64_t>(d + 1);
  AppendU64(out, coeff_count);
  for (long i = 0; i <= d; ++i) {
    const Bytes coeff_bytes = SerializeFieldElement(NTL::coeff(poly, i));
    AppendU64(out, static_cast<std::uint64_t>(coeff_bytes.size()));
    out.insert(out.end(), coeff_bytes.begin(), coeff_bytes.end());
  }
  return out;
}

void AbsorbChallengeConfig(HashTranscript &transcript,
                           const BaseFoldPCSChallengeConfig &challenge_cfg) {
  Bytes tag;
  tag.push_back(
      static_cast<Byte>(challenge_cfg.use_extension_challenges ? 1 : 0));
  transcript.AbsorbBytes(tag);
  if (!challenge_cfg.use_extension_challenges) {
    return;
  }
  transcript.AbsorbBytes(
      SerializeExtensionPolynomial(challenge_cfg.challenge_extension_modulus));
}

long ExtensionDegreeOrThrow(const ZZ_pEX &extension_modulus,
                            const char *func_name) {
  const long ext_degree = NTL::deg(extension_modulus);
  if (ext_degree <= 0) {
    const std::string msg =
        std::string(func_name) + ": invalid extension degree";
    LogicError(msg.c_str());
  }
  return ext_degree;
}

const NTL::ZZ_pEXModulus &ExtensionModulusContextOrThrow(
    const ZZ_pEX &extension_modulus, const char *func_name) {
  (void)ExtensionDegreeOrThrow(extension_modulus, func_name);

  struct CachedModulusContext {
    bool initialized = false;
    ZZ base_modulus;
    long base_degree = 0;
    ZZ_pEX modulus_poly;
    NTL::ZZ_pEXModulus modulus_ctx;
  };

  thread_local CachedModulusContext cache;

  const ZZ &cur_base_modulus = ZZ_p::modulus();
  const long cur_base_degree = ZZ_pE::degree();
  if (!cache.initialized || cache.base_modulus != cur_base_modulus ||
      cache.base_degree != cur_base_degree ||
      cache.modulus_poly != extension_modulus) {
    cache.base_modulus = cur_base_modulus;
    cache.base_degree = cur_base_degree;
    cache.modulus_poly = extension_modulus;
    NTL::build(cache.modulus_ctx, cache.modulus_poly);
    cache.initialized = true;
  }
  return cache.modulus_ctx;
}

FieldElement BaseRingOne() {
  FieldElement one;
  NTL::set(one);
  return one;
}

ZZ_pEX LiftBaseToExtension(const FieldElement &x) {
  ZZ_pEX out;
  NTL::clear(out);
  NTL::SetCoeff(out, 0, x);
  return out;
}

void ReduceExtensionElementInPlace(ZZ_pEX &x, const ZZ_pEX &extension_modulus) {
  const long extension_degree = ExtensionDegreeOrThrow(
      extension_modulus, "ReduceExtensionElementInPlace");
  if (NTL::deg(x) >= extension_degree) {
    const NTL::ZZ_pEXModulus &mod_ctx = ExtensionModulusContextOrThrow(
        extension_modulus, "ReduceExtensionElementInPlace");
    NTL::rem(x, x, mod_ctx);
    x.normalize();
  }
}

ZZ_pEX ExtensionZero() {
  ZZ_pEX out;
  NTL::clear(out);
  return out;
}

ZZ_pEX ExtensionOne() {
  return LiftBaseToExtension(BaseRingOne());
}

ZZ_pEX MulExtensionByBaseConstant(const ZZ_pEX &a, const FieldElement &scalar) {
  ZZ_pEX out;
  if (NTL::deg(a) < 0 || scalar == 0) {
    NTL::clear(out);
    return out;
  }
  NTL::mul(out, a, scalar);
  out.normalize();
  return out;
}

ZZ_pEX SubBaseConstantFromExtension(const ZZ_pEX &a, const FieldElement &c) {
  ZZ_pEX out = a;
  NTL::SetCoeff(out, 0, NTL::coeff(out, 0) - c);
  out.normalize();
  return out;
}

ZZ_pEX AddExtension(const ZZ_pEX &a, const ZZ_pEX &b,
                    const ZZ_pEX &extension_modulus) {
  ZZ_pEX out = a + b;
  const long extension_degree =
      ExtensionDegreeOrThrow(extension_modulus, "AddExtension");
  if (NTL::deg(out) >= extension_degree) {
    const NTL::ZZ_pEXModulus &mod_ctx =
        ExtensionModulusContextOrThrow(extension_modulus, "AddExtension");
    NTL::rem(out, out, mod_ctx);
    out.normalize();
  }
  return out;
}

ZZ_pEX SubExtension(const ZZ_pEX &a, const ZZ_pEX &b,
                    const ZZ_pEX &extension_modulus) {
  ZZ_pEX out = a - b;
  const long extension_degree =
      ExtensionDegreeOrThrow(extension_modulus, "SubExtension");
  if (NTL::deg(out) >= extension_degree) {
    const NTL::ZZ_pEXModulus &mod_ctx =
        ExtensionModulusContextOrThrow(extension_modulus, "SubExtension");
    NTL::rem(out, out, mod_ctx);
    out.normalize();
  }
  return out;
}

ZZ_pEX MulExtension(const ZZ_pEX &a, const ZZ_pEX &b,
                    const ZZ_pEX &extension_modulus) {
  const long extension_degree =
      ExtensionDegreeOrThrow(extension_modulus, "MulExtension");

  const long deg_a = NTL::deg(a);
  const long deg_b = NTL::deg(b);

  ZZ_pEX out;
  if (deg_a <= 0) {
    out = MulExtensionByBaseConstant(b, NTL::coeff(a, 0));
  } else if (deg_b <= 0) {
    out = MulExtensionByBaseConstant(a, NTL::coeff(b, 0));
  } else {
    const NTL::ZZ_pEXModulus &mod_ctx =
        ExtensionModulusContextOrThrow(extension_modulus, "MulExtension");
    NTL::MulMod(out, a, b, mod_ctx);
  }

  if (NTL::deg(out) >= extension_degree) {
    const NTL::ZZ_pEXModulus &mod_ctx =
        ExtensionModulusContextOrThrow(extension_modulus, "MulExtension");
    NTL::rem(out, out, mod_ctx);
    out.normalize();
  }
  return out;
}

ZZ_pEX EqFactorExtension(const ZZ_pEX &z_i, const ZZ_pEX &x_i,
                         const ZZ_pEX &extension_modulus) {
  const ZZ_pEX one = ExtensionOne();
  const ZZ_pEX zx = MulExtension(z_i, x_i, extension_modulus);
  const ZZ_pEX two_zx = AddExtension(zx, zx, extension_modulus);
  const ZZ_pEX linear = SubExtension(SubExtension(one, z_i, extension_modulus),
                                     x_i, extension_modulus);
  return AddExtension(linear, two_zx, extension_modulus);
}

ZZ_pEX EqFactorExtensionFromBase(const FieldElement &z_i, const ZZ_pEX &x_i,
                                 const ZZ_pEX &extension_modulus) {
  const FieldElement one = BaseRingOne();
  const FieldElement factor0 = one - z_i;
  const FieldElement delta_factor = z_i - factor0;
  return AddExtension(LiftBaseToExtension(factor0),
                      MulExtensionByBaseConstant(x_i, delta_factor),
                      extension_modulus);
}

ZZ_pEX EvalExtensionQuadraticPoly(const ExtensionQuadraticPoly &p,
                                  const ZZ_pEX &x,
                                  const ZZ_pEX &extension_modulus) {
  const ZZ_pEX t = AddExtension(p.a1, MulExtension(p.a2, x, extension_modulus),
                                extension_modulus);
  return AddExtension(p.a0, MulExtension(t, x, extension_modulus),
                      extension_modulus);
}

void AppendSerializedExtensionElement(Bytes &out, const ZZ_pEX &x,
                                      const ZZ_pEX &extension_modulus,
                                      const char *func_name) {
  const long ext_degree = ExtensionDegreeOrThrow(extension_modulus, func_name);
  ZZ_pEX reduced = x;
  ReduceExtensionElementInPlace(reduced, extension_modulus);

  AppendU64(out, static_cast<std::uint64_t>(ext_degree));
  for (long i = 0; i < ext_degree; ++i) {
    AppendSerializedFieldElement(out, NTL::coeff(reduced, i), func_name);
  }
}

Bytes SerializeExtensionElement(const ZZ_pEX &x,
                                const ZZ_pEX &extension_modulus) {
  Bytes out;
  AppendSerializedExtensionElement(out, x, extension_modulus,
                                   "SerializeExtensionElement");
  return out;
}

void AbsorbExtensionElement(HashTranscript &transcript, const ZZ_pEX &x,
                            const ZZ_pEX &extension_modulus) {
  transcript.AbsorbBytes(SerializeExtensionElement(x, extension_modulus));
}

void AbsorbExtensionQuadraticPoly(HashTranscript &transcript,
                                  const ExtensionQuadraticPoly &p,
                                  const ZZ_pEX &extension_modulus) {
  AbsorbExtensionElement(transcript, p.a0, extension_modulus);
  AbsorbExtensionElement(transcript, p.a1, extension_modulus);
  AbsorbExtensionElement(transcript, p.a2, extension_modulus);
}

ZZ_pEX SampleExtensionChallenge(const HashTranscript &transcript,
                                const std::string &label,
                                const ZZ_pEX &extension_modulus) {
  const long ext_degree =
      ExtensionDegreeOrThrow(extension_modulus, "SampleExtensionChallenge");

  ZZ_pEX sampled;
  NTL::clear(sampled);
  for (long i = 0; i < ext_degree; ++i) {
    const FieldElement coeff = transcript.ChallengeFieldElement(
        "ext/" + label + "/coeff/" + std::to_string(i));
    NTL::SetCoeff(sampled, i, coeff);
  }
  ReduceExtensionElementInPlace(sampled, extension_modulus);
  return sampled;
}

ZZ_pEX EvalLineAtExtensionWithInvDenom(const ZZ_pEX &x, const FieldElement &x1,
                                       const ZZ_pEX &y1, const ZZ_pEX &y2,
                                       const FieldElement &inv_denom,
                                       const ZZ_pEX &extension_modulus) {
  const ZZ_pEX delta_y = SubExtension(y2, y1, extension_modulus);
  const ZZ_pEX slope = MulExtensionByBaseConstant(delta_y, inv_denom);
  const ZZ_pEX delta_x = SubBaseConstantFromExtension(x, x1);
  const ZZ_pEX correction = MulExtension(slope, delta_x, extension_modulus);
  return AddExtension(y1, correction, extension_modulus);
}

ZZ_pEX EvalLineAtExtension(const ZZ_pEX &x, const FieldElement &x1,
                           const ZZ_pEX &y1, const FieldElement &x2,
                           const ZZ_pEX &y2, const ZZ_pEX &extension_modulus) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->ext_eval_line_at_ns : nullptr,
                    prof ? &prof->ext_eval_line_at_calls : nullptr);

  const FieldElement denom = x2 - x1;
  FieldElement inv_denom;
  if (!basefold_pcs_internal::TryInvertBaseUnit(inv_denom, denom)) {
    LogicError("EvalLineAtExtension: denominator is not invertible");
  }
  return EvalLineAtExtensionWithInvDenom(x, x1, y1, y2, inv_denom,
                                         extension_modulus);
}

std::vector<ZZ_pEX> LiftOracleToExtension(const Oracle &oracle) {
  std::vector<ZZ_pEX> out(static_cast<std::size_t>(oracle.length()));
  for (long i = 0; i < oracle.length(); ++i) {
    out[static_cast<std::size_t>(i)] =
        LiftBaseToExtension(oracle[static_cast<std::size_t>(i)]);
  }
  return out;
}

std::vector<ZZ_pEX> BooleanEvalTableFromMonomialCoeffsExtension(
    const std::vector<ZZ_pEX> &coeffs, long k,
    const ZZ_pEX &extension_modulus) {
  if (k < 0) {
    LogicError(
        "BooleanEvalTableFromMonomialCoeffsExtension: negative dimension");
  }
  if (static_cast<long>(coeffs.size()) != (1L << k)) {
    LogicError("BooleanEvalTableFromMonomialCoeffsExtension: length mismatch");
  }

  std::vector<ZZ_pEX> eval = coeffs;
  for (long bit = 0; bit < k; ++bit) {
    const long step = 1L << bit;
    for (long mask = 0; mask < static_cast<long>(eval.size()); ++mask) {
      if (mask & step) {
        eval[static_cast<std::size_t>(mask)] = AddExtension(
            eval[static_cast<std::size_t>(mask)],
            eval[static_cast<std::size_t>(mask ^ step)], extension_modulus);
      }
    }
  }
  return eval;
}

std::vector<ZZ_pEX> LiftFieldVecToExtension(const FieldVec &values) {
  std::vector<ZZ_pEX> out(static_cast<std::size_t>(values.length()));
  for (long i = 0; i < values.length(); ++i) {
    out[static_cast<std::size_t>(i)] =
        LiftBaseToExtension(values[static_cast<std::size_t>(i)]);
  }
  return out;
}

std::vector<ZZ_pEX> Msg0CoeffsAtSuffixChallenges(
    const vec_ZZ_pE &f_coeffs, long kappa, const std::vector<ZZ_pEX> &r_by_level,
    const ZZ_pEX &extension_modulus) {
  if (kappa < 0) {
    LogicError("Msg0CoeffsAtSuffixChallenges: negative kappa");
  }
  const long d = static_cast<long>(r_by_level.size());
  const long point_dim = kappa + d;
  if (f_coeffs.length() != (1L << point_dim)) {
    LogicError("Msg0CoeffsAtSuffixChallenges: f_coeffs length mismatch");
  }

  std::vector<ZZ_pEX> cur(static_cast<std::size_t>(f_coeffs.length()));
  for (long i = 0; i < f_coeffs.length(); ++i) {
    cur[static_cast<std::size_t>(i)] =
        LiftBaseToExtension(f_coeffs[static_cast<std::size_t>(i)]);
  }

  long cur_len = static_cast<long>(cur.size());
  for (long var = point_dim; var-- > kappa;) {
    const long half = cur_len / 2;
    const ZZ_pEX &r_var = r_by_level[static_cast<std::size_t>(var - kappa)];
    for (long i = 0; i < half; ++i) {
      cur[static_cast<std::size_t>(i)] = AddExtension(
          cur[static_cast<std::size_t>(i)],
          MulExtension(cur[static_cast<std::size_t>(i + half)], r_var,
                       extension_modulus),
          extension_modulus);
    }
    cur.resize(static_cast<std::size_t>(half));
    cur_len = half;
  }
  return cur;
}

std::vector<ZZ_pEX> EncodeC0Extension(const std::vector<ZZ_pEX> &msg0_coeffs,
                                      const FoldableCodeParams &params,
                                      const ZZ_pEX &extension_modulus) {
  if (static_cast<long>(msg0_coeffs.size()) != params.k0) {
    LogicError("EncodeC0Extension: msg0_coeffs has wrong length");
  }

  const long n0 = basefold_pcs_internal::CodewordLengthAtLevelNoValidate(params, 0);
  std::vector<ZZ_pEX> out(static_cast<std::size_t>(n0), ExtensionZero());
  for (long j = 0; j < n0; ++j) {
    ZZ_pEX acc = ExtensionZero();
    for (long row = 0; row < params.k0; ++row) {
      const ZZ_pEX g =
          LiftBaseToExtension(params.G0[static_cast<std::size_t>(row)][j]);
      acc = AddExtension(
          acc,
          MulExtension(msg0_coeffs[static_cast<std::size_t>(row)], g,
                       extension_modulus),
          extension_modulus);
    }
    out[static_cast<std::size_t>(j)] = acc;
  }
  return out;
}

class ExtensionSumcheckProver {
 public:
  ExtensionSumcheckProver(const FieldVec &f_coeffs,
                          const std::vector<FieldElement> &z,
                          const ZZ_pEX &extension_modulus)
      : extension_modulus_(extension_modulus) {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_sumcheck_init_ns : nullptr,
                      prof ? &prof->ext_sumcheck_init_calls : nullptr);

    const long n = f_coeffs.length();
    if (!basefold_pcs_internal::IsPowerOfTwoLong(n)) {
      LogicError("ExtensionSumcheckProver: f_coeffs length must be 2^d");
    }
    d_ = basefold_pcs_internal::Log2ExactPowerOfTwoLong(n);
    InitializePointAndPrefixProducts(z);

    const std::vector<ZZ_pEX> lifted_coeffs = LiftFieldVecToExtension(f_coeffs);
    f_eval_table_ = BooleanEvalTableFromMonomialCoeffsExtension(
        lifted_coeffs, d_, extension_modulus_);
    suffix_eq_prod_ = ExtensionOne();
  }

  ExtensionSumcheckProver(
      const SumcheckMonomialPrecomputation &precomputation,
      const std::vector<FieldElement> &z, const ZZ_pEX &extension_modulus)
      : extension_modulus_(extension_modulus) {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_sumcheck_init_ns : nullptr,
                      prof ? &prof->ext_sumcheck_init_calls : nullptr);

    if (!precomputation.valid) {
      LogicError("ExtensionSumcheckProver: precomputation must be valid");
    }

    const long n = precomputation.f_eval_table.length();
    if (!basefold_pcs_internal::IsPowerOfTwoLong(n)) {
      LogicError(
          "ExtensionSumcheckProver: precomputation f_eval_table length must be 2^d");
    }
    d_ = basefold_pcs_internal::Log2ExactPowerOfTwoLong(n);
    if (precomputation.d != d_) {
      LogicError("ExtensionSumcheckProver: precomputation dimension mismatch");
    }

    InitializePointAndPrefixProducts(z);
    f_eval_table_ = LiftFieldVecToExtension(precomputation.f_eval_table);
    suffix_eq_prod_ = ExtensionOne();
  }

 private:
  void InitializePointAndPrefixProducts(const std::vector<FieldElement> &z) {
    if (static_cast<long>(z.size()) != d_) {
      LogicError("ExtensionSumcheckProver: z dimension mismatch");
    }

    cur_k_ = d_;
    z_ = z;

    prefix_eq_by_vars_.resize(static_cast<std::size_t>(d_));
    if (d_ > 0) {
      const FieldElement one = BaseRingOne();
      prefix_eq_by_vars_[0].SetLength(1);
      prefix_eq_by_vars_[0][0] = one;

      for (long t = 1; t < d_; ++t) {
        const FieldElement z_var = z_[static_cast<std::size_t>(t - 1)];
        const FieldElement factor0 = one - z_var;
        const FieldElement factor1 = z_var;

        const FieldVec &prev = prefix_eq_by_vars_[static_cast<std::size_t>(t - 1)];
        const long old = prev.length();
        prefix_eq_by_vars_[static_cast<std::size_t>(t)].SetLength(2 * old);
        FieldVec &cur = prefix_eq_by_vars_[static_cast<std::size_t>(t)];
        for (long mask = 0; mask < old; ++mask) {
          const FieldElement &base = prev[mask];
          cur[mask] = base * factor0;
          cur[mask + old] = base * factor1;
        }
      }
    }
  }

 public:
  ExtensionQuadraticPoly CurrentPolynomial() const {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_sumcheck_current_poly_ns : nullptr,
                      prof ? &prof->ext_sumcheck_current_poly_calls : nullptr);

    if (cur_k_ <= 0) {
      LogicError("ExtensionSumcheckProver::CurrentPolynomial: no variables");
    }

    const long k = cur_k_;
    const long n = static_cast<long>(f_eval_table_.size());
    if (n != (1L << k)) {
      LogicError(
          "ExtensionSumcheckProver::CurrentPolynomial: internal length mismatch");
    }

    const long half = 1L << (k - 1);
    const FieldVec &prefix =
        prefix_eq_by_vars_[static_cast<std::size_t>(k - 1)];
    if (prefix.length() != half) {
      LogicError(
          "ExtensionSumcheckProver::CurrentPolynomial: prefix size mismatch");
    }

    const FieldElement one_base = BaseRingOne();
    const FieldElement z_k_base = z_[static_cast<std::size_t>(k - 1)];
    const FieldElement factor0_base = one_base - z_k_base;
    const FieldElement delta_factor_base = z_k_base - factor0_base;

    ExtensionQuadraticPoly out;
    out.a0 = ExtensionZero();
    out.a1 = ExtensionZero();
    out.a2 = ExtensionZero();

    for (long mask = 0; mask < half; ++mask) {
      const ZZ_pEX common =
          MulExtensionByBaseConstant(suffix_eq_prod_, prefix[mask]);

      const ZZ_pEX eq0 = MulExtensionByBaseConstant(common, factor0_base);
      const ZZ_pEX delta_eq =
          MulExtensionByBaseConstant(common, delta_factor_base);

      const ZZ_pEX &f0 = f_eval_table_[static_cast<std::size_t>(mask)];
      const ZZ_pEX &f1 = f_eval_table_[static_cast<std::size_t>(mask + half)];
      const ZZ_pEX delta_f = SubExtension(f1, f0, extension_modulus_);

      out.a0 = AddExtension(out.a0, MulExtension(f0, eq0, extension_modulus_),
                            extension_modulus_);
      const ZZ_pEX term1 = MulExtension(f0, delta_eq, extension_modulus_);
      const ZZ_pEX term2 = MulExtension(delta_f, eq0, extension_modulus_);
      out.a1 =
          AddExtension(out.a1, AddExtension(term1, term2, extension_modulus_),
                       extension_modulus_);
      out.a2 = AddExtension(out.a2,
                            MulExtension(delta_f, delta_eq, extension_modulus_),
                            extension_modulus_);
    }

    return out;
  }

  void ReceiveChallenge(const ZZ_pEX &r_kminus1) {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_sumcheck_receive_challenge_ns : nullptr,
                      prof ? &prof->ext_sumcheck_receive_challenge_calls
                           : nullptr);

    if (cur_k_ <= 0) {
      LogicError("ExtensionSumcheckProver::ReceiveChallenge: no variables");
    }

    const long k = cur_k_;
    const long n = static_cast<long>(f_eval_table_.size());
    if (n != (1L << k)) {
      LogicError(
          "ExtensionSumcheckProver::ReceiveChallenge: internal length mismatch");
    }

    const ZZ_pEX eq = EqFactorExtensionFromBase(
        z_[static_cast<std::size_t>(k - 1)], r_kminus1, extension_modulus_);
    suffix_eq_prod_ = MulExtension(suffix_eq_prod_, eq, extension_modulus_);

    const long half = n / 2;
    for (long i = 0; i < half; ++i) {
      const ZZ_pEX &f0 = f_eval_table_[static_cast<std::size_t>(i)];
      const ZZ_pEX &f1 = f_eval_table_[static_cast<std::size_t>(i + half)];
      const ZZ_pEX delta_f = SubExtension(f1, f0, extension_modulus_);
      f_eval_table_[static_cast<std::size_t>(i)] =
          AddExtension(f0, MulExtension(delta_f, r_kminus1, extension_modulus_),
                       extension_modulus_);
    }
    f_eval_table_.resize(static_cast<std::size_t>(half));
    --cur_k_;
  }

 private:
  long d_ = 0;
  long cur_k_ = 0;
  ZZ_pEX extension_modulus_;
  std::vector<FieldElement> z_;
  std::vector<ZZ_pEX> f_eval_table_;
  std::vector<FieldVec> prefix_eq_by_vars_;
  ZZ_pEX suffix_eq_prod_;
};

bool CheckExtensionSumcheckRelations(
    const std::vector<ExtensionQuadraticPoly> &h_by_level,
    const std::vector<ZZ_pEX> &r, const FieldElement &claimed_y,
    const ZZ_pEX &extension_modulus) {
  const long d = static_cast<long>(h_by_level.size());
  if (static_cast<long>(r.size()) != d) {
    return false;
  }
  if (d == 0) {
    return true;
  }

  const ZZ_pEX zero = ExtensionZero();
  const ZZ_pEX one = ExtensionOne();
  const ZZ_pEX claimed_y_ext = LiftBaseToExtension(claimed_y);

  const ExtensionQuadraticPoly &h_d =
      h_by_level[static_cast<std::size_t>(d - 1)];
  const ZZ_pEX hd0 = EvalExtensionQuadraticPoly(h_d, zero, extension_modulus);
  const ZZ_pEX hd1 = EvalExtensionQuadraticPoly(h_d, one, extension_modulus);
  if (AddExtension(hd0, hd1, extension_modulus) != claimed_y_ext) {
    return false;
  }

  for (long k = 1; k < d; ++k) {
    const ExtensionQuadraticPoly &h_k =
        h_by_level[static_cast<std::size_t>(k - 1)];
    const ExtensionQuadraticPoly &h_kp1 =
        h_by_level[static_cast<std::size_t>(k)];
    const ZZ_pEX lhs =
        AddExtension(EvalExtensionQuadraticPoly(h_k, zero, extension_modulus),
                     EvalExtensionQuadraticPoly(h_k, one, extension_modulus),
                     extension_modulus);
    const ZZ_pEX rhs =
        EvalExtensionQuadraticPoly(h_kp1, r[static_cast<std::size_t>(k)],
                                   extension_modulus);
    if (lhs != rhs) {
      return false;
    }
  }

  return true;
}

const ExtensionCommitRoundLevelPrecomputation *
LookupExtensionCommitRoundLevelPrecomputation(
    const ExtensionCommitRoundPrecomputation *precomputation, long level_i,
    long n_i) {
  if (precomputation == nullptr || level_i < 0 ||
      level_i >= static_cast<long>(precomputation->levels.size())) {
    return nullptr;
  }
  const ExtensionCommitRoundLevelPrecomputation &level_cache =
      precomputation->levels[static_cast<std::size_t>(level_i)];
  const long inv_len = level_cache.inv_denoms.length();
  if (inv_len == 1 || inv_len == n_i) {
    return &level_cache;
  }
  return nullptr;
}

ExtensionSumcheckProver BuildExtensionSumcheckProver(
    const FieldVec &f_coeffs, const std::vector<FieldElement> &z,
    const ZZ_pEX &extension_modulus,
    const BaseFoldPCSCommitArtifacts &commit_artifacts) {
  if (commit_artifacts.base_sumcheck_precomputation.valid) {
    return ExtensionSumcheckProver(commit_artifacts.base_sumcheck_precomputation,
                                   z, extension_modulus);
  }
  return ExtensionSumcheckProver(f_coeffs, z, extension_modulus);
}

void ProverCommitRoundExtensionNoValidate(std::vector<ZZ_pEX> &pi_i,
                                          const std::vector<ZZ_pEX> &pi_ip1,
                                          const ZZ_pEX &alpha_i, long level_i,
                                          const FoldableCodeParams &params,
                                          const ZZ_pEX &extension_modulus,
                                          const ExtensionCommitRoundPrecomputation
                                              *precomputation = nullptr) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->ext_prover_commit_round_ns : nullptr,
                    prof ? &prof->ext_prover_commit_round_calls : nullptr);

  const long n_i =
      basefold_pcs_internal::CodewordLengthAtLevelNoValidate(params, level_i);
  if (static_cast<long>(pi_ip1.size()) != 2 * n_i) {
    LogicError("ProverCommitRoundExtensionNoValidate: pi_ip1 has wrong length");
  }
  const Oracle &diag = params.diag_T[static_cast<std::size_t>(level_i)];
  if (diag.length() != n_i) {
    LogicError("ProverCommitRoundExtensionNoValidate: diag_T length mismatch");
  }
  pi_i.resize(static_cast<std::size_t>(n_i));
  if (n_i == 0) {
    return;
  }

  const FieldElement one = BaseRingOne();
  const FieldElement zeta_minus_one = params.zeta - one;
  if (zeta_minus_one == 0) {
    LogicError("ProverCommitRoundExtensionNoValidate: zeta must not be 1");
  }
  constexpr long kParallelThreshold = 4096;

  const ExtensionCommitRoundLevelPrecomputation *level_cache =
      LookupExtensionCommitRoundLevelPrecomputation(precomputation, level_i, n_i);
  if (level_cache != nullptr) {
    const long inv_len = level_cache->inv_denoms.length();
    if (inv_len == 1) {
      const FieldElement &x1 = diag[0];
      const FieldElement &inv_denom = level_cache->inv_denoms[0];
      basefold_pcs_internal::ForEachIndexMaybeParallel(
          0, n_i, kParallelThreshold, [&](long j) {
            pi_i[static_cast<std::size_t>(j)] = EvalLineAtExtensionWithInvDenom(
                alpha_i, x1, pi_ip1[static_cast<std::size_t>(j)],
                pi_ip1[static_cast<std::size_t>(j + n_i)], inv_denom,
                extension_modulus);
          });
      return;
    }
    if (inv_len == n_i) {
      basefold_pcs_internal::ForEachIndexMaybeParallel(
          0, n_i, kParallelThreshold, [&](long j) {
            const FieldElement &x1 = diag[static_cast<std::size_t>(j)];
            pi_i[static_cast<std::size_t>(j)] = EvalLineAtExtensionWithInvDenom(
                alpha_i, x1, pi_ip1[static_cast<std::size_t>(j)],
                pi_ip1[static_cast<std::size_t>(j + n_i)],
                level_cache->inv_denoms[static_cast<std::size_t>(j)],
                extension_modulus);
          });
      return;
    }
  }

  const FieldElement first_t = diag[0];
  bool all_equal = true;
  for (long j = 1; j < n_i; ++j) {
    if (diag[static_cast<std::size_t>(j)] != first_t) {
      all_equal = false;
      break;
    }
  }
  if (all_equal) {
    FieldElement inv_denom;
    if (!basefold_pcs_internal::TryInvertBaseUnit(inv_denom,
                                                  zeta_minus_one * first_t)) {
      LogicError(
          "ProverCommitRoundExtensionNoValidate: denominator is not invertible");
    }
    basefold_pcs_internal::ForEachIndexMaybeParallel(
        0, n_i, kParallelThreshold, [&](long j) {
          pi_i[static_cast<std::size_t>(j)] = EvalLineAtExtensionWithInvDenom(
              alpha_i, first_t, pi_ip1[static_cast<std::size_t>(j)],
              pi_ip1[static_cast<std::size_t>(j + n_i)], inv_denom,
              extension_modulus);
        });
    return;
  }

  std::vector<FieldElement> denoms(static_cast<std::size_t>(n_i));
  basefold_pcs_internal::ForEachIndexMaybeParallel(
      0, n_i, kParallelThreshold, [&](long j) {
        denoms[static_cast<std::size_t>(j)] =
            zeta_minus_one * diag[static_cast<std::size_t>(j)];
      });

  std::vector<FieldElement> inv_denoms;
  if (!basefold_pcs_internal::BatchInvertBaseUnits(inv_denoms, denoms)) {
    inv_denoms.resize(static_cast<std::size_t>(n_i));
    basefold_pcs_internal::ForEachIndexMaybeParallel(
        0, n_i, kParallelThreshold, [&](long j) {
          if (!basefold_pcs_internal::TryInvertBaseUnit(
                  inv_denoms[static_cast<std::size_t>(j)],
                  denoms[static_cast<std::size_t>(j)])) {
            LogicError(
                "ProverCommitRoundExtensionNoValidate: denominator is not invertible");
          }
        });
  }

  basefold_pcs_internal::ForEachIndexMaybeParallel(
      0, n_i, kParallelThreshold, [&](long j) {
        const FieldElement &x1 = diag[static_cast<std::size_t>(j)];
        pi_i[static_cast<std::size_t>(j)] = EvalLineAtExtensionWithInvDenom(
            alpha_i, x1, pi_ip1[static_cast<std::size_t>(j)],
            pi_ip1[static_cast<std::size_t>(j + n_i)],
            inv_denoms[static_cast<std::size_t>(j)], extension_modulus);
      });
}

std::size_t MaxSerializedFieldElementSizeOrThrow(const char *func_name) {
  const long r = ZZ_pE::degree();
  if (r <= 0) {
    const std::string msg =
        std::string(func_name) + ": invalid extension degree";
    LogicError(msg.c_str());
  }
  const long coeff_max_bytes = std::max<long>(1, NumBytes(ZZ_p::modulus()));
  return static_cast<std::size_t>(8) +
         static_cast<std::size_t>(r) *
             (static_cast<std::size_t>(8) +
              static_cast<std::size_t>(coeff_max_bytes));
}

std::size_t MaxSerializedExtensionElementSizeOrThrow(
    const ZZ_pEX &extension_modulus, const char *func_name) {
  const long ext_degree = ExtensionDegreeOrThrow(extension_modulus, func_name);
  return static_cast<std::size_t>(8) +
         static_cast<std::size_t>(ext_degree) *
             MaxSerializedFieldElementSizeOrThrow(func_name);
}

Digest HashExtensionLeaf(long index, const ZZ_pEX &value,
                         const ZZ_pEX &extension_modulus,
                         Bytes *scratch_payload = nullptr) {
  Bytes local_payload;
  Bytes &payload =
      (scratch_payload != nullptr) ? *scratch_payload : local_payload;
  payload.clear();
  payload.push_back(static_cast<Byte>(0x30));
  AppendU64(payload, static_cast<std::uint64_t>(index));

  const std::size_t encoded_len_pos = payload.size();
  AppendU64(payload, 0);
  const std::size_t encoded_start = payload.size();
  AppendSerializedExtensionElement(payload, value, extension_modulus,
                                   "HashExtensionLeaf");
  const std::uint64_t encoded_len =
      static_cast<std::uint64_t>(payload.size() - encoded_start);
  StoreU64BigEndian(payload.data() + encoded_len_pos, encoded_len);

  return HashDigest(payload.data(), payload.size(), "HashExtensionLeaf");
}

Digest HashExtensionNode(const Digest &left, const Digest &right) {
  std::array<Byte, 1 + 32 + 32> payload{};
  payload[0] = static_cast<Byte>(0x31);
  std::copy(left.begin(), left.end(), payload.begin() + 1);
  std::copy(right.begin(), right.end(), payload.begin() + 1 + left.size());
  return HashDigest(payload.data(), payload.size(), "HashExtensionNode");
}

Digest HashExtensionRawRootEmpty() {
  const Byte payload[1] = {static_cast<Byte>(0x32)};
  return HashDigest(payload, sizeof(payload), "HashExtensionRawRootEmpty");
}

Digest HashExtensionRootWithCount(long leaf_count, const Digest &raw_root) {
  if (leaf_count <= 0) {
    LogicError("HashExtensionRootWithCount: invalid leaf_count");
  }
  std::array<Byte, 1 + 8 + 32> payload{};
  payload[0] = static_cast<Byte>(0x33);
  StoreU64BigEndian(payload.data() + 1, static_cast<std::uint64_t>(leaf_count));
  std::copy(raw_root.begin(), raw_root.end(), payload.begin() + 1 + 8);
  return HashDigest(payload.data(), payload.size(),
                    "HashExtensionRootWithCount");
}

class ExtensionMerkleTree {
 public:
  static ExtensionMerkleTree Build(const std::vector<ZZ_pEX> &oracle,
                                   const ZZ_pEX &extension_modulus) {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_merkle_tree_build_ns : nullptr,
                      prof ? &prof->ext_merkle_tree_build_calls : nullptr);

    if (oracle.empty()) {
      LogicError("ExtensionMerkleTree::Build: oracle must be non-empty");
    }

    ExtensionMerkleTree tree;
    tree.leaf_count_ = static_cast<long>(oracle.size());

    Bytes leaf_payload;
    leaf_payload.reserve(static_cast<std::size_t>(1 + 8 + 8) +
                         MaxSerializedExtensionElementSizeOrThrow(
                             extension_modulus, "ExtensionMerkleTree::Build"));

    std::vector<Digest> level(oracle.size());
    for (long i = 0; i < tree.leaf_count_; ++i) {
      level[static_cast<std::size_t>(i)] =
          HashExtensionLeaf(i, oracle[static_cast<std::size_t>(i)],
                            extension_modulus, &leaf_payload);
    }

    tree.levels_.push_back(level);
    while (level.size() > 1) {
      if ((level.size() % 2U) == 1U) {
        level.push_back(level.back());
      }
      tree.levels_.back() = level;

      std::vector<Digest> next;
      next.reserve(level.size() / 2U);
      for (std::size_t i = 0; i < level.size(); i += 2U) {
        next.push_back(HashExtensionNode(level[i], level[i + 1U]));
      }
      tree.levels_.push_back(next);
      level = std::move(next);
    }

    tree.raw_root_ = tree.levels_.empty() ? HashExtensionRawRootEmpty()
                                          : tree.levels_.back()[0];
    return tree;
  }

  MerkleRoot Root() const {
    return HashExtensionRootWithCount(leaf_count_, raw_root_);
  }

  ExtensionMerkleMultiproof OpenMany(
      const std::vector<ZZ_pEX> &oracle,
      const std::vector<long> &queried_indices) const {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->ext_merkle_tree_open_ns : nullptr,
                      prof ? &prof->ext_merkle_tree_open_calls : nullptr);

    if (static_cast<long>(oracle.size()) != leaf_count_) {
      LogicError("ExtensionMerkleTree::OpenMany: oracle length mismatch");
    }
    const std::vector<long> unique =
        multiproof_planner::SortAndValidateIndicesOrThrow(
            leaf_count_, queried_indices, "ExtensionMerkleTree::OpenMany");

    ExtensionMerkleMultiproof proof;
    proof.queried_indices = unique;
    proof.values.resize(unique.size());
    if (unique.empty()) {
      return proof;
    }

    const LocalMultiproofPlan plan =
        multiproof_planner::BuildPlanFromSortedUnique(leaf_count_, unique);
    for (std::size_t i = 0; i < unique.size(); ++i) {
      proof.values[i] = oracle[static_cast<std::size_t>(unique[i])];
    }
    proof.sibling_hashes = multiproof_replay::CollectSiblingHashesForPlan(
        plan, [&](std::size_t level_index, long sibling) {
          return levels_[level_index][static_cast<std::size_t>(sibling)];
        });
    return proof;
  }

 private:
  long leaf_count_ = 0;
  Digest raw_root_{};
  std::vector<std::vector<Digest>> levels_;
};

bool VerifyExtensionMerkleMultiproofNoProfile(
    const MerkleRoot &root, long leaf_count,
    const std::vector<long> &queried_indices,
    const ExtensionMerkleMultiproof &proof, const ZZ_pEX &extension_modulus) {
  if (leaf_count <= 0) {
    return false;
  }
  if (proof.values.size() != queried_indices.size()) {
    return false;
  }
  if (!multiproof_planner::IsSortedUniqueIndicesInRange(leaf_count,
                                                        queried_indices)) {
    return false;
  }

  const LocalMultiproofPlan plan =
      multiproof_planner::BuildPlanFromSortedUnique(leaf_count, queried_indices);
  if (proof.sibling_hashes.size() !=
      static_cast<std::size_t>(plan.stats.unique_sibling_count)) {
    return false;
  }
  if (queried_indices.empty()) {
    return proof.sibling_hashes.empty();
  }

  std::vector<std::pair<long, Digest>> current(queried_indices.size());
  for (std::size_t i = 0; i < queried_indices.size(); ++i) {
    current[i] = {queried_indices[i],
                  HashExtensionLeaf(queried_indices[i], proof.values[i],
                                    extension_modulus)};
  }

  Digest raw_root;
  if (!multiproof_replay::ReplayPlanToRawRoot(
          plan, proof.sibling_hashes, std::move(current),
          [](const Digest &lhs, const Digest &rhs) {
            return HashExtensionNode(lhs, rhs);
          },
          &raw_root)) {
    return false;
  }
  return HashExtensionRootWithCount(leaf_count, raw_root) == root;
}

bool VerifyExtensionMerkleMultiproof(const MerkleRoot &root, long leaf_count,
                                     const std::vector<long> &queried_indices,
                                     const ExtensionMerkleMultiproof &proof,
                                     const ZZ_pEX &extension_modulus) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->ext_merkle_verify_opening_ns : nullptr,
                    prof ? &prof->ext_merkle_verify_opening_calls : nullptr);
  return VerifyExtensionMerkleMultiproofNoProfile(
      root, leaf_count, queried_indices, proof, extension_modulus);
}

BaseFoldPCSEvalProof
ProveEvalWithExtensionChallengesFromCommittedTopOracleUnchecked(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  basefold_pcs_internal::ValidateCommittedTopOracleArtifactsOrThrow(
      commit_artifacts,
      basefold_pcs_internal::CodewordLengthAtLevelNoValidate(params, params.d),
      "BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked");

  if (params.d == 0) {
    return BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
        f_coeffs, z, claimed_y, num_queries, params, commit_artifacts);
  }

  const ZZ_pEX &extension_modulus = challenge_cfg.challenge_extension_modulus;
  ExtensionDegreeOrThrow(
      extension_modulus,
      "BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked");

  BaseFoldPCSEvalProof proof;
  proof.extension.has_extension_payload = true;
  proof.commitments.roots_by_level.clear();
  proof.extension.roots_by_level.resize(static_cast<std::size_t>(params.d));
  proof.commitments.roots_by_level.push_back(commit_artifacts.root_d);

  HashTranscript transcript = basefold_pcs_internal::MakeBaseFoldTranscript();
  basefold_pcs_internal::AbsorbPublicInput(transcript, commit_artifacts.root_d, z,
                                           claimed_y);
  AbsorbChallengeConfig(transcript, challenge_cfg);

  proof.extension.h_by_level.resize(static_cast<std::size_t>(params.d));
  std::vector<ZZ_pEX> r_by_level(static_cast<std::size_t>(params.d));
  const ExtensionCommitRoundPrecomputation *commit_round_precomputation =
      commit_artifacts.extension_commit_precomputation.levels.empty()
          ? nullptr
          : &commit_artifacts.extension_commit_precomputation;

  std::vector<std::vector<ZZ_pEX>> ext_oracles(static_cast<std::size_t>(params.d + 1));
  ext_oracles[static_cast<std::size_t>(params.d)] =
      LiftOracleToExtension(commit_artifacts.pi_d);

  std::vector<ExtensionMerkleTree> ext_merkle_trees(
      static_cast<std::size_t>(params.d));

  ExtensionSumcheckProver sumcheck_ext =
      BuildExtensionSumcheckProver(f_coeffs, z, extension_modulus,
                                   commit_artifacts);
  const ExtensionQuadraticPoly h_d_ext = sumcheck_ext.CurrentPolynomial();
  proof.extension.h_by_level[static_cast<std::size_t>(params.d - 1)] = h_d_ext;
  AbsorbExtensionQuadraticPoly(transcript, h_d_ext, extension_modulus);

  for (long i = params.d; i-- > 0;) {
    const ZZ_pEX r_i_ext =
        SampleExtensionChallenge(transcript, "r/" + std::to_string(i),
                                 extension_modulus);
    r_by_level[static_cast<std::size_t>(i)] = r_i_ext;

    ProverCommitRoundExtensionNoValidate(
        ext_oracles[static_cast<std::size_t>(i)],
        ext_oracles[static_cast<std::size_t>(i + 1)], r_i_ext, i, params,
        extension_modulus, commit_round_precomputation);
    ext_merkle_trees[static_cast<std::size_t>(i)] = ExtensionMerkleTree::Build(
        ext_oracles[static_cast<std::size_t>(i)], extension_modulus);
    proof.extension.roots_by_level[static_cast<std::size_t>(i)] =
        ext_merkle_trees[static_cast<std::size_t>(i)].Root();
    transcript.AbsorbDigest(
        proof.extension.roots_by_level[static_cast<std::size_t>(i)]);

    sumcheck_ext.ReceiveChallenge(r_i_ext);
    if (i > 0) {
      const ExtensionQuadraticPoly h_i_ext = sumcheck_ext.CurrentPolynomial();
      proof.extension.h_by_level[static_cast<std::size_t>(i - 1)] = h_i_ext;
      AbsorbExtensionQuadraticPoly(transcript, h_i_ext, extension_modulus);
    }
  }

  proof.extension.pi0_codeword = ext_oracles[0];

  const long kappa = basefold_pcs_internal::Log2ExactPowerOfTwoLong(params.k0);
  proof.extension.msg0_coeffs = Msg0CoeffsAtSuffixChallenges(
      f_coeffs, kappa, r_by_level, extension_modulus);
  const std::vector<ZZ_pEX> expected_pi0 =
      EncodeC0Extension(proof.extension.msg0_coeffs, params, extension_modulus);
  if (expected_pi0 != proof.extension.pi0_codeword) {
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked: internal pi0 mismatch");
  }

  proof.extension.base_top_query_multiproof = MerkleMultiproof{};
  proof.extension.query_multiproofs.resize(static_cast<std::size_t>(params.d));

  const long n_last =
      basefold_pcs_internal::CodewordLengthAtLevelNoValidate(params, params.d - 1);
  std::vector<IOPPQueryPlan> query_plans(static_cast<std::size_t>(num_queries));
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    query_plans[static_cast<std::size_t>(q)] =
        basefold_pcs_internal::MakeQueryPlanNoValidate(mu, params);
  }

  if (num_queries > 0) {
    const std::vector<std::vector<long>> requested_indices_by_tree =
        basefold_pcs_internal::CollectBaseQueryIndicesByTree(query_plans, params);
    proof.extension.base_top_query_multiproof = commit_artifacts.merkle_d.OpenMany(
        commit_artifacts.pi_d,
        requested_indices_by_tree[static_cast<std::size_t>(params.d)]);
    for (long tree_level = 0; tree_level < params.d; ++tree_level) {
      proof.extension.query_multiproofs[static_cast<std::size_t>(tree_level)] =
          ext_merkle_trees[static_cast<std::size_t>(tree_level)].OpenMany(
              ext_oracles[static_cast<std::size_t>(tree_level)],
              requested_indices_by_tree[static_cast<std::size_t>(tree_level)]);
    }
  }

  return proof;
}

bool VerifyEvalWithExtensionChallenges(
    const MerkleRoot &commitment_C, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const BaseFoldPCSEvalProof &proof, const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->pcs_verify_ns : nullptr,
                    prof ? &prof->pcs_verify_calls : nullptr);

  if (params.d == 0) {
    return BaseFoldPCSVerifyEval(commitment_C, z, claimed_y, num_queries, proof,
                                 params);
  }

  const ZZ_pEX &extension_modulus = challenge_cfg.challenge_extension_modulus;
  ExtensionDegreeOrThrow(extension_modulus,
                         "VerifyEvalWithExtensionChallenges");

  basefold_pcs_internal::ValidateParamsOrThrow(params);
  if (!basefold_pcs_internal::IsPowerOfTwoLong(params.k0)) {
    return false;
  }
  const long kappa = basefold_pcs_internal::Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim || num_queries < 0 ||
      !proof.extension.has_extension_payload) {
    return false;
  }
  if (static_cast<long>(proof.extension.roots_by_level.size()) != params.d ||
      static_cast<long>(proof.extension.h_by_level.size()) != params.d ||
      static_cast<long>(proof.extension.msg0_coeffs.size()) != params.k0) {
    return false;
  }
  if (!proof.query_multiproofs.empty()) {
    return false;
  }
  const long n0 = CodewordLengthAtLevel(params, 0);
  if (static_cast<long>(proof.extension.pi0_codeword.size()) != n0) {
    return false;
  }
  if (num_queries > 0 &&
      !basefold_pcs_internal::HasMerkleMultiproofPayload(
          proof.extension.base_top_query_multiproof)) {
    return false;
  }
  if (static_cast<long>(proof.extension.query_multiproofs.size()) != params.d) {
    return false;
  }

  const std::size_t base_root_count = proof.commitments.roots_by_level.size();
  if (base_root_count == 1U) {
    if (proof.commitments.roots_by_level[0] != commitment_C) {
      return false;
    }
  } else if (base_root_count == static_cast<std::size_t>(params.d + 1)) {
    if (proof.commitments.roots_by_level[static_cast<std::size_t>(params.d)] !=
        commitment_C) {
      return false;
    }
  } else if (base_root_count != 0U) {
    return false;
  }

  HashTranscript transcript = basefold_pcs_internal::MakeBaseFoldTranscript();
  basefold_pcs_internal::AbsorbPublicInput(transcript, commitment_C, z, claimed_y);
  AbsorbChallengeConfig(transcript, challenge_cfg);
  AbsorbExtensionQuadraticPoly(
      transcript,
      proof.extension.h_by_level[static_cast<std::size_t>(params.d - 1)],
      extension_modulus);

  std::vector<ZZ_pEX> r_by_level(static_cast<std::size_t>(params.d));
  for (long i = params.d; i-- > 0;) {
    const ZZ_pEX r_i =
        SampleExtensionChallenge(transcript, "r/" + std::to_string(i),
                                 extension_modulus);
    r_by_level[static_cast<std::size_t>(i)] = r_i;
    transcript.AbsorbDigest(
        proof.extension.roots_by_level[static_cast<std::size_t>(i)]);
    if (i > 0) {
      AbsorbExtensionQuadraticPoly(
          transcript,
          proof.extension.h_by_level[static_cast<std::size_t>(i - 1)],
          extension_modulus);
    }
  }

  if (!CheckExtensionSumcheckRelations(proof.extension.h_by_level, r_by_level,
                                       claimed_y, extension_modulus)) {
    return false;
  }

  const ZZ_pEX r0 = r_by_level[0];
  const ZZ_pEX h1_r0 = EvalExtensionQuadraticPoly(proof.extension.h_by_level[0],
                                                  r0, extension_modulus);

  ZZ_pEX suffix_eq = ExtensionOne();
  for (long i = 0; i < params.d; ++i) {
    const ZZ_pEX zi_ext =
        LiftBaseToExtension(z[static_cast<std::size_t>(kappa + i)]);
    suffix_eq = MulExtension(
        suffix_eq,
        EqFactorExtension(zi_ext, r_by_level[static_cast<std::size_t>(i)],
                          extension_modulus),
        extension_modulus);
  }

  const std::vector<ZZ_pEX> f_eval =
      BooleanEvalTableFromMonomialCoeffsExtension(proof.extension.msg0_coeffs,
                                                  kappa, extension_modulus);

  std::vector<ZZ_pEX> prefix_eq(f_eval.size());
  if (kappa == 0) {
    prefix_eq[0] = ExtensionOne();
  } else {
    prefix_eq.resize(static_cast<std::size_t>(1L << kappa));
    prefix_eq[0] = ExtensionOne();
    const ZZ_pEX one = ExtensionOne();
    for (long var = 0; var < kappa; ++var) {
      const long old = 1L << var;
      const ZZ_pEX zi = LiftBaseToExtension(z[static_cast<std::size_t>(var)]);
      const ZZ_pEX f0 = SubExtension(one, zi, extension_modulus);
      const ZZ_pEX f1 = zi;
      for (long mask = 0; mask < old; ++mask) {
        const ZZ_pEX base = prefix_eq[static_cast<std::size_t>(mask)];
        prefix_eq[static_cast<std::size_t>(mask)] =
            MulExtension(base, f0, extension_modulus);
        prefix_eq[static_cast<std::size_t>(mask + old)] =
            MulExtension(base, f1, extension_modulus);
      }
    }
  }

  ZZ_pEX sum = ExtensionZero();
  for (long mask = 0; mask < static_cast<long>(f_eval.size()); ++mask) {
    sum = AddExtension(sum,
                       MulExtension(f_eval[static_cast<std::size_t>(mask)],
                                    prefix_eq[static_cast<std::size_t>(mask)],
                                    extension_modulus),
                       extension_modulus);
  }
  if (MulExtension(suffix_eq, sum, extension_modulus) != h1_r0) {
    return false;
  }

  const std::vector<ZZ_pEX> expected_pi0 =
      EncodeC0Extension(proof.extension.msg0_coeffs, params, extension_modulus);
  if (expected_pi0 != proof.extension.pi0_codeword) {
    return false;
  }
  if (proof.extension.roots_by_level[0] !=
      ExtensionMerkleTree::Build(proof.extension.pi0_codeword, extension_modulus)
          .Root()) {
    return false;
  }

  const long n_last = CodewordLengthAtLevel(params, params.d - 1);
  const long n_d = CodewordLengthAtLevel(params, params.d);
  std::vector<IOPPQueryPlan> query_plans(static_cast<std::size_t>(num_queries));
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    query_plans[static_cast<std::size_t>(q)] = MakeQueryPlan(mu, params);
  }

  const std::vector<std::vector<long>> requested_indices_by_tree =
      basefold_pcs_internal::CollectBaseQueryIndicesByTree(query_plans, params);
  const MerkleMultiproof &top_query_multiproof =
      proof.extension.base_top_query_multiproof;
  const std::vector<long> &expected_top_indices =
      requested_indices_by_tree[static_cast<std::size_t>(params.d)];
  if (static_cast<long>(top_query_multiproof.values.length()) !=
      static_cast<long>(expected_top_indices.size())) {
    return false;
  }
  if (!MerkleVerifyMultiproof(commitment_C, n_d, expected_top_indices,
                              top_query_multiproof)) {
    return false;
  }

  ScopedTimer query_timer(
      prof ? &prof->ext_verify_query_merkle_ns : nullptr,
      prof ? &prof->ext_verify_query_merkle_calls : nullptr);
  multiproof_replay::ValuePositionCache top_value_positions;
  if (!multiproof_replay::BuildValuePositionCache(
          expected_top_indices, top_query_multiproof.values.length(),
          &top_value_positions)) {
    return false;
  }
  std::vector<multiproof_replay::ValuePositionCache> value_positions_by_tree(
      static_cast<std::size_t>(params.d));
  for (long tree_level = 0; tree_level < params.d; ++tree_level) {
    const ExtensionMerkleMultiproof &multiproof =
        proof.extension.query_multiproofs[static_cast<std::size_t>(tree_level)];
    const std::vector<long> &expected_indices =
        requested_indices_by_tree[static_cast<std::size_t>(tree_level)];
    if (multiproof.values.size() != expected_indices.size()) {
      return false;
    }
    if (!multiproof_replay::BuildValuePositionCache(
            expected_indices, static_cast<long>(multiproof.values.size()),
            &value_positions_by_tree[static_cast<std::size_t>(tree_level)])) {
      return false;
    }
    const long leaf_count = CodewordLengthAtLevel(params, tree_level);
    if (!VerifyExtensionMerkleMultiproof(
            proof.extension.roots_by_level[static_cast<std::size_t>(tree_level)],
            leaf_count, expected_indices, multiproof, extension_modulus)) {
      return false;
    }
  }

  auto verify_one_query = [&](long q) -> bool {
    Profile *query_prof = ActiveProfile();
    ScopedTimer query_timer(
        query_prof ? &query_prof->ext_verify_query_merkle_ns : nullptr,
        query_prof ? &query_prof->ext_verify_query_merkle_calls : nullptr);

    const IOPPQueryPlan &plan = query_plans[static_cast<std::size_t>(q)];
    const long top_i = params.d - 1;
    const long mu_top = plan.mu_by_level[static_cast<std::size_t>(top_i)];
    const long n_top = CodewordLengthAtLevel(params, top_i);

    const FieldElement *left_top = multiproof_replay::FindMultiproofValue(
        top_value_positions, mu_top,
        [&](std::size_t pos) {
          return &top_query_multiproof.values[static_cast<long>(pos)];
        });
    const FieldElement *right_top = multiproof_replay::FindMultiproofValue(
        top_value_positions, mu_top + n_top, [&](std::size_t pos) {
          return &top_query_multiproof.values[static_cast<long>(pos)];
        });
    if (left_top == nullptr || right_top == nullptr) {
      return false;
    }

    const FieldElement left_top_value = *left_top;
    const FieldElement right_top_value = *right_top;

    for (long i = params.d; i-- > 0;) {
      const long mu_i = plan.mu_by_level[static_cast<std::size_t>(i)];
      const long n_i = CodewordLengthAtLevel(params, i);
      if (mu_i < 0 || mu_i >= n_i) {
        return false;
      }

      ZZ_pEX left_ext_value;
      ZZ_pEX right_ext_value;
      ZZ_pEX folded_ext_value;

      const ExtensionMerkleMultiproof &folded_multiproof =
          proof.extension.query_multiproofs[static_cast<std::size_t>(i)];
      const multiproof_replay::ValuePositionCache &folded_positions =
          value_positions_by_tree[static_cast<std::size_t>(i)];
      const ZZ_pEX *folded_ext = multiproof_replay::FindMultiproofValue(
          folded_positions, mu_i,
          [&](std::size_t pos) { return &folded_multiproof.values[pos]; });
      if (folded_ext == nullptr) {
        return false;
      }
      folded_ext_value = *folded_ext;

      if (i < params.d - 1) {
        const ExtensionMerkleMultiproof &next_multiproof =
            proof.extension.query_multiproofs[static_cast<std::size_t>(i + 1)];
        const multiproof_replay::ValuePositionCache &next_positions =
            value_positions_by_tree[static_cast<std::size_t>(i + 1)];
        const ZZ_pEX *left_ext = multiproof_replay::FindMultiproofValue(
            next_positions, mu_i,
            [&](std::size_t pos) { return &next_multiproof.values[pos]; });
        const ZZ_pEX *right_ext = multiproof_replay::FindMultiproofValue(
            next_positions, mu_i + n_i,
            [&](std::size_t pos) { return &next_multiproof.values[pos]; });
        if (left_ext == nullptr || right_ext == nullptr) {
          return false;
        }
        left_ext_value = *left_ext;
        right_ext_value = *right_ext;
      } else {
        left_ext_value = LiftBaseToExtension(left_top_value);
        right_ext_value = LiftBaseToExtension(right_top_value);
      }

      const FieldElement &t = params.diag_T[static_cast<std::size_t>(i)][mu_i];
      const FieldElement x1 = t;
      const FieldElement x2 = params.zeta * t;
      const ZZ_pEX expected_folded = EvalLineAtExtension(
          r_by_level[static_cast<std::size_t>(i)], x1, left_ext_value, x2,
          right_ext_value, extension_modulus);
      if (expected_folded != folded_ext_value) {
        return false;
      }

      if (i > 0) {
        const long n_prev = CodewordLengthAtLevel(params, i - 1);
        const long mu_prev = plan.mu_by_level[static_cast<std::size_t>(i - 1)];
        const long left_carry_index = mu_prev;
        const long right_carry_index = mu_prev + n_prev;
        const bool carries_left = (mu_i == left_carry_index);
        const bool carries_right = (mu_i == right_carry_index);
        if (!carries_left && !carries_right) {
          return false;
        }

        const ZZ_pEX *next_value = multiproof_replay::FindMultiproofValue(
            folded_positions, carries_left ? left_carry_index : right_carry_index,
            [&](std::size_t pos) { return &folded_multiproof.values[pos]; });
        if (next_value == nullptr || folded_ext_value != *next_value) {
          return false;
        }
      } else if (folded_ext_value !=
                 proof.extension.pi0_codeword[static_cast<std::size_t>(mu_i)]) {
        return false;
      }
    }

    return true;
  };

  return basefold_pcs_internal::VerifyQueriesMaybeParallel(num_queries, prof,
                                                           verify_one_query);
}

}  // namespace

BaseFoldPCSEvalProof BaseFoldPCSProveEvalWithChallengeConfig(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  ValidateChallengeConfigOrThrow(challenge_cfg, params);
  if (!challenge_cfg.use_extension_challenges) {
    return BaseFoldPCSProveEval(f_coeffs, z, claimed_y, num_queries, params);
  }

  basefold_pcs_internal::ValidateParamsOrThrow(params);
  const long kappa = basefold_pcs_internal::Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim) {
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfig: z has wrong dimension");
  }
  if (f_coeffs.length() != MessageLength(params)) {
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfig: f_coeffs has wrong length");
  }
  if (num_queries < 0) {
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfig: num_queries must be non-negative");
  }
  if (EvalMultilinearMonomialCoeffs(f_coeffs, z) != claimed_y) {
    LogicError("BaseFoldPCSProveEvalWithChallengeConfig: claimed_y != f(z)");
  }

  const BaseFoldPCSCommitArtifacts commit_artifacts =
      BaseFoldPCSBuildCommitArtifactsUnchecked(f_coeffs, params);
  return BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, commit_artifacts,
      challenge_cfg);
}

BaseFoldPCSEvalProof BaseFoldPCSProveEvalWithChallengeConfigUnchecked(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  ValidateChallengeConfigOrThrow(challenge_cfg, params);
  if (!challenge_cfg.use_extension_challenges) {
    return BaseFoldPCSProveEvalUnchecked(f_coeffs, z, claimed_y, num_queries,
                                         params);
  }

  const BaseFoldPCSCommitArtifacts commit_artifacts =
      BaseFoldPCSBuildCommitArtifactsUnchecked(f_coeffs, params);
  return BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, commit_artifacts,
      challenge_cfg);
}

BaseFoldPCSEvalProof
BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  ValidateChallengeConfigOrThrow(challenge_cfg, params);
  if (!challenge_cfg.use_extension_challenges) {
    return BaseFoldPCSProveEvalFromCommittedTopOracle(
        f_coeffs, z, claimed_y, num_queries, params, commit_artifacts);
  }

  basefold_pcs_internal::ValidateParamsOrThrow(params);
  const long kappa = basefold_pcs_internal::Log2ExactPowerOfTwoLong(params.k0);
  const long point_dim = params.d + kappa;
  if (static_cast<long>(z.size()) != point_dim) {
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle: z has wrong dimension");
  }
  if (f_coeffs.length() != MessageLength(params)) {
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle: f_coeffs has wrong length");
  }
  if (num_queries < 0) {
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle: num_queries must be non-negative");
  }
  if (EvalMultilinearMonomialCoeffs(f_coeffs, z) != claimed_y) {
    LogicError(
        "BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle: claimed_y != f(z)");
  }

  return BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, commit_artifacts,
      challenge_cfg);
}

BaseFoldPCSEvalProof
BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params,
    const BaseFoldPCSCommitArtifacts &commit_artifacts,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  ValidateChallengeConfigOrThrow(challenge_cfg, params);
  if (!challenge_cfg.use_extension_challenges) {
    return BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
        f_coeffs, z, claimed_y, num_queries, params, commit_artifacts);
  }
  return ProveEvalWithExtensionChallengesFromCommittedTopOracleUnchecked(
      f_coeffs, z, claimed_y, num_queries, params, commit_artifacts,
      challenge_cfg);
}

bool BaseFoldPCSVerifyEvalWithChallengeConfig(
    const MerkleRoot &commitment_C, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const BaseFoldPCSEvalProof &proof, const FoldableCodeParams &params,
    const BaseFoldPCSChallengeConfig &challenge_cfg) {
  ValidateChallengeConfigOrThrow(challenge_cfg, params);
  if (!challenge_cfg.use_extension_challenges) {
    return BaseFoldPCSVerifyEval(commitment_C, z, claimed_y, num_queries, proof,
                                 params);
  }
  return VerifyEvalWithExtensionChallenges(commitment_C, z, claimed_y,
                                           num_queries, proof, params,
                                           challenge_cfg);
}

}  // namespace basefold
