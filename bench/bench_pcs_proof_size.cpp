#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pX.h>

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "BaseFold/BaseFoldPCS.hpp"
#include "BaseFold/Multilinear.hpp"
#include "BaseFold/ProofSize.hpp"
#include "GaloisRing/PrimitiveElement.hpp"

using NTL::conv;
using NTL::LogicError;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEX;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;

namespace {

struct ContextSpec {
  std::string label;
  ZZ mod = ZZ(0);      // ZZ_p modulus (p for fields, p^s for rings)
  ZZ prime_p = ZZ(0);  // optional: the prime p (only used by unit checks)
  std::vector<ZZ> F_coeffs;     // extension modulus polynomial coefficients
  std::vector<ZZ> zeta_coeffs;  // ζ element coefficients
  // Coefficients for challenge extension modulus E(U), represented as
  // "a0;a1;...;ad", where each ai is a ZZ_pE element written "c0,c1,...".
  // Empty means "use default E(U)"; degree is 2 unless challenge_ext_degree is set.
  std::vector<std::vector<ZZ>> challenge_ext_coeffs;
  long challenge_ext_degree = 0;  // optional default degree for E(U)
};

long Pow2Checked(long e) {
  if (e < 0) LogicError("Pow2Checked: negative exponent");
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
}

bool IsPowerOfTwoLong(long n) { return n > 0 && (n & (n - 1)) == 0; }

long Log2ExactPowerOfTwoLong(long n) {
  if (!IsPowerOfTwoLong(n)) LogicError("Log2ExactPowerOfTwoLong: not a power of two");
  long d = 0;
  while (n > 1) {
    n >>= 1;
    ++d;
  }
  return d;
}

bool ParseZZString(const std::string &s, ZZ &out) {
  if (s.empty())
    return false;
  std::size_t pos = 0;
  bool neg = false;
  if (s[pos] == '+' || s[pos] == '-') {
    neg = (s[pos] == '-');
    ++pos;
  }
  if (pos >= s.size())
    return false;

  ZZ v(0);
  for (; pos < s.size(); ++pos) {
    const unsigned char ch = static_cast<unsigned char>(s[pos]);
    if (!std::isdigit(ch))
      return false;
    v *= 10;
    v += static_cast<long>(ch - static_cast<unsigned char>('0'));
  }
  out = neg ? -v : v;
  return true;
}

std::vector<ZZ> ParseCoeffList(const std::string &s) {
  std::vector<ZZ> out;
  std::size_t pos = 0;
  while (pos < s.size()) {
    const std::size_t comma = s.find(',', pos);
    const std::size_t end = (comma == std::string::npos) ? s.size() : comma;
    std::string token = s.substr(pos, end - pos);

    const std::size_t first = token.find_first_not_of(" \t");
    const std::size_t last = token.find_last_not_of(" \t");
    if (first == std::string::npos)
      LogicError("ParseCoeffList: empty coefficient");
    token = token.substr(first, last - first + 1);

    ZZ v;
    if (!ParseZZString(token, v)) LogicError("ParseCoeffList: bad integer token");

    out.push_back(v);
    pos = (comma == std::string::npos) ? s.size() : (comma + 1);
  }
  if (out.empty())
    LogicError("ParseCoeffList: empty list");
  return out;
}

std::vector<std::vector<ZZ>> ParseNestedCoeffList(const std::string &s) {
  std::vector<std::vector<ZZ>> out;
  std::size_t pos = 0;
  while (pos < s.size()) {
    const std::size_t semi = s.find(';', pos);
    const std::size_t end = (semi == std::string::npos) ? s.size() : semi;
    std::string token = s.substr(pos, end - pos);

    const std::size_t first = token.find_first_not_of(" \t");
    const std::size_t last = token.find_last_not_of(" \t");
    if (first == std::string::npos)
      LogicError("ParseNestedCoeffList: empty coefficient block");
    token = token.substr(first, last - first + 1);

    out.push_back(ParseCoeffList(token));
    pos = (semi == std::string::npos) ? s.size() : (semi + 1);
  }
  if (out.empty())
    LogicError("ParseNestedCoeffList: empty list");
  return out;
}

ZZ_pX BuildZZpX(const std::vector<ZZ> &coeffs) {
  ZZ_pX poly;
  NTL::clear(poly);
  for (std::size_t i = 0; i < coeffs.size(); ++i) {
    if (coeffs[i] != 0) {
      SetCoeff(poly, static_cast<long>(i), conv<ZZ_p>(coeffs[i]));
    }
  }
  return poly;
}

ZZ_pE BuildZZpE(const std::vector<ZZ> &coeffs) {
  const ZZ_pX poly = BuildZZpX(coeffs);
  ZZ_pE out;
  conv(out, poly);
  return out;
}

ZZ_pEX BuildZZpEX(const std::vector<std::vector<ZZ>> &coeffs) {
  if (coeffs.empty()) {
    LogicError("BuildZZpEX: empty coefficient list");
  }
  ZZ_pEX poly;
  NTL::clear(poly);
  for (std::size_t i = 0; i < coeffs.size(); ++i) {
    const ZZ_pE c = BuildZZpE(coeffs[i]);
    if (c != 0) {
      NTL::SetCoeff(poly, static_cast<long>(i), c);
    }
  }
  poly.normalize();
  return poly;
}

ZZ NormalizeMod(const ZZ &x, const ZZ &mod) {
  if (mod <= 0) LogicError("NormalizeMod: mod must be positive");
  ZZ r = x % mod;
  if (r < 0) r += mod;
  return r;
}

void DeduceBasePrimeAndExponent(const ContextSpec &spec, ZZ &p_out, long &k_out) {
  if (spec.mod <= 1) LogicError("DeduceBasePrimeAndExponent: mod must be > 1");

  if (spec.prime_p > 1) {
    ZZ m = spec.mod;
    long k = 0;
    while ((m % spec.prime_p) == 0) {
      m /= spec.prime_p;
      ++k;
    }
    if (k <= 0 || m != 1) {
      LogicError("DeduceBasePrimeAndExponent: mod must equal prime_p^k");
    }
    p_out = spec.prime_p;
    k_out = k;
    return;
  }

  p_out = spec.mod;
  k_out = 1;
}

void ValidateMonic(const std::vector<ZZ> &coeffs, const ZZ &mod,
                   const char *what) {
  if (coeffs.empty()) {
    const std::string msg = std::string(what) + ": empty polynomial";
    LogicError(msg.c_str());
  }
  long last = static_cast<long>(coeffs.size()) - 1;
  while (last > 0 &&
         NormalizeMod(coeffs[static_cast<std::size_t>(last)], mod) == 0) {
    --last;
  }
  if (last <= 0) {
    const std::string msg = std::string(what) + ": degree must be >= 1";
    LogicError(msg.c_str());
  }
  const ZZ lead = NormalizeMod(coeffs[static_cast<std::size_t>(last)], mod);
  if (lead != 1) {
    const std::string msg =
        std::string(what) + ": leading coefficient must be 1 (monic)";
    LogicError(msg.c_str());
  }
}

std::uint64_t FixedCoeffByteWidth(const ZZ &mod) {
  if (mod <= 1)
    LogicError("FixedCoeffByteWidth: mod must be > 1");
  const long bits = NTL::NumBits(mod - 1);
  return static_cast<std::uint64_t>((bits + 7) / 8);
}

std::uint64_t FixedFieldElementBytes(const ZZ &mod, long ext_degree) {
  if (ext_degree <= 0)
    LogicError("FixedFieldElementBytes: ext_degree must be > 0");
  const std::uint64_t coeff_bytes = FixedCoeffByteWidth(mod);
  unsigned __int128 total = 0;
  total += static_cast<unsigned __int128>(ext_degree) * coeff_bytes;
  if (total > std::numeric_limits<std::uint64_t>::max()) {
    LogicError("FixedFieldElementBytes: overflow");
  }
  return static_cast<std::uint64_t>(total);
}

std::uint64_t FixedExtensionElementBytes(std::uint64_t field_elem_bytes,
                                         long challenge_ext_degree) {
  if (challenge_ext_degree <= 0) {
    LogicError("FixedExtensionElementBytes: challenge_ext_degree must be > 0");
  }
  unsigned __int128 total = 0;
  total += static_cast<unsigned __int128>(challenge_ext_degree) *
           field_elem_bytes;
  if (total > std::numeric_limits<std::uint64_t>::max()) {
    LogicError("FixedExtensionElementBytes: overflow");
  }
  return static_cast<std::uint64_t>(total);
}

std::size_t MerkleHeight(std::uint64_t leaf_count) {
  if (leaf_count <= 1)
    return 0;
  std::size_t height = 0;
  std::uint64_t n = leaf_count;
  while (n > 1) {
    if (n & 1ULL)
      n += 1;
    n /= 2;
    ++height;
  }
  return height;
}

std::uint64_t MerkleOpeningBytes(std::uint64_t leaf_count,
                                 std::uint64_t field_elem_bytes) {
  static constexpr std::uint64_t kHashBytes = 32;

  const std::size_t height = MerkleHeight(leaf_count);
  unsigned __int128 total = 0;
  total += field_elem_bytes;
  total += static_cast<unsigned __int128>(height) * kHashBytes;  // sibling hashes

  if (total > std::numeric_limits<std::uint64_t>::max())
    LogicError("MerkleOpeningBytes: overflow");
  return static_cast<std::uint64_t>(total);
}

std::uint64_t EstimateEvalProofSizeFormulaBytes(const ZZ &mod, long c, long d,
                                                long k0, long num_queries,
                                                long ext_degree) {
  static constexpr std::uint64_t kHashBytes = 32;
  if (c <= 0)
    LogicError("EstimateEvalProofSizeFormulaBytes: c must be > 0");
  if (d < 0)
    LogicError("EstimateEvalProofSizeFormulaBytes: d must be >= 0");
  if (k0 <= 0)
    LogicError("EstimateEvalProofSizeFormulaBytes: k0 must be > 0");
  if (!IsPowerOfTwoLong(k0))
    LogicError("EstimateEvalProofSizeFormulaBytes: k0 must be a power of two");
  if (num_queries < 0)
    LogicError("EstimateEvalProofSizeFormulaBytes: queries must be >= 0");

  const std::uint64_t fe_bytes = FixedFieldElementBytes(mod, ext_degree);

  unsigned __int128 total = 0;
  total += static_cast<unsigned __int128>(d + 1) * kHashBytes;          // roots
  total += static_cast<unsigned __int128>(d) * 3ULL * fe_bytes;         // sumcheck polys
  total += static_cast<unsigned __int128>(c) *
           static_cast<unsigned __int128>(k0) * fe_bytes;               // pi0_full

  std::vector<std::uint64_t> open_bytes;
  open_bytes.resize(static_cast<std::size_t>(d + 1));

  std::uint64_t leaf_count =
      static_cast<std::uint64_t>(c) * static_cast<std::uint64_t>(k0);
  for (long level = 0; level <= d; ++level) {
    open_bytes[static_cast<std::size_t>(level)] =
        MerkleOpeningBytes(leaf_count, fe_bytes);
    if (level < d) {
      if (leaf_count > std::numeric_limits<std::uint64_t>::max() / 2ULL) {
        LogicError("EstimateEvalProofSizeFormulaBytes: overflow in n_i");
      }
      leaf_count *= 2ULL;
    }
  }

  unsigned __int128 per_query = 0;
  for (long i = 0; i < d; ++i) {
    per_query += open_bytes[static_cast<std::size_t>(i)];          // folded at level i
    per_query += 2ULL * open_bytes[static_cast<std::size_t>(i + 1)];  // left+right at level i+1
  }

  total += static_cast<unsigned __int128>(num_queries) * per_query;

  if (total > std::numeric_limits<std::uint64_t>::max())
    LogicError("EstimateEvalProofSizeFormulaBytes: overflow");
  return static_cast<std::uint64_t>(total);
}

std::uint64_t EstimateEvalProofSizeFormulaBytesExtensionChallenges(
    const ZZ &mod, long c, long d, long k0, long num_queries,
    long base_ext_degree, long challenge_ext_degree) {
  static constexpr std::uint64_t kHashBytes = 32;
  if (c <= 0)
    LogicError(
        "EstimateEvalProofSizeFormulaBytesExtensionChallenges: c must be > 0");
  if (d < 0)
    LogicError(
        "EstimateEvalProofSizeFormulaBytesExtensionChallenges: d must be >= 0");
  if (k0 <= 0)
    LogicError(
        "EstimateEvalProofSizeFormulaBytesExtensionChallenges: k0 must be > 0");
  if (!IsPowerOfTwoLong(k0))
    LogicError("EstimateEvalProofSizeFormulaBytesExtensionChallenges: k0 must "
               "be a power of two");
  if (num_queries < 0)
    LogicError("EstimateEvalProofSizeFormulaBytesExtensionChallenges: queries "
               "must be >= 0");
  if (base_ext_degree <= 0)
    LogicError("EstimateEvalProofSizeFormulaBytesExtensionChallenges: "
               "base_ext_degree must be > 0");
  if (challenge_ext_degree <= 0)
    LogicError("EstimateEvalProofSizeFormulaBytesExtensionChallenges: "
               "challenge_ext_degree must be > 0");

  const std::uint64_t fe_bytes = FixedFieldElementBytes(mod, base_ext_degree);
  const std::uint64_t ext_fe_bytes =
      FixedExtensionElementBytes(fe_bytes, challenge_ext_degree);

  unsigned __int128 total = 0;
  const std::uint64_t n0 =
      static_cast<std::uint64_t>(c) * static_cast<std::uint64_t>(k0);

  std::vector<std::uint64_t> n_by_level;
  n_by_level.resize(static_cast<std::size_t>(d + 1));
  std::uint64_t leaf_count = n0;
  for (long level = 0; level <= d; ++level) {
    n_by_level[static_cast<std::size_t>(level)] = leaf_count;
    if (level < d) {
      if (leaf_count > std::numeric_limits<std::uint64_t>::max() / 2ULL) {
        LogicError(
            "EstimateEvalProofSizeFormulaBytesExtensionChallenges: overflow in "
            "n_i");
      }
      leaf_count *= 2ULL;
    }
  }

  total += static_cast<unsigned __int128>(d + 1) * kHashBytes;  // roots
  total += static_cast<unsigned __int128>(d) * 3ULL *
           ext_fe_bytes;  // sumcheck polys
  total += static_cast<unsigned __int128>(n0) * ext_fe_bytes;  // pi0_full

  unsigned __int128 per_query_ext = 0;
  for (long i = 0; i < d; ++i) {
    const std::uint64_t n_i = n_by_level[static_cast<std::size_t>(i)];
    const std::uint64_t n_ip1 = n_by_level[static_cast<std::size_t>(i + 1)];
    per_query_ext += MerkleOpeningBytes(n_i, ext_fe_bytes);  // folded at level i
    per_query_ext +=
        2ULL * static_cast<unsigned __int128>(MerkleOpeningBytes(
                   n_ip1, ext_fe_bytes));  // left+right at level i+1
  }
  total += static_cast<unsigned __int128>(num_queries) * per_query_ext;

  if (total > std::numeric_limits<std::uint64_t>::max())
    LogicError("EstimateEvalProofSizeFormulaBytesExtensionChallenges: overflow");
  return static_cast<std::uint64_t>(total);
}

basefold::FoldableCodeParams BuildParams_k0_1(long c, long d, const ZZ &prime_p,
                                             const ZZ_pE &zeta) {
  if (c <= 0) LogicError("BuildParams_k0_1: c must be positive");
  if (d < 0) LogicError("BuildParams_k0_1: d must be non-negative");

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = 1;
  params.d = d;
  params.p = prime_p;
  params.zeta = zeta;

  const long n0 = c;
  params.G0.SetDims(/*rows=*/1, /*cols=*/n0);
  params.G0[0][0] = ZZ_pE(1);
  if (n0 > 1) {
    params.G0[0][1] = zeta;
  }
  for (long j = 2; j < n0; ++j) {
    params.G0[0][j] = ZZ_pE(1);
  }

  params.diag_T.resize(static_cast<std::size_t>(d));
  for (long level = 0; level < d; ++level) {
    const long pow2 = Pow2Checked(level);
    if (c > std::numeric_limits<long>::max() / pow2) {
      LogicError("BuildParams_k0_1: overflow in n_i");
    }
    const long ni = c * pow2;
    params.diag_T[static_cast<std::size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      params.diag_T[static_cast<std::size_t>(level)][i] = ZZ_pE(1);
    }
  }

  return params;
}

NTL::mat_ZZ_pE BuildSystematicG0(long c, long k0) {
  if (c <= 0) LogicError("BuildSystematicG0: c must be positive");
  if (k0 <= 0) LogicError("BuildSystematicG0: k0 must be positive");
  if (c > std::numeric_limits<long>::max() / k0)
    LogicError("BuildSystematicG0: overflow in n0");
  const long n0 = c * k0;

  NTL::mat_ZZ_pE G0;
  G0.SetDims(k0, n0);

  const ZZ_pE one = ZZ_pE(1);
  for (long block = 0; block < c; ++block) {
    const long base = block * k0;
    for (long r = 0; r < k0; ++r) {
      G0[r][base + r] = one;
    }
  }
  return G0;
}

basefold::FoldableCodeParams BuildParams_k0_pow2(long c, long k0, long d,
                                                 const ZZ &prime_p,
                                                 const ZZ_pE &zeta) {
  if (c <= 0) LogicError("BuildParams_k0_pow2: c must be positive");
  if (k0 <= 0) LogicError("BuildParams_k0_pow2: k0 must be positive");
  if (!IsPowerOfTwoLong(k0))
    LogicError("BuildParams_k0_pow2: k0 must be a power of two");
  if (d < 0) LogicError("BuildParams_k0_pow2: d must be non-negative");

  if (c > std::numeric_limits<long>::max() / k0)
    LogicError("BuildParams_k0_pow2: overflow in n0");
  const long n0 = c * k0;

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = prime_p;
  params.zeta = zeta;

  params.G0 = BuildSystematicG0(c, k0);

  params.diag_T.resize(static_cast<std::size_t>(d));
  for (long level = 0; level < d; ++level) {
    const long pow2 = Pow2Checked(level);
    if (n0 > std::numeric_limits<long>::max() / pow2) {
      LogicError("BuildParams_k0_pow2: overflow in n_i");
    }
    const long ni = n0 * pow2;
    params.diag_T[static_cast<std::size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      params.diag_T[static_cast<std::size_t>(level)][i] = ZZ_pE(1);
    }
  }

  return params;
}

ZZ_pEX BuildChallengeExtensionModulus(const ContextSpec &spec,
                                      const ZZ_pE &zeta) {
  auto BuildDefault = [&](long degree) -> ZZ_pEX {
    if (degree <= 0) {
      LogicError("BuildChallengeExtensionModulus: degree must be > 0");
    }
    ZZ_pEX challenge_modulus;
    NTL::clear(challenge_modulus);
    NTL::SetCoeff(challenge_modulus, 0, zeta);
    if (degree > 1) {
      NTL::SetCoeff(challenge_modulus, 1, ZZ_pE(1));
    }
    NTL::SetCoeff(challenge_modulus, degree, ZZ_pE(1));
    return challenge_modulus;
  };

  if (!spec.challenge_ext_coeffs.empty()) {
    return BuildZZpEX(spec.challenge_ext_coeffs);
  }
  if (spec.challenge_ext_degree > 0) {
    return BuildDefault(spec.challenge_ext_degree);
  }
  return BuildDefault(2);
}

std::uint64_t SplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

ZZ ZZFromU64(std::uint64_t x) {
  static const ZZ kTwo32 = NTL::power2_ZZ(32);
  const long lo = static_cast<long>(x & 0xffffffffULL);
  const long hi = static_cast<long>((x >> 32) & 0xffffffffULL);
  ZZ out = to_ZZ(hi);
  out *= kTwo32;
  out += to_ZZ(lo);
  return out;
}

ZZ DeterministicResidue(std::uint64_t seed, const ZZ &modulus) {
  const long bits = NTL::NumBits(modulus);
  const long blocks = std::max<long>(1, (bits + 63) / 64);
  ZZ v(0);
  std::uint64_t x = seed;
  for (long i = 0; i < blocks; ++i) {
    x = SplitMix64(x + 0x9e3779b97f4a7c15ULL +
                   static_cast<std::uint64_t>(i));
    v <<= 64;
    v += ZZFromU64(x);
  }
  return v % modulus;
}

ZZ_pE MakeDeterministicElement(std::uint64_t seed) {
  const long r = ZZ_pE::degree();
  if (r <= 0) LogicError("MakeDeterministicElement: invalid extension degree");
  const ZZ modulus = NTL::ZZ_p::modulus();
  if (modulus <= 1) LogicError("MakeDeterministicElement: invalid modulus");

  ZZ_pX poly;
  NTL::clear(poly);
  std::uint64_t x = SplitMix64(seed);
  for (long j = 0; j < r; ++j) {
    x = SplitMix64(x + static_cast<std::uint64_t>(j));
    const ZZ cj = DeterministicResidue(
        x ^ (static_cast<std::uint64_t>(j) << 32), modulus);
    SetCoeff(poly, j, conv<ZZ_p>(cj));
  }
  ZZ_pE elem;
  conv(elem, poly);
  return elem;
}

vec_ZZ_pE MakeDeterministicCoefficients(long coeff_count, std::uint64_t seed) {
  vec_ZZ_pE coeffs;
  coeffs.SetLength(coeff_count);
  for (long i = 0; i < coeff_count; ++i) {
    coeffs[i] = MakeDeterministicElement(seed ^ static_cast<std::uint64_t>(i));
  }
  return coeffs;
}

std::vector<ZZ_pE> MakeDeterministicPoint(long d, std::uint64_t seed) {
  if (d < 0) LogicError("MakeDeterministicPoint: negative d");
  std::vector<ZZ_pE> z;
  z.resize(static_cast<std::size_t>(d));
  for (long i = 0; i < d; ++i) {
    z[static_cast<std::size_t>(i)] =
        MakeDeterministicElement(seed + 0x12345678ULL +
                                 static_cast<std::uint64_t>(i) * 0x9e3779b9ULL);
  }
  return z;
}

bool ParseLong(const char *s, long &out) {
  if (!s) return false;
  try {
    std::size_t idx = 0;
    const long v = std::stol(std::string(s), &idx, 10);
    if (idx != std::string(s).size()) return false;
    out = v;
    return true;
  } catch (...) {
    return false;
  }
}

bool ParseZZ(const char *s, ZZ &out) {
  if (!s) return false;
  return ParseZZString(std::string(s), out);
}

std::uint64_t ParseU64OrDie(const char *s, const char *flag) {
  try {
    std::size_t idx = 0;
    const std::uint64_t v = std::stoull(std::string(s), &idx, 10);
    if (idx != std::string(s).size()) {
      std::cerr << "Invalid " << flag << "\n";
      std::exit(2);
    }
    return v;
  } catch (...) {
    std::cerr << "Invalid " << flag << "\n";
    std::exit(2);
  }
}

void PrintHelp() {
  std::cout
      << "bench_pcs_proof_size (estimate eval proof size)\n\n"
      << "Usage:\n"
      << "  bench_pcs_proof_size [--mode field|ring|both] [--c <int>] [--k0 <int>] [--d <int>]\n"
      << "                     [--queries <int>] [--seed <u64>] [--formula]\n"
      << "                     [--use-extension-challenges]\n"
      << "                     [--field-challenge-ext <a0;a1;...>] [--ring-challenge-ext <a0;a1;...>]\n"
      << "                     [--field-challenge-degree <int>] [--ring-challenge-degree <int>]\n"
      << "                     [--auto-zeta teich]\n"
      << "                     [--field-mod <decimal-int>] [--field-F <a0,a1,...>] [--field-zeta <b0,b1,...>]\n"
      << "                     [--ring-mod <decimal-int>]  [--ring-p <decimal-int>] [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n\n"
      << "Notes:\n"
      << "  KB is KiB (1024 bytes).\n"
      << "  By default, prover uses BaseFoldPCSProveEval (includes parameter/length checks and claimed_y == f(z)).\n"
      << "  With --use-extension-challenges, prover switches to BaseFoldPCSProveEvalWithChallengeConfig.\n\n"
      << "  With --auto-zeta teich, zeta is derived as a Teichmuller generator from (p,k,F);\n"
      << "  then --field-zeta/--ring-zeta are ignored.\n\n"
      << "  With --formula, it does NOT run the prover. It estimates proof size from (c,d,queries)\n"
      << "  assuming fixed-size serialization and sha256-based Merkle hashing.\n"
      << "  For extension-challenge mode, formula uses extension element width and the same\n"
      << "  roots/sumcheck/pi0/query-opening structure.\n\n"
      << "  --field-challenge-ext / --ring-challenge-ext use ';' to separate ZZ_pE\n"
      << "  coefficients and ',' for each ZZ_pE coefficient polynomial.\n"
      << "  Example: '0,1;1;1' means E(U)=x + U + U^2.\n\n"
      << "  --field-challenge-degree / --ring-challenge-degree set only the degree m,\n"
      << "  and auto-build default E(U)=zeta + U + U^m (m=1 uses E(U)=zeta + U).\n\n"
      << "  If both --*-challenge-ext and --*-challenge-degree are provided, it's an error.\n\n"
      << "  If --*-challenge-ext is omitted, default is E(U)=zeta + U + U^2.\n\n"
      << "  PCS Eval supports k0 = 2^κ. The multilinear point dimension is (d + κ).\n\n"
      << "Examples:\n"
      << "  # GF(2^2) with F(x)=x^2+x+1 and zeta=x\n"
      << "  bench_pcs_proof_size --mode field --field-mod 2 --field-F 1,1,1 --field-zeta 0,1 --d 16 --queries 4\n"
      << "  # GR(4,2) with the same extension polynomial and zeta=x\n"
      << "  bench_pcs_proof_size --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 16 --queries 4\n"
      << "  # Auto zeta from Teichmuller subgroup generator\n"
      << "  bench_pcs_proof_size --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --auto-zeta teich --d 16 --queries 4\n"
      << "  # Extension-challenge formula estimate\n"
      << "  bench_pcs_proof_size --mode field --field-mod 2 --field-F 1,1,1 --field-zeta 0,1 --use-extension-challenges --field-challenge-ext '0,1;1;1' --d 16 --queries 4 --formula\n"
      << "  # Extension-challenge with degree-only default modulus\n"
      << "  bench_pcs_proof_size --mode field --field-mod 2 --field-F 1,1,1 --field-zeta 0,1 --use-extension-challenges --field-challenge-degree 4 --d 16 --queries 4 --formula\n";
}

void RunOneContext(const ContextSpec &spec, long c, long k0, long d,
                   long num_queries, bool use_extension_challenges,
                   bool formula_only, bool auto_zeta_teich, std::uint64_t seed) {
  if (spec.mod <= 1) LogicError("RunOneContext: modulus must be > 1");
  if (c <= 0) LogicError("RunOneContext: c must be > 0");
  if (k0 <= 0) LogicError("RunOneContext: k0 must be > 0");
  if (!IsPowerOfTwoLong(k0)) LogicError("RunOneContext: k0 must be a power of two");
  if (d < 0) LogicError("RunOneContext: d must be >= 0");
  if (num_queries < 0) LogicError("RunOneContext: queries must be >= 0");

  const ZZ modulus = spec.mod;
  ZZ_pPush mod_push(modulus);

  ValidateMonic(spec.F_coeffs, spec.mod, "F");
  const ZZ_pX F = BuildZZpX(spec.F_coeffs);
  const long base_ext_degree = NTL::deg(F);
  ZZ_pEPush e_push(F);

  ZZ_pE zeta;
  if (auto_zeta_teich) {
    ZZ p_base;
    long k_base = 0;
    DeduceBasePrimeAndExponent(spec, p_base, k_base);
    const long s = NTL::deg(F);
    zeta = FindTeichmullerGenerator(p_base, k_base, s, F);
  } else {
    const ZZ_pX zpoly = BuildZZpX(spec.zeta_coeffs);
    conv(zeta, zpoly);
  }

  basefold::BaseFoldPCSChallengeConfig challenge_cfg;
  long challenge_degree = 0;
  const basefold::BaseFoldPCSChallengeConfig *challenge_cfg_ptr = nullptr;
  if (use_extension_challenges) {
    challenge_cfg.use_extension_challenges = true;
    challenge_cfg.challenge_extension_modulus =
        BuildChallengeExtensionModulus(spec, zeta);
    challenge_degree = NTL::deg(challenge_cfg.challenge_extension_modulus);
    if (challenge_degree <= 0) {
      LogicError(
          "RunOneContext: challenge extension modulus degree must be > 0");
    }
    challenge_cfg_ptr = &challenge_cfg;
  }

  if (formula_only) {
    const std::uint64_t bytes =
        use_extension_challenges
            ? EstimateEvalProofSizeFormulaBytesExtensionChallenges(
                  spec.mod, c, d, k0, num_queries, base_ext_degree,
                  challenge_degree)
            : EstimateEvalProofSizeFormulaBytes(spec.mod, c, d, k0, num_queries,
                                                base_ext_degree);
    const double kb = static_cast<double>(bytes) / 1024.0;

    std::cout << "\n[" << spec.label << "] c=" << c << " k0=" << k0 << " d=" << d
              << "  mod=" << spec.mod << "  queries=" << num_queries
              << "  (formula)";
    if (use_extension_challenges) {
      std::cout << "  ext_challenges=on"
                << "  ext_deg=" << challenge_degree;
    }
    std::cout << "\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  proof size  " << kb << " KB  (" << bytes << " B)\n";
    return;
  }

  const basefold::FoldableCodeParams params = (k0 == 1)
                                                  ? BuildParams_k0_1(
                                                        c, d,
                                                        (spec.prime_p > 1)
                                                            ? spec.prime_p
                                                            : spec.mod,
                                                        zeta)
                                                  : BuildParams_k0_pow2(
                                                        c, k0, d,
                                                        (spec.prime_p > 1)
                                                            ? spec.prime_p
                                                            : spec.mod,
                                                        zeta);

  const long pow2_d = Pow2Checked(d);
  if (k0 > std::numeric_limits<long>::max() / pow2_d) {
    LogicError("RunOneContext: overflow in k_d");
  }
  const long k_d = k0 * pow2_d;
  const vec_ZZ_pE f_coeffs = MakeDeterministicCoefficients(k_d, seed);
  const long point_dim = d + Log2ExactPowerOfTwoLong(k0);
  const std::vector<ZZ_pE> z =
      MakeDeterministicPoint(point_dim, seed ^ 0xdeadbeefULL);

  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  const basefold::BaseFoldPCSEvalProof proof =
      (challenge_cfg_ptr != nullptr)
          ? basefold::BaseFoldPCSProveEvalWithChallengeConfig(
                f_coeffs, z, y, num_queries, params, *challenge_cfg_ptr)
          : basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);

  const std::uint64_t bytes = basefold::BaseFoldPCSEvalProofSizeBytes(proof);
  const double kb = static_cast<double>(bytes) / 1024.0;

  std::cout << "\n[" << spec.label << "] c=" << c << " k0=" << k0 << " d=" << d
            << "  mod=" << spec.mod << "  queries=" << num_queries;
  if (use_extension_challenges) {
    std::cout << "  ext_challenges=on"
              << "  ext_deg=" << challenge_degree;
  }
  std::cout << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  proof size  " << kb << " KB  (" << bytes << " B)\n";
}

}  // namespace

int main(int argc, char **argv) {
  long d = 16;
  long c = 2;
  long k0 = 1;
  long num_queries = 4;
  bool use_extension_challenges = false;
  bool formula_only = false;
  std::uint64_t seed = 0;
  bool auto_zeta_teich = false;

  bool do_field = true;
  bool do_ring = true;

  ContextSpec field;
  field.label = "Field";
  field.mod = to_ZZ(2);
  field.prime_p = ZZ(0);
  field.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};   // x^2 + x + 1
  field.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};   // x

  ContextSpec ring;
  ring.label = "Ring";
  ring.mod = to_ZZ(4);
  ring.prime_p = to_ZZ(2);
  ring.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};    // x^2 + x + 1
  ring.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};    // x

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    auto NeedValue = [&](const char *flag) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << flag << "\n";
        std::exit(2);
      }
      return argv[++i];
    };

    if (arg == "--mode") {
      const std::string m = NeedValue("--mode");
      if (m == "field") {
        do_field = true;
        do_ring = false;
      } else if (m == "ring") {
        do_field = false;
        do_ring = true;
      } else if (m == "both") {
        do_field = true;
        do_ring = true;
      } else {
        std::cerr << "Invalid --mode (expected field|ring|both)\n";
        return 2;
      }
    } else if (arg == "--d") {
      if (!ParseLong(NeedValue("--d"), d) || d < 0) {
        std::cerr << "Invalid --d\n";
        return 2;
      }
    } else if (arg == "--c") {
      if (!ParseLong(NeedValue("--c"), c) || c <= 0) {
        std::cerr << "Invalid --c\n";
        return 2;
      }
    } else if (arg == "--k0") {
      if (!ParseLong(NeedValue("--k0"), k0) || k0 <= 0) {
        std::cerr << "Invalid --k0\n";
        return 2;
      }
    } else if (arg == "--queries") {
      if (!ParseLong(NeedValue("--queries"), num_queries) || num_queries < 0) {
        std::cerr << "Invalid --queries\n";
        return 2;
      }
    } else if (arg == "--formula") {
      formula_only = true;
    } else if (arg == "--use-extension-challenges") {
      use_extension_challenges = true;
    } else if (arg == "--seed") {
      seed = ParseU64OrDie(NeedValue("--seed"), "--seed");
    } else if (arg == "--auto-zeta") {
      const std::string mode = NeedValue("--auto-zeta");
      if (mode == "teich") {
        auto_zeta_teich = true;
      } else {
        std::cerr << "Invalid --auto-zeta (expected teich)\n";
        return 2;
      }
    } else if (arg == "--field-mod") {
      if (!ParseZZ(NeedValue("--field-mod"), field.mod) || field.mod <= 1) {
        std::cerr << "Invalid --field-mod\n";
        return 2;
      }
    } else if (arg == "--field-F") {
      field.F_coeffs = ParseCoeffList(NeedValue("--field-F"));
    } else if (arg == "--field-zeta") {
      field.zeta_coeffs = ParseCoeffList(NeedValue("--field-zeta"));
    } else if (arg == "--field-challenge-ext") {
      field.challenge_ext_coeffs =
          ParseNestedCoeffList(NeedValue("--field-challenge-ext"));
    } else if (arg == "--field-challenge-degree") {
      if (!ParseLong(NeedValue("--field-challenge-degree"),
                     field.challenge_ext_degree) ||
          field.challenge_ext_degree <= 0) {
        std::cerr << "Invalid --field-challenge-degree\n";
        return 2;
      }
    } else if (arg == "--ring-mod") {
      if (!ParseZZ(NeedValue("--ring-mod"), ring.mod) || ring.mod <= 1) {
        std::cerr << "Invalid --ring-mod\n";
        return 2;
      }
    } else if (arg == "--ring-p") {
      if (!ParseZZ(NeedValue("--ring-p"), ring.prime_p) || ring.prime_p <= 1) {
        std::cerr << "Invalid --ring-p\n";
        return 2;
      }
    } else if (arg == "--ring-F") {
      ring.F_coeffs = ParseCoeffList(NeedValue("--ring-F"));
    } else if (arg == "--ring-zeta") {
      ring.zeta_coeffs = ParseCoeffList(NeedValue("--ring-zeta"));
    } else if (arg == "--ring-challenge-ext") {
      ring.challenge_ext_coeffs =
          ParseNestedCoeffList(NeedValue("--ring-challenge-ext"));
    } else if (arg == "--ring-challenge-degree") {
      if (!ParseLong(NeedValue("--ring-challenge-degree"),
                     ring.challenge_ext_degree) ||
          ring.challenge_ext_degree <= 0) {
        std::cerr << "Invalid --ring-challenge-degree\n";
        return 2;
      }
    } else if (arg == "--help" || arg == "-h") {
      PrintHelp();
      return 0;
    } else {
      std::cerr << "Unknown arg: " << arg << "\n";
      return 2;
    }
  }

  try {
    if (!field.challenge_ext_coeffs.empty() && field.challenge_ext_degree > 0) {
      std::cerr << "Invalid options: --field-challenge-ext and "
                   "--field-challenge-degree are mutually exclusive\n";
      return 2;
    }
    if (!ring.challenge_ext_coeffs.empty() && ring.challenge_ext_degree > 0) {
      std::cerr << "Invalid options: --ring-challenge-ext and "
                   "--ring-challenge-degree are mutually exclusive\n";
      return 2;
    }
    if (!do_field && !do_ring) {
      std::cerr << "Nothing to do: --mode disabled both field and ring\n";
      return 2;
    }
    if (do_field)
      RunOneContext(field, c, k0, d, num_queries, use_extension_challenges,
                    formula_only,
                    auto_zeta_teich, seed);
    if (do_ring)
      RunOneContext(ring, c, k0, d, num_queries, use_extension_challenges,
                    formula_only,
                    auto_zeta_teich, seed);
  } catch (const std::exception &e) {
    std::cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    std::cerr << "Unhandled non-std exception\n";
    return 2;
  }

  return 0;
}
