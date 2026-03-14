#ifndef BASEFOLD_BENCH_COMMON_HELPERS_HPP_
#define BASEFOLD_BENCH_COMMON_HELPERS_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "PCS/BaseFold/FoldableCode.hpp"

namespace basefold_bench_common {

using NTL::LogicError;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEX;
using NTL::ZZ_pX;
using NTL::conv;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;

struct BasicContextSpec {
  ZZ scalar_modulus = ZZ(0);
  ZZ base_prime = ZZ(0);
  std::vector<ZZ> F_coeffs;
  std::vector<ZZ> zeta_coeffs;
};

struct Stats {
  double mean_ms = 0.0;
  double min_ms = 0.0;
  double max_ms = 0.0;
};

inline Stats ComputeStats(const std::vector<double> &xs) {
  Stats out;
  if (xs.empty()) {
    return out;
  }
  double sum = 0.0;
  out.min_ms = xs[0];
  out.max_ms = xs[0];
  for (double v : xs) {
    sum += v;
    out.min_ms = std::min(out.min_ms, v);
    out.max_ms = std::max(out.max_ms, v);
  }
  out.mean_ms = sum / static_cast<double>(xs.size());
  return out;
}

inline double MsSince(const std::chrono::steady_clock::time_point &a,
                      const std::chrono::steady_clock::time_point &b) {
  return std::chrono::duration_cast<
             std::chrono::duration<double, std::milli>>(b - a)
      .count();
}

inline long Pow2Checked(long e) {
  if (e < 0) {
    LogicError("Pow2Checked: negative exponent");
  }
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
}

inline bool IsPowerOfTwoLong(long n) { return n > 0 && (n & (n - 1)) == 0; }

inline long Log2ExactPowerOfTwoLong(long n) {
  if (!IsPowerOfTwoLong(n)) {
    LogicError("Log2ExactPowerOfTwoLong: not a power of two");
  }
  long d = 0;
  while (n > 1) {
    n >>= 1;
    ++d;
  }
  return d;
}

inline bool ParseZZString(const std::string &s, ZZ &out) {
  if (s.empty()) {
    return false;
  }
  std::size_t pos = 0;
  bool neg = false;
  if (s[pos] == '+' || s[pos] == '-') {
    neg = (s[pos] == '-');
    ++pos;
  }
  if (pos >= s.size()) {
    return false;
  }

  ZZ v(0);
  for (; pos < s.size(); ++pos) {
    const unsigned char ch = static_cast<unsigned char>(s[pos]);
    if (!std::isdigit(ch)) {
      return false;
    }
    v *= 10;
    v += static_cast<long>(ch - static_cast<unsigned char>('0'));
  }

  out = neg ? -v : v;
  return true;
}

inline bool ParseLong(const char *s, long &out) {
  if (s == nullptr) {
    return false;
  }
  try {
    std::size_t idx = 0;
    const long value = std::stol(std::string(s), &idx, 10);
    if (idx != std::string(s).size()) {
      return false;
    }
    out = value;
    return true;
  } catch (...) {
    return false;
  }
}

inline bool ParseZZ(const char *s, ZZ &out) {
  if (s == nullptr) {
    return false;
  }
  return ParseZZString(std::string(s), out);
}

inline bool ParseInt(const char *s, int &out) {
  long value = 0;
  if (!ParseLong(s, value)) {
    return false;
  }
  if (value < std::numeric_limits<int>::min() ||
      value > std::numeric_limits<int>::max()) {
    return false;
  }
  out = static_cast<int>(value);
  return true;
}

inline std::uint64_t ParseU64OrDie(const char *s, const char *flag) {
  try {
    std::size_t idx = 0;
    const std::uint64_t value = std::stoull(std::string(s), &idx, 10);
    if (idx != std::string(s).size()) {
      std::cerr << "Invalid " << flag << "\n";
      std::exit(2);
    }
    return value;
  } catch (...) {
    std::cerr << "Invalid " << flag << "\n";
    std::exit(2);
  }
}

inline std::vector<ZZ> ParseCoeffList(const std::string &s) {
  std::vector<ZZ> out;
  std::size_t pos = 0;
  while (pos < s.size()) {
    const std::size_t comma = s.find(',', pos);
    const std::size_t end = (comma == std::string::npos) ? s.size() : comma;
    std::string token = s.substr(pos, end - pos);

    const std::size_t first = token.find_first_not_of(" \t");
    const std::size_t last = token.find_last_not_of(" \t");
    if (first == std::string::npos) {
      LogicError("ParseCoeffList: empty coefficient");
    }
    token = token.substr(first, last - first + 1);

    ZZ value;
    if (!ParseZZString(token, value)) {
      LogicError("ParseCoeffList: bad integer token");
    }
    out.push_back(value);
    pos = (comma == std::string::npos) ? s.size() : (comma + 1);
  }
  if (out.empty()) {
    LogicError("ParseCoeffList: empty list");
  }
  return out;
}

inline std::vector<std::vector<ZZ>> ParseNestedCoeffList(
    const std::string &s) {
  std::vector<std::vector<ZZ>> out;
  std::size_t pos = 0;
  while (pos < s.size()) {
    const std::size_t semi = s.find(';', pos);
    const std::size_t end = (semi == std::string::npos) ? s.size() : semi;
    std::string token = s.substr(pos, end - pos);

    const std::size_t first = token.find_first_not_of(" \t");
    const std::size_t last = token.find_last_not_of(" \t");
    if (first == std::string::npos) {
      LogicError("ParseNestedCoeffList: empty coefficient block");
    }
    token = token.substr(first, last - first + 1);

    out.push_back(ParseCoeffList(token));
    pos = (semi == std::string::npos) ? s.size() : (semi + 1);
  }
  if (out.empty()) {
    LogicError("ParseNestedCoeffList: empty list");
  }
  return out;
}

inline ZZ_pX BuildZZpX(const std::vector<ZZ> &coeffs) {
  ZZ_pX poly;
  NTL::clear(poly);
  for (std::size_t i = 0; i < coeffs.size(); ++i) {
    if (coeffs[i] != 0) {
      SetCoeff(poly, static_cast<long>(i), conv<ZZ_p>(coeffs[i]));
    }
  }
  return poly;
}

inline ZZ_pE BuildZZpE(const std::vector<ZZ> &coeffs) {
  const ZZ_pX poly = BuildZZpX(coeffs);
  ZZ_pE out;
  conv(out, poly);
  return out;
}

inline ZZ_pEX BuildZZpEX(const std::vector<std::vector<ZZ>> &coeffs) {
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

inline ZZ NormalizeMod(const ZZ &x, const ZZ &mod) {
  if (mod <= 0) {
    LogicError("NormalizeMod: mod must be positive");
  }
  ZZ out = x % mod;
  if (out < 0) {
    out += mod;
  }
  return out;
}

template <typename ContextSpecLike>
inline void DeduceBasePrimeAndExponent(const ContextSpecLike &spec, ZZ &p_out,
                                       long &k_out) {
  if (spec.scalar_modulus <= 1) {
    LogicError("DeduceBasePrimeAndExponent: scalar_modulus must be > 1");
  }
  if (spec.base_prime > 1) {
    ZZ m = spec.scalar_modulus;
    long k = 0;
    while ((m % spec.base_prime) == 0) {
      m /= spec.base_prime;
      ++k;
    }
    if (k <= 0 || m != 1) {
      LogicError(
          "DeduceBasePrimeAndExponent: scalar_modulus must equal base_prime^k");
    }
    p_out = spec.base_prime;
    k_out = k;
    return;
  }
  p_out = spec.scalar_modulus;
  k_out = 1;
}

inline void ValidateMonic(const std::vector<ZZ> &coeffs, const ZZ &mod,
                          const char *what) {
  if (coeffs.empty()) {
    LogicError((std::string(what) + ": empty polynomial").c_str());
  }
  long last = static_cast<long>(coeffs.size()) - 1;
  while (last > 0 &&
         NormalizeMod(coeffs[static_cast<std::size_t>(last)], mod) == 0) {
    --last;
  }
  if (last <= 0) {
    LogicError((std::string(what) + ": degree must be >= 1").c_str());
  }
  if (NormalizeMod(coeffs[static_cast<std::size_t>(last)], mod) != 1) {
    LogicError((std::string(what) + ": leading coefficient must be 1").c_str());
  }
}

inline std::uint64_t SplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

inline ZZ ZZFromU64(std::uint64_t x) {
  static const ZZ kTwo32 = NTL::power2_ZZ(32);
  const long lo = static_cast<long>(x & 0xffffffffULL);
  const long hi = static_cast<long>((x >> 32) & 0xffffffffULL);
  ZZ out = to_ZZ(hi);
  out *= kTwo32;
  out += to_ZZ(lo);
  return out;
}

inline ZZ DeterministicResidue(std::uint64_t seed, const ZZ &modulus) {
  const long bits = NTL::NumBits(modulus);
  const long blocks = std::max<long>(1, (bits + 63) / 64);
  ZZ value(0);
  std::uint64_t x = seed;
  for (long i = 0; i < blocks; ++i) {
    x = SplitMix64(x + 0x9e3779b97f4a7c15ULL +
                   static_cast<std::uint64_t>(i));
    value <<= 64;
    value += ZZFromU64(x);
  }
  return value % modulus;
}

inline ZZ_pE BaseRingConstant(const ZZ &value) {
  ZZ_pE out;
  conv(out, BuildZZpX({value}));
  return out;
}

inline vec_ZZ_pE MakeDeterministicBaseRingTable(long count,
                                                std::uint64_t seed) {
  if (count <= 0) {
    LogicError("MakeDeterministicBaseRingTable: count must be positive");
  }
  const ZZ modulus = NTL::ZZ_p::modulus();
  vec_ZZ_pE out;
  out.SetLength(count);
  for (long i = 0; i < count; ++i) {
    const ZZ value =
        DeterministicResidue(seed ^ static_cast<std::uint64_t>(i), modulus);
    out[i] = BaseRingConstant(value);
  }
  return out;
}

inline ZZ_pE MakeDeterministicElement(std::uint64_t seed) {
  const long r = ZZ_pE::degree();
  if (r <= 0) {
    LogicError("MakeDeterministicElement: invalid extension degree");
  }
  const ZZ modulus = NTL::ZZ_p::modulus();
  if (modulus <= 1) {
    LogicError("MakeDeterministicElement: invalid modulus");
  }
  ZZ_pX poly;
  NTL::clear(poly);
  std::uint64_t x = SplitMix64(seed);
  for (long j = 0; j < r; ++j) {
    x = SplitMix64(x + static_cast<std::uint64_t>(j));
    const ZZ coeff = DeterministicResidue(
        x ^ (static_cast<std::uint64_t>(j) << 32), modulus);
    SetCoeff(poly, j, conv<ZZ_p>(coeff));
  }
  ZZ_pE out;
  conv(out, poly);
  return out;
}

inline vec_ZZ_pE MakeDeterministicCoefficients(long coeff_count,
                                               std::uint64_t seed) {
  vec_ZZ_pE coeffs;
  coeffs.SetLength(coeff_count);
  for (long i = 0; i < coeff_count; ++i) {
    coeffs[i] = MakeDeterministicElement(seed ^ static_cast<std::uint64_t>(i));
  }
  return coeffs;
}

inline std::vector<ZZ_pE> MakeDeterministicPoint(long dimension,
                                                 std::uint64_t seed) {
  if (dimension < 0) {
    LogicError("MakeDeterministicPoint: dimension must be non-negative");
  }
  std::vector<ZZ_pE> out(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    out[static_cast<std::size_t>(i)] =
        MakeDeterministicElement(seed + 0x12345678ULL +
                                 static_cast<std::uint64_t>(i) *
                                     0x9e3779b9ULL);
  }
  return out;
}

inline NTL::mat_ZZ_pE BuildSystematicG0(long c, long k0) {
  if (c <= 0) {
    LogicError("BuildSystematicG0: c must be positive");
  }
  if (k0 <= 0) {
    LogicError("BuildSystematicG0: k0 must be positive");
  }
  if (c > std::numeric_limits<long>::max() / k0) {
    LogicError("BuildSystematicG0: overflow in n0");
  }
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

inline basefold::FoldableCodeParams BuildFoldableParamsK0Eq1(
    long c, long d, const ZZ &prime_p, const ZZ_pE &zeta,
    const char *func_name) {
  if (c <= 0) {
    LogicError((std::string(func_name) + ": c must be positive").c_str());
  }
  if (d < 0) {
    LogicError((std::string(func_name) + ": d must be non-negative").c_str());
  }

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = 1;
  params.d = d;
  params.p = prime_p;
  params.zeta = zeta;

  params.G0.SetDims(1, c);
  params.G0[0][0] = ZZ_pE(1);
  for (long j = 1; j < c; ++j) {
    params.G0[0][j] = (j == 1) ? zeta : ZZ_pE(1);
  }

  params.diag_T.resize(static_cast<std::size_t>(d));
  for (long level = 0; level < d; ++level) {
    const long pow2 = Pow2Checked(level);
    if (c > std::numeric_limits<long>::max() / pow2) {
      LogicError((std::string(func_name) + ": overflow in n_i").c_str());
    }
    const long ni = c * pow2;
    params.diag_T[static_cast<std::size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      params.diag_T[static_cast<std::size_t>(level)][i] = ZZ_pE(1);
    }
  }
  return params;
}

inline basefold::FoldableCodeParams BuildFoldableParamsK0Pow2(
    long c, long k0, long d, const ZZ &prime_p, const ZZ_pE &zeta,
    const char *func_name) {
  if (c <= 0) {
    LogicError((std::string(func_name) + ": c must be positive").c_str());
  }
  if (k0 <= 0) {
    LogicError((std::string(func_name) + ": k0 must be positive").c_str());
  }
  if (!IsPowerOfTwoLong(k0)) {
    LogicError(
        (std::string(func_name) + ": k0 must be a power of two").c_str());
  }
  if (d < 0) {
    LogicError((std::string(func_name) + ": d must be non-negative").c_str());
  }
  if (c > std::numeric_limits<long>::max() / k0) {
    LogicError((std::string(func_name) + ": overflow in n0").c_str());
  }
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
      LogicError((std::string(func_name) + ": overflow in n_i").c_str());
    }
    const long ni = n0 * pow2;
    params.diag_T[static_cast<std::size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      params.diag_T[static_cast<std::size_t>(level)][i] = ZZ_pE(1);
    }
  }
  return params;
}

inline std::uint64_t FixedCoeffBytesOrThrow() {
  const ZZ mod_minus_one = NTL::ZZ_p::modulus() - 1;
  const long bits = NTL::NumBits(mod_minus_one);
  if (bits <= 0) {
    LogicError("FixedCoeffBytesOrThrow: invalid base modulus bit width");
  }
  return static_cast<std::uint64_t>((bits + 7) / 8);
}

inline std::uint64_t MulU64OrThrow(std::uint64_t lhs, std::uint64_t rhs,
                                   const char *what) {
  if (lhs == 0 || rhs == 0) {
    return 0;
  }
  if (lhs > std::numeric_limits<std::uint64_t>::max() / rhs) {
    LogicError(what);
  }
  return lhs * rhs;
}

inline std::uint64_t ExtensionElementBytesOrThrow() {
  const long degree = NTL::ZZ_pE::degree();
  if (degree <= 0) {
    LogicError("ExtensionElementBytesOrThrow: invalid extension degree");
  }
  return MulU64OrThrow(FixedCoeffBytesOrThrow(),
                       static_cast<std::uint64_t>(degree),
                       "ExtensionElementBytesOrThrow: byte width overflow");
}

inline std::uint64_t PackedVectorFixedBytesOrThrow(const vec_ZZ_pE &values) {
  return MulU64OrThrow(
      ExtensionElementBytesOrThrow(), static_cast<std::uint64_t>(values.length()),
      "PackedVectorFixedBytesOrThrow: total byte width overflow");
}

}  // namespace basefold_bench_common

#endif  // BASEFOLD_BENCH_COMMON_HELPERS_HPP_
