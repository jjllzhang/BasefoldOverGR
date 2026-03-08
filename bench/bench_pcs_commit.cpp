#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "BaseFold/BaseFoldPCS.hpp"
#include "BaseFold/Hash.hpp"
#include "GaloisRing/PrimitiveElement.hpp"

using NTL::conv;
using NTL::LogicError;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;

namespace {

struct ContextSpec {
  std::string label;
  ZZ mod = ZZ(0);      // ZZ_p modulus (p for fields, p^s for rings)
  ZZ prime_p = ZZ(0);  // optional: the prime p (only used by checked paths)
  std::vector<ZZ> F_coeffs;     // extension modulus polynomial coefficients
  std::vector<ZZ> zeta_coeffs;  // ζ element coefficients
};

struct Stats {
  double mean_ms = 0.0;
  double min_ms = 0.0;
  double max_ms = 0.0;
};

Stats ComputeStats(const std::vector<double> &xs) {
  Stats out;
  if (xs.empty())
    return out;
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

double MsSince(const std::chrono::steady_clock::time_point &a,
               const std::chrono::steady_clock::time_point &b) {
  return std::chrono::duration_cast<
             std::chrono::duration<double, std::milli>>(b - a)
      .count();
}

long Pow2Checked(long e) {
  if (e < 0) LogicError("Pow2Checked: negative exponent");
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
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
    if (first == std::string::npos) LogicError("ParseCoeffList: empty coefficient");
    token = token.substr(first, last - first + 1);

    ZZ v;
    if (!ParseZZString(token, v)) LogicError("ParseCoeffList: bad integer token");

    out.push_back(v);
    pos = (comma == std::string::npos) ? s.size() : (comma + 1);
  }
  if (out.empty()) LogicError("ParseCoeffList: empty list");
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
  while (last > 0 && NormalizeMod(coeffs[static_cast<std::size_t>(last)], mod) == 0) {
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

struct BenchResult {
  Stats encode;
  Stats top_commit;
  Stats commit;
  std::uint64_t sink = 0;
};

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

vec_ZZ_pE MakeDeterministicCoefficients(long coeff_count, std::uint64_t seed) {
  const long r = ZZ_pE::degree();
  if (r <= 0) LogicError("MakeDeterministicCoefficients: invalid extension degree");
  const ZZ modulus = NTL::ZZ_p::modulus();
  if (modulus <= 1) LogicError("MakeDeterministicCoefficients: invalid modulus");

  vec_ZZ_pE coeffs;
  coeffs.SetLength(coeff_count);
  for (long i = 0; i < coeff_count; ++i) {
    ZZ_pX poly;
    NTL::clear(poly);
    std::uint64_t x = SplitMix64(seed ^ static_cast<std::uint64_t>(i));
    for (long j = 0; j < r; ++j) {
      x = SplitMix64(x + static_cast<std::uint64_t>(j));
      const ZZ cj = DeterministicResidue(
          x ^ (static_cast<std::uint64_t>(j) << 32), modulus);
      SetCoeff(poly, j, conv<ZZ_p>(cj));
    }
    ZZ_pE elem;
    conv(elem, poly);
    coeffs[i] = elem;
  }
  return coeffs;
}

NTL::mat_ZZ_pE BuildVandermondeG0(long k0, const vec_ZZ_pE &points) {
  if (k0 <= 0) LogicError("BuildVandermondeG0: k0 must be positive");
  const long n0 = points.length();
  if (n0 <= 0) LogicError("BuildVandermondeG0: points must be non-empty");

  NTL::mat_ZZ_pE G0;
  G0.SetDims(k0, n0);

  for (long j = 0; j < n0; ++j) {
    ZZ_pE cur;
    NTL::set(cur);  // 1
    for (long r = 0; r < k0; ++r) {
      G0[r][j] = cur;
      cur *= points[j];
    }
  }

  return G0;
}

basefold::FoldableCodeParams BuildParams_k0_gt1(long c, long k0, long d,
                                               const ZZ &prime_p,
                                               const ZZ_pE &zeta,
                                               std::uint64_t seed) {
  if (c <= 0) LogicError("BuildParams_k0_gt1: c must be positive");
  if (k0 <= 1) LogicError("BuildParams_k0_gt1: k0 must be > 1");
  if (d < 0) LogicError("BuildParams_k0_gt1: d must be non-negative");

  if (c > std::numeric_limits<long>::max() / k0) {
    LogicError("BuildParams_k0_gt1: overflow in n0");
  }
  const long n0 = c * k0;

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = prime_p;
  params.zeta = zeta;

  const vec_ZZ_pE points =
      MakeDeterministicCoefficients(n0, seed ^ 0x7f4a7c15ULL);
  params.G0 = BuildVandermondeG0(k0, points);

  params.diag_T.resize(static_cast<std::size_t>(d));
  for (long level = 0; level < d; ++level) {
    const long pow2 = Pow2Checked(level);
    if (n0 > std::numeric_limits<long>::max() / pow2) {
      LogicError("BuildParams_k0_gt1: overflow in n_i");
    }
    const long ni = n0 * pow2;
    params.diag_T[static_cast<std::size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      params.diag_T[static_cast<std::size_t>(level)][i] = ZZ_pE(1);
    }
  }

  return params;
}

BenchResult RunCommitBenchmark(const vec_ZZ_pE &f_coeffs,
                               const basefold::FoldableCodeParams &params,
                               int warmup, int reps) {
  if (warmup < 0) LogicError("RunCommitBenchmark: warmup must be >= 0");
  if (reps <= 0) LogicError("RunCommitBenchmark: reps must be > 0");

  std::vector<double> encode_ms;
  std::vector<double> top_commit_ms;
  std::vector<double> commit_ms;
  encode_ms.reserve(static_cast<std::size_t>(reps));
  top_commit_ms.reserve(static_cast<std::size_t>(reps));
  commit_ms.reserve(static_cast<std::size_t>(reps));

  std::uint64_t sink = 0;
  vec_ZZ_pE pi_d;
  basefold::MerkleTree merkle_d;

  for (int iter = -warmup; iter < reps; ++iter) {
    const auto t0 = std::chrono::steady_clock::now();
    basefold::EncodeFoldableUnchecked(pi_d, f_coeffs, params);
    const auto t1 = std::chrono::steady_clock::now();
    merkle_d = basefold::MerkleTree::Build(pi_d);
    const basefold::MerkleRoot root_d = merkle_d.Root();
    const auto t2 = std::chrono::steady_clock::now();

    if (pi_d.length() > 0) {
      // Touch a single output element to discourage aggressive dead-code
      // elimination without adding heavy conversions that slow the benchmark.
      sink ^= static_cast<std::uint64_t>(pi_d[0] != ZZ_pE(0));
    }
    if (!root_d.empty()) {
      sink ^= static_cast<std::uint64_t>(root_d[0]);
    }

    if (iter >= 0) {
      encode_ms.push_back(MsSince(t0, t1));
      top_commit_ms.push_back(MsSince(t1, t2));
      commit_ms.push_back(MsSince(t0, t2));
    }
  }

  BenchResult out;
  out.encode = ComputeStats(encode_ms);
  out.top_commit = ComputeStats(top_commit_ms);
  out.commit = ComputeStats(commit_ms);
  out.sink = sink;
  return out;
}

void PrintResult(const std::string &label, const ZZ &mod, long c, long k0, long d,
                 long k_d, long n_d, int warmup, int reps,
                 const BenchResult &r) {
  std::cout << "\n[" << label << "] c=" << c << " k0=" << k0 << " d=" << d
            << "  mod=" << mod << "  k_d=" << k_d << "  n_d=" << n_d
            << "  warmup=" << warmup << " reps=" << reps << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  hash backend " << basefold::SelectedHashBackendName() << "\n";
  std::cout << "  encode-only mean " << r.encode.mean_ms << " ms  (min "
            << r.encode.min_ms << ", max " << r.encode.max_ms << ")\n";
  std::cout << "  top-commit mean " << r.top_commit.mean_ms << " ms  (min "
            << r.top_commit.min_ms << ", max " << r.top_commit.max_ms << ")\n";
  std::cout << "  commit     mean " << r.commit.mean_ms << " ms  (min "
            << r.commit.min_ms << ", max " << r.commit.max_ms << ")\n";
  std::cout << "  sink    " << r.sink << "\n";
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

bool ParseInt(const char *s, int &out) {
  long v = 0;
  if (!ParseLong(s, v)) return false;
  if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max())
    return false;
  out = static_cast<int>(v);
  return true;
}

void PrintHelp() {
  std::cout
      << "bench_pcs_commit (top commit benchmark, unchecked)\n\n"
      << "Usage:\n"
      << "  bench_pcs_commit [--mode field|ring|both] [--c <int>] [--k0 <int>] [--d <int>]\n"
      << "                 [--warmup <int>] [--reps <int>] [--seed <u64>]\n"
      << "                 [--auto-zeta teich]\n"
      << "                 [--field-mod <decimal-int>] [--field-F <a0,a1,...>] [--field-zeta <b0,b1,...>]\n"
      << "                 [--ring-mod <decimal-int>]  [--ring-p <decimal-int>] [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n\n"
      << "Notes:\n"
      << "  Headline commit time includes top-level EncodeFoldable + MerkleTree::Build.\n"
      << "  Auxiliary lines print encode-only and top-commit-only splits.\n\n"
      << "  With --auto-zeta teich, zeta is derived as a Teichmuller generator from (p,k,F);\n"
      << "  then --field-zeta/--ring-zeta are ignored.\n\n"
      << "Examples:\n"
      << "  # GF(2^2) with F(x)=x^2+x+1 and zeta=x\n"
      << "  bench_pcs_commit --mode field --field-mod 2 --field-F 1,1,1 --field-zeta 0,1 --d 16\n"
      << "  # GR(4,2) with the same extension polynomial and zeta=x\n"
      << "  bench_pcs_commit --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 16\n"
      << "  # k0>1 (exercise the general-k0 encoder path)\n"
      << "  bench_pcs_commit --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --k0 2 --d 16\n"
      << "  # Auto zeta from Teichmuller subgroup generator\n"
      << "  bench_pcs_commit --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --auto-zeta teich --d 16\n";
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

void RunOneContext(const ContextSpec &spec, long c, long k0, long d, int warmup,
                   int reps, bool auto_zeta_teich, std::uint64_t seed) {
  if (spec.mod <= 1) LogicError("RunOneContext: modulus must be > 1");

  const ZZ modulus = spec.mod;
  ZZ_pPush mod_push(modulus);

  ValidateMonic(spec.F_coeffs, spec.mod, "F");
  const ZZ_pX F = BuildZZpX(spec.F_coeffs);
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

  const basefold::FoldableCodeParams params = [&] {
    if (k0 == 1) {
      return BuildParams_k0_1(
          c, d, (spec.prime_p > 1) ? spec.prime_p : spec.mod, zeta);
    }
    return BuildParams_k0_gt1(
        c, k0, d, (spec.prime_p > 1) ? spec.prime_p : spec.mod, zeta, seed);
  }();

  const long pow2_d = Pow2Checked(d);
  if (k0 > std::numeric_limits<long>::max() / pow2_d) {
    LogicError("RunOneContext: overflow in k_d");
  }
  const long k_d = k0 * pow2_d;
  if (c > std::numeric_limits<long>::max() / k_d) {
    LogicError("RunOneContext: overflow in n_d");
  }
  const long n_d = c * k_d;

  const vec_ZZ_pE f_coeffs = MakeDeterministicCoefficients(k_d, seed);
  const BenchResult r = RunCommitBenchmark(f_coeffs, params, warmup, reps);
  PrintResult(spec.label, spec.mod, c, k0, d, k_d, n_d, warmup, reps, r);
}

}  // namespace

int main(int argc, char **argv) {
  long d = 16;
  long c = 2;
  long k0 = 1;
  int warmup = 1;
  int reps = 5;
  std::uint64_t seed = 0;
  bool auto_zeta_teich = false;

  bool do_field = true;
  bool do_ring = true;

  ContextSpec field;
  field.label = "Field";
  field.mod = to_ZZ(2);
  field.prime_p = ZZ(0);
  field.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};      // x^2 + x + 1
  field.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};      // x

  ContextSpec ring;
  ring.label = "Ring";
  ring.mod = to_ZZ(4);
  ring.prime_p = to_ZZ(2);
  ring.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};       // x^2 + x + 1
  ring.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};       // x

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
    } else if (arg == "--warmup") {
      if (!ParseInt(NeedValue("--warmup"), warmup) || warmup < 0) {
        std::cerr << "Invalid --warmup\n";
        return 2;
      }
    } else if (arg == "--reps") {
      if (!ParseInt(NeedValue("--reps"), reps) || reps <= 0) {
        std::cerr << "Invalid --reps\n";
        return 2;
      }
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
    } else if (arg == "--help" || arg == "-h") {
      PrintHelp();
      return 0;
    } else {
      std::cerr << "Unknown arg: " << arg << "\n";
      return 2;
    }
  }

  try {
    if (!do_field && !do_ring) {
      std::cerr << "Nothing to do: --mode disabled both field and ring\n";
      return 2;
    }
    if (do_field)
      RunOneContext(field, c, k0, d, warmup, reps, auto_zeta_teich, seed);
    if (do_ring)
      RunOneContext(ring, c, k0, d, warmup, reps, auto_zeta_teich, seed);
  } catch (const std::exception &e) {
    std::cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    std::cerr << "Unhandled non-std exception\n";
    return 2;
  }

  return 0;
}
