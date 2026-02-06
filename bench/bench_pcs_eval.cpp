#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
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
#include "BaseFold/Multilinear.hpp"
#include "BaseFold/Profile.hpp"

using NTL::conv;
using NTL::LogicError;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;

namespace {

struct ContextSpec {
  std::string label;
  long mod = 0;      // ZZ_p modulus (p for fields, p^s for rings)
  long prime_p = 0;  // optional: the prime p (only used by unit checks)
  std::vector<long> F_coeffs;     // extension modulus polynomial coefficients
  std::vector<long> zeta_coeffs;  // ζ element coefficients
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

std::vector<long> ParseCoeffList(const std::string &s) {
  std::vector<long> out;
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

    std::size_t idx = 0;
    long v = 0;
    try {
      v = std::stol(token, &idx, 10);
    } catch (...) {
      LogicError("ParseCoeffList: bad integer token");
    }
    if (idx != token.size())
      LogicError("ParseCoeffList: bad integer token");

    out.push_back(v);
    pos = (comma == std::string::npos) ? s.size() : (comma + 1);
  }
  if (out.empty())
    LogicError("ParseCoeffList: empty list");
  return out;
}

ZZ_pX BuildZZpX(const std::vector<long> &coeffs) {
  ZZ_pX poly;
  NTL::clear(poly);
  for (std::size_t i = 0; i < coeffs.size(); ++i) {
    if (coeffs[i] != 0) {
      SetCoeff(poly, static_cast<long>(i), coeffs[i]);
    }
  }
  return poly;
}

long NormalizeMod(long x, long mod) {
  if (mod <= 0) LogicError("NormalizeMod: mod must be positive");
  long r = x % mod;
  if (r < 0) r += mod;
  return r;
}

void ValidateMonic(const std::vector<long> &coeffs, long mod,
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
  const long lead = NormalizeMod(coeffs[static_cast<std::size_t>(last)], mod);
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

std::uint64_t SplitMix64(std::uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

long BaseModulusAsLong() {
  const ZZ modulus_zz = NTL::ZZ_p::modulus();
  long modulus = 0;
  NTL::conv(modulus, modulus_zz);
  if (modulus <= 1)
    LogicError("BaseModulusAsLong: invalid modulus");
  return modulus;
}

ZZ_pE MakeDeterministicElement(std::uint64_t seed) {
  const long r = ZZ_pE::degree();
  if (r <= 0) LogicError("MakeDeterministicElement: invalid extension degree");
  const long modulus = BaseModulusAsLong();

  ZZ_pX poly;
  NTL::clear(poly);
  std::uint64_t x = SplitMix64(seed);
  for (long j = 0; j < r; ++j) {
    x = SplitMix64(x + static_cast<std::uint64_t>(j));
    const long cj = static_cast<long>(x % static_cast<std::uint64_t>(modulus));
    SetCoeff(poly, j, cj);
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

struct BenchResult {
  Stats prover;
  Stats verifier;
  std::uint64_t sink = 0;
  basefold::Profile prover_profile;
  basefold::Profile verifier_profile;
  bool has_profile = false;
};

BenchResult RunEvalBenchmark(const vec_ZZ_pE &f_coeffs,
                             const std::vector<ZZ_pE> &z,
                             const ZZ_pE &y,
                             long num_queries,
                             const basefold::FoldableCodeParams &params,
                             bool checked_prover,
                             bool enable_profile,
                             int warmup, int reps) {
  if (warmup < 0) LogicError("RunEvalBenchmark: warmup must be >= 0");
  if (reps <= 0) LogicError("RunEvalBenchmark: reps must be > 0");
  if (num_queries < 0) LogicError("RunEvalBenchmark: num_queries must be >= 0");

  std::vector<double> prover_ms;
  std::vector<double> verifier_ms;
  prover_ms.reserve(static_cast<std::size_t>(reps));
  verifier_ms.reserve(static_cast<std::size_t>(reps));

  basefold::Profile prover_prof;
  basefold::Profile verifier_prof;
  basefold::ResetProfile(prover_prof);
  basefold::ResetProfile(verifier_prof);

  std::uint64_t sink = 0;

  for (int iter = -warmup; iter < reps; ++iter) {
    const auto t0 = std::chrono::steady_clock::now();
    const basefold::BaseFoldPCSEvalProof proof = [&] {
      if (enable_profile && iter >= 0) {
        basefold::ProfileGuard guard(&prover_prof);
        basefold::ScopedTimer timer(&prover_prof.pcs_prove_ns,
                                    &prover_prof.pcs_prove_calls);
        return checked_prover
                   ? basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries,
                                                   params)
                   : basefold::BaseFoldPCSProveEvalUnchecked(
                         f_coeffs, z, y, num_queries, params);
      }
      return checked_prover
                 ? basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries,
                                                 params)
                 : basefold::BaseFoldPCSProveEvalUnchecked(f_coeffs, z, y,
                                                           num_queries, params);
    }();
    const auto t1 = std::chrono::steady_clock::now();

    const basefold::MerkleRoot &C =
        proof.commitments.roots_by_level[static_cast<std::size_t>(params.d)];

    const auto t2 = std::chrono::steady_clock::now();
    bool ok = false;
    if (enable_profile && iter >= 0) {
      basefold::ProfileGuard guard(&verifier_prof);
      ok = basefold::BaseFoldPCSVerifyEval(C, z, y, num_queries, proof, params);
    } else {
      ok = basefold::BaseFoldPCSVerifyEval(C, z, y, num_queries, proof, params);
    }
    const auto t3 = std::chrono::steady_clock::now();

    if (!ok) {
      LogicError("RunEvalBenchmark: verification failed");
    }

    // Prevent over-optimization: fold in a few bytes.
    if (!C.empty())
      sink ^= static_cast<std::uint64_t>(C[0]);
    sink ^= static_cast<std::uint64_t>(ok);

    if (iter >= 0) {
      prover_ms.push_back(MsSince(t0, t1));
      verifier_ms.push_back(MsSince(t2, t3));
    }
  }

  BenchResult out;
  out.prover = ComputeStats(prover_ms);
  out.verifier = ComputeStats(verifier_ms);
  out.sink = sink;
  out.prover_profile = prover_prof;
  out.verifier_profile = verifier_prof;
  out.has_profile = enable_profile;
  return out;
}

void PrintResult(const std::string &label, long mod, long c, long d,
                 long k0,
                 long num_queries, int warmup, int reps,
                 const BenchResult &r) {
  const long pow2_d = Pow2Checked(d);
  if (k0 > std::numeric_limits<long>::max() / pow2_d) {
    LogicError("PrintResult: overflow in k_d");
  }
  const long k_d = k0 * pow2_d;
  if (c > std::numeric_limits<long>::max() / k_d) {
    LogicError("PrintResult: overflow in n_d");
  }
  const long n_d = c * k_d;

  std::cout << "\n[" << label << "] c=" << c << " k0=" << k0 << " d=" << d
            << "  mod=" << mod << "  k_d=" << k_d << "  n_d=" << n_d
            << "  queries=" << num_queries << "  warmup=" << warmup
            << " reps=" << reps << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  prover   mean " << r.prover.mean_ms << " ms  (min "
            << r.prover.min_ms << ", max " << r.prover.max_ms << ")\n";
  std::cout << "  verifier mean " << r.verifier.mean_ms << " ms  (min "
            << r.verifier.min_ms << ", max " << r.verifier.max_ms << ")\n";

  std::cout << "  sink    " << r.sink << "\n";
  if (r.has_profile) {
    basefold::PrintProfile(std::cout, r.prover_profile);
    basefold::PrintProfile(std::cout, r.verifier_profile);
  }
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

bool ParseInt(const char *s, int &out) {
  long v = 0;
  if (!ParseLong(s, v)) return false;
  if (v < std::numeric_limits<int>::min() || v > std::numeric_limits<int>::max())
    return false;
  out = static_cast<int>(v);
  return true;
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
      << "bench_pcs_eval (PCS prove+verify, includes Merkle+FS)\n\n"
      << "Usage:\n"
      << "  bench_pcs_eval [--mode field|ring|both] [--c <int>] [--k0 <int>] [--d <int>]\n"
      << "               [--queries <int>] [--checked] [--profile] [--warmup <int>] [--reps <int>] [--seed <u64>]\n"
      << "               [--field-mod <int>] [--field-F <a0,a1,...>] [--field-zeta <b0,b1,...>]\n"
      << "               [--ring-mod <int>]  [--ring-p <int>] [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n\n"
      << "Notes:\n"
      << "  By default, prover uses BaseFoldPCSProveEvalUnchecked (skips validation and claimed_y check).\n\n"
      << "  PCS Eval supports k0 = 2^κ. The multilinear point dimension is (d + κ).\n\n"
      << "Examples:\n"
      << "  # GF(2^2) with F(x)=x^2+x+1 and zeta=x\n"
      << "  bench_pcs_eval --mode field --field-mod 2 --field-F 1,1,1 --field-zeta 0,1 --d 16 --queries 4\n"
      << "  # GR(4,2) with the same extension polynomial and zeta=x\n"
      << "  bench_pcs_eval --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 16 --queries 4\n";
}

void RunOneContext(const ContextSpec &spec, long c, long k0, long d,
                   long num_queries,
                   bool checked_prover, bool enable_profile, int warmup,
                   int reps,
                   std::uint64_t seed) {
  if (spec.mod <= 1) LogicError("RunOneContext: modulus must be > 1");
  if (k0 <= 0) LogicError("RunOneContext: k0 must be positive");
  if (!IsPowerOfTwoLong(k0)) LogicError("RunOneContext: k0 must be a power of two");

  const ZZ modulus = to_ZZ(spec.mod);
  ZZ_pPush mod_push(modulus);

  ValidateMonic(spec.F_coeffs, spec.mod, "F");
  const ZZ_pX F = BuildZZpX(spec.F_coeffs);
  ZZ_pEPush e_push(F);

  const ZZ_pX zpoly = BuildZZpX(spec.zeta_coeffs);
  ZZ_pE zeta;
  conv(zeta, zpoly);

  const basefold::FoldableCodeParams params = (k0 == 1)
                                                  ? BuildParams_k0_1(
                                                        c, d, to_ZZ(spec.prime_p), zeta)
                                                  : BuildParams_k0_pow2(
                                                        c, k0, d, to_ZZ(spec.prime_p),
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

  const BenchResult r =
      RunEvalBenchmark(f_coeffs, z, y, num_queries, params, checked_prover,
                       enable_profile, warmup, reps);
  PrintResult(spec.label, spec.mod, c, d, k0, num_queries, warmup, reps, r);
}

}  // namespace

int main(int argc, char **argv) {
  long d = 16;
  long c = 2;
  long k0 = 1;
  long num_queries = 4;
  bool checked_prover = false;
  bool enable_profile = false;
  int warmup = 1;
  int reps = 3;
  std::uint64_t seed = 0;

  bool do_field = true;
  bool do_ring = true;

  ContextSpec field;
  field.label = "Field";
  field.mod = 2;
  field.prime_p = 0;
  field.F_coeffs = {1, 1, 1};   // x^2 + x + 1
  field.zeta_coeffs = {0, 1};   // x

  ContextSpec ring;
  ring.label = "Ring";
  ring.mod = 4;
  ring.prime_p = 2;
  ring.F_coeffs = {1, 1, 1};    // x^2 + x + 1
  ring.zeta_coeffs = {0, 1};    // x

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
    } else if (arg == "--checked") {
      checked_prover = true;
    } else if (arg == "--profile") {
      enable_profile = true;
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
    } else if (arg == "--field-mod") {
      if (!ParseLong(NeedValue("--field-mod"), field.mod) || field.mod <= 1) {
        std::cerr << "Invalid --field-mod\n";
        return 2;
      }
    } else if (arg == "--field-F") {
      field.F_coeffs = ParseCoeffList(NeedValue("--field-F"));
    } else if (arg == "--field-zeta") {
      field.zeta_coeffs = ParseCoeffList(NeedValue("--field-zeta"));
    } else if (arg == "--ring-mod") {
      if (!ParseLong(NeedValue("--ring-mod"), ring.mod) || ring.mod <= 1) {
        std::cerr << "Invalid --ring-mod\n";
        return 2;
      }
    } else if (arg == "--ring-p") {
      if (!ParseLong(NeedValue("--ring-p"), ring.prime_p) || ring.prime_p <= 1) {
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
      RunOneContext(field, c, k0, d, num_queries, checked_prover, enable_profile,
                    warmup, reps, seed);
    if (do_ring)
      RunOneContext(ring, c, k0, d, num_queries, checked_prover, enable_profile,
                    warmup, reps, seed);
  } catch (const std::exception &e) {
    std::cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    std::cerr << "Unhandled non-std exception\n";
    return 2;
  }

  return 0;
}
