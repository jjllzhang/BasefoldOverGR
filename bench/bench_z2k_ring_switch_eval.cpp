#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "PCS/Common/Hash.hpp"
#include "PCS/Common/Multilinear.hpp"
#include "PCS/Common/Profile.hpp"
#include "Compiler/Z2k/BaseFoldBackendAdapter.hpp"
#include "Compiler/Z2k/RingSwitchPCS.hpp"
#include "Compiler/Z2k/RingSwitchProofSerialize.hpp"
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

struct BenchResult {
  Stats prove_total;
  Stats prove_outer;
  Stats prove_backend;
  Stats verify_total;
  Stats verify_outer;
  Stats verify_backend;
  std::uint64_t outer_proof_size_bytes = 0;
  std::uint64_t proof_size_bytes = 0;
  double outer_proof_size_kb = 0.0;
  double proof_size_kb = 0.0;
  std::uint64_t anti_opt_checksum = 0;
};

Stats ComputeStats(const std::vector<double> &xs) {
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

double MsSince(const std::chrono::steady_clock::time_point &a,
               const std::chrono::steady_clock::time_point &b) {
  return std::chrono::duration_cast<
             std::chrono::duration<double, std::milli>>(b - a)
      .count();
}

long Pow2Checked(long e) {
  if (e < 0) {
    LogicError("Pow2Checked: negative exponent");
  }
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
}

bool ParseZZString(const std::string &s, ZZ &out) {
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

bool ParseLong(const char *s, long &out) {
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

bool ParseZZ(const char *s, ZZ &out) {
  if (s == nullptr) {
    return false;
  }
  return ParseZZString(std::string(s), out);
}

bool ParseInt(const char *s, int &out) {
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

std::uint64_t ParseU64OrDie(const char *s, const char *flag) {
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

std::vector<ZZ> ParseCoeffList(const std::string &s) {
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

ZZ NormalizeMod(const ZZ &x, const ZZ &mod) {
  if (mod <= 0) {
    LogicError("NormalizeMod: mod must be positive");
  }
  ZZ out = x % mod;
  if (out < 0) {
    out += mod;
  }
  return out;
}

void ValidateMonic(const std::vector<ZZ> &coeffs, const ZZ &mod,
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

ZZ_pE BaseRingConstant(const ZZ &value) {
  ZZ_pE out;
  conv(out, BuildZZpX({value}));
  return out;
}

vec_ZZ_pE MakeDeterministicBaseRingTable(long count, std::uint64_t seed) {
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

ZZ_pE MakeDeterministicElement(std::uint64_t seed) {
  const long r = ZZ_pE::degree();
  if (r <= 0) {
    LogicError("MakeDeterministicElement: invalid extension degree");
  }
  const ZZ modulus = NTL::ZZ_p::modulus();
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

std::vector<ZZ_pE> MakeDeterministicPoint(long dimension, std::uint64_t seed) {
  if (dimension < 0) {
    LogicError("MakeDeterministicPoint: dimension must be non-negative");
  }
  std::vector<ZZ_pE> out(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    out[static_cast<std::size_t>(i)] =
        MakeDeterministicElement(seed + 0x12345678ULL +
                                 static_cast<std::uint64_t>(i) * 0x9e3779b9ULL);
  }
  return out;
}

void DeduceBasePrimeAndExponent(const ContextSpec &spec, ZZ &p_out, long &k_out) {
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

basefold::FoldableCodeParams BuildBackendParams(long c, long d, const ZZ &p,
                                                const ZZ_pE &zeta) {
  if (c <= 0) {
    LogicError("BuildBackendParams: c must be positive");
  }
  if (d < 0) {
    LogicError("BuildBackendParams: d must be non-negative");
  }

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = 1;
  params.d = d;
  params.p = p;
  params.zeta = zeta;

  params.G0.SetDims(1, c);
  params.G0[0][0] = ZZ_pE(1);
  for (long j = 1; j < c; ++j) {
    params.G0[0][j] = (j == 1) ? zeta : ZZ_pE(1);
  }

  params.diag_T.resize(static_cast<std::size_t>(d));
  for (long level = 0; level < d; ++level) {
    const long ni = c * Pow2Checked(level);
    params.diag_T[static_cast<std::size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      params.diag_T[static_cast<std::size_t>(level)][i] = ZZ_pE(1);
    }
  }
  return params;
}

basefold::RingSwitchPCSParams BuildRingSwitchParams(
    long c, long ell, long kappa, const ContextSpec &spec, const ZZ_pX &F,
    const ZZ_pE &zeta) {
  ZZ p_base;
  long k_base = 0;
  DeduceBasePrimeAndExponent(spec, p_base, k_base);
  (void)k_base;

  basefold::RingSwitchPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = spec.scalar_modulus;
  input.extension_modulus = F;
  input.alpha_basis = basefold::ActivePolynomialBasisDescriptor();
  input.beta_basis = basefold::ActivePolynomialBasisDescriptor();
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(
      BuildBackendParams(c, ell - kappa, p_base, zeta));
  return basefold::RingSwitchPCSSetup(input);
}

BenchResult RunEvalBenchmark(const basefold::RingSwitchPCSParams &params,
                             const vec_ZZ_pE &t_table,
                             const std::vector<ZZ_pE> &z,
                             const ZZ_pE &claimed_s, long num_queries,
                             int warmup, int reps) {
  if (warmup < 0) {
    LogicError("RunEvalBenchmark: warmup must be >= 0");
  }
  if (reps <= 0) {
    LogicError("RunEvalBenchmark: reps must be > 0");
  }
  if (num_queries < 0) {
    LogicError("RunEvalBenchmark: num_queries must be >= 0");
  }

  std::vector<double> prove_total_ms;
  std::vector<double> prove_outer_ms;
  std::vector<double> prove_backend_ms;
  std::vector<double> verify_total_ms;
  std::vector<double> verify_outer_ms;
  std::vector<double> verify_backend_ms;
  prove_total_ms.reserve(static_cast<std::size_t>(reps));
  prove_outer_ms.reserve(static_cast<std::size_t>(reps));
  prove_backend_ms.reserve(static_cast<std::size_t>(reps));
  verify_total_ms.reserve(static_cast<std::size_t>(reps));
  verify_outer_ms.reserve(static_cast<std::size_t>(reps));
  verify_backend_ms.reserve(static_cast<std::size_t>(reps));

  std::uint64_t anti_opt_checksum = 0;
  std::uint64_t outer_proof_size_bytes = 0;
  std::uint64_t proof_size_bytes = 0;

  for (int iter = -warmup; iter < reps; ++iter) {
    const basefold::RingSwitchPCSCommitArtifacts commit_artifacts =
        basefold::RingSwitchPCSBuildCommitArtifacts(params, t_table);

    basefold::Profile prover_prof;
    basefold::Profile verifier_prof;
    basefold::ResetProfile(prover_prof);
    basefold::ResetProfile(verifier_prof);

    const auto t0 = std::chrono::steady_clock::now();
    basefold::RingSwitchPCSEvalProof proof;
    {
      basefold::ProfileGuard guard(&prover_prof);
      proof = basefold::RingSwitchPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, num_queries, commit_artifacts);
    }
    const auto t1 = std::chrono::steady_clock::now();

    outer_proof_size_bytes =
        basefold::RingSwitchPCSOuterProofSizeBytes(params, proof);
    proof_size_bytes = basefold::RingSwitchPCSEvalProofSizeBytes(params, proof);

    const auto t2 = std::chrono::steady_clock::now();
    bool ok = false;
    {
      basefold::ProfileGuard guard(&verifier_prof);
      ok = basefold::RingSwitchPCSVerifyEval(params, commit_artifacts.commitment,
                                             z, claimed_s, num_queries, proof);
    }
    const auto t3 = std::chrono::steady_clock::now();
    if (!ok) {
      LogicError("RunEvalBenchmark: verification failed");
    }

    const double prove_total = MsSince(t0, t1);
    const double prove_backend =
        basefold::NsToMs(prover_prof.z2k_backend_prove_ns);
    const double verify_total = MsSince(t2, t3);
    const double verify_backend =
        basefold::NsToMs(verifier_prof.z2k_backend_verify_ns);
    const double prove_outer = std::max(0.0, prove_total - prove_backend);
    const double verify_outer = std::max(0.0, verify_total - verify_backend);

    anti_opt_checksum ^= static_cast<std::uint64_t>(ok);
    anti_opt_checksum ^= proof_size_bytes;
    if (!commit_artifacts.commitment.empty()) {
      anti_opt_checksum ^=
          static_cast<std::uint64_t>(commit_artifacts.commitment[0]);
    }

    if (iter >= 0) {
      prove_total_ms.push_back(prove_total);
      prove_outer_ms.push_back(prove_outer);
      prove_backend_ms.push_back(prove_backend);
      verify_total_ms.push_back(verify_total);
      verify_outer_ms.push_back(verify_outer);
      verify_backend_ms.push_back(verify_backend);
    }
  }

  BenchResult out;
  out.prove_total = ComputeStats(prove_total_ms);
  out.prove_outer = ComputeStats(prove_outer_ms);
  out.prove_backend = ComputeStats(prove_backend_ms);
  out.verify_total = ComputeStats(verify_total_ms);
  out.verify_outer = ComputeStats(verify_outer_ms);
  out.verify_backend = ComputeStats(verify_backend_ms);
  out.outer_proof_size_bytes = outer_proof_size_bytes;
  out.proof_size_bytes = proof_size_bytes;
  out.outer_proof_size_kb =
      static_cast<double>(outer_proof_size_bytes) / 1024.0;
  out.proof_size_kb = static_cast<double>(proof_size_bytes) / 1024.0;
  out.anti_opt_checksum = anti_opt_checksum;
  return out;
}

void PrintResult(long c, long ell, long kappa, long queries, int warmup,
                 int reps, const BenchResult &result) {
  std::cout << "\n[ring-switch eval]"
            << " c=" << c << " ell=" << ell << " kappa=" << kappa
            << " ell'=" << (ell - kappa) << " queries=" << queries
            << " warmup=" << warmup << " reps=" << reps << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  hash backend " << basefold::SelectedHashBackendName() << "\n";
  std::cout << "  prove-phase mean " << result.prove_total.mean_ms << " ms  (min "
            << result.prove_total.min_ms << ", max "
            << result.prove_total.max_ms << ")\n";
  std::cout << "  outer prover mean " << result.prove_outer.mean_ms << " ms  (min "
            << result.prove_outer.min_ms << ", max "
            << result.prove_outer.max_ms << ")\n";
  std::cout << "  backend prover mean " << result.prove_backend.mean_ms
            << " ms  (min " << result.prove_backend.min_ms << ", max "
            << result.prove_backend.max_ms << ")\n";
  std::cout << "  verifier mean " << result.verify_total.mean_ms << " ms  (min "
            << result.verify_total.min_ms << ", max "
            << result.verify_total.max_ms << ")\n";
  std::cout << "  outer verifier mean " << result.verify_outer.mean_ms
            << " ms  (min " << result.verify_outer.min_ms << ", max "
            << result.verify_outer.max_ms << ")\n";
  std::cout << "  backend verifier mean " << result.verify_backend.mean_ms
            << " ms  (min " << result.verify_backend.min_ms << ", max "
            << result.verify_backend.max_ms << ")\n";
  std::cout << "  outer proof size " << result.outer_proof_size_kb << " KB  ("
            << result.outer_proof_size_bytes << " B)\n";
  std::cout << "  proof size  " << result.proof_size_kb << " KB  ("
            << result.proof_size_bytes << " B)\n";
  std::cout << "  anti-opt checksum " << result.anti_opt_checksum << "\n";
}

void PrintHelp() {
  std::cout
      << "bench_z2k_ring_switch_eval\n\n"
      << "Usage:\n"
      << "  bench_z2k_ring_switch_eval [--c <int>] [--ell <int>] [--kappa <int>]\n"
      << "                             [--queries <int>] [--warmup <int>] [--reps <int>]\n"
      << "                             [--seed <u64>] [--auto-zeta teich]\n"
      << "                             [--ring-mod <decimal-int>] [--ring-p <decimal-int>]\n"
      << "                             [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n\n"
      << "Notes:\n"
      << "  Commit/artifact construction happens before timed prove, matching current bench_pcs_eval semantics.\n"
      << "  Outer prover/verifier times are total times with the backend prove/verify subcall removed.\n"
      << "  Proof size reports exact serializer-backed bytes for the public RingSwitchPCSEvalProof.\n";
}

}  // namespace

int main(int argc, char **argv) {
  long c = 2;
  long ell = 3;
  long kappa = 1;
  long queries = 2;
  int warmup = 1;
  int reps = 3;
  std::uint64_t seed = 0x5eed5678ULL;
  bool auto_zeta_teich = false;

  ContextSpec spec;
  spec.scalar_modulus = to_ZZ(4);
  spec.base_prime = to_ZZ(2);
  spec.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};
  spec.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    const auto need_value = [&](const char *flag) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << flag << "\n";
        std::exit(2);
      }
      return argv[++i];
    };

    if (arg == "--help" || arg == "-h") {
      PrintHelp();
      return 0;
    } else if (arg == "--c") {
      if (!ParseLong(need_value("--c"), c)) {
        std::cerr << "Invalid --c\n";
        return 2;
      }
    } else if (arg == "--ell") {
      if (!ParseLong(need_value("--ell"), ell)) {
        std::cerr << "Invalid --ell\n";
        return 2;
      }
    } else if (arg == "--kappa") {
      if (!ParseLong(need_value("--kappa"), kappa)) {
        std::cerr << "Invalid --kappa\n";
        return 2;
      }
    } else if (arg == "--queries") {
      if (!ParseLong(need_value("--queries"), queries)) {
        std::cerr << "Invalid --queries\n";
        return 2;
      }
    } else if (arg == "--warmup") {
      if (!ParseInt(need_value("--warmup"), warmup)) {
        std::cerr << "Invalid --warmup\n";
        return 2;
      }
    } else if (arg == "--reps") {
      if (!ParseInt(need_value("--reps"), reps)) {
        std::cerr << "Invalid --reps\n";
        return 2;
      }
    } else if (arg == "--seed") {
      seed = ParseU64OrDie(need_value("--seed"), "--seed");
    } else if (arg == "--auto-zeta") {
      const std::string mode = need_value("--auto-zeta");
      if (mode != "teich") {
        std::cerr << "Unsupported --auto-zeta mode\n";
        return 2;
      }
      auto_zeta_teich = true;
    } else if (arg == "--ring-mod") {
      if (!ParseZZ(need_value("--ring-mod"), spec.scalar_modulus)) {
        std::cerr << "Invalid --ring-mod\n";
        return 2;
      }
    } else if (arg == "--ring-p") {
      if (!ParseZZ(need_value("--ring-p"), spec.base_prime)) {
        std::cerr << "Invalid --ring-p\n";
        return 2;
      }
    } else if (arg == "--ring-F") {
      spec.F_coeffs = ParseCoeffList(need_value("--ring-F"));
    } else if (arg == "--ring-zeta") {
      spec.zeta_coeffs = ParseCoeffList(need_value("--ring-zeta"));
    } else {
      std::cerr << "Unknown argument: " << arg << "\n";
      return 2;
    }
  }

  try {
    if (ell < kappa) {
      LogicError("main: ell must be >= kappa");
    }
    if (queries < 0) {
      LogicError("main: queries must be non-negative");
    }
    if (c <= 0) {
      LogicError("main: c must be positive");
    }

    ZZ_pPush mod_push(spec.scalar_modulus);
    ValidateMonic(spec.F_coeffs, spec.scalar_modulus, "ring F");
    const ZZ_pX F = BuildZZpX(spec.F_coeffs);
    ZZ_pEPush ext_push(F);

    ZZ p_base;
    long k_base = 0;
    DeduceBasePrimeAndExponent(spec, p_base, k_base);
    ZZ_pE zeta;
    if (auto_zeta_teich) {
      zeta = FindTeichmullerGenerator(p_base, k_base, NTL::deg(F), F);
    } else {
      zeta = BuildZZpE(spec.zeta_coeffs);
    }

    const basefold::RingSwitchPCSParams params =
        BuildRingSwitchParams(c, ell, kappa, spec, F, zeta);
    const vec_ZZ_pE t_table =
        MakeDeterministicBaseRingTable(Pow2Checked(ell), seed);
    const std::vector<ZZ_pE> z =
        MakeDeterministicPoint(ell, seed ^ 0xabcddcbaULL);
    const ZZ_pE claimed_s = basefold::EvalMultilinearMonomialCoeffs(
        basefold::BooleanHypercubeTableToMonomialCoeffs(t_table), z);

    const BenchResult result =
        RunEvalBenchmark(params, t_table, z, claimed_s, queries, warmup, reps);
    PrintResult(c, ell, kappa, queries, warmup, reps, result);
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << "\n";
    return 2;
  }
}
