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

#include "bench_common_helpers.hpp"
#include "GaloisRing/PrimitiveElement.hpp"
#include "PCS/BaseFold/BaseFoldPCS.hpp"
#include "PCS/Common/Hash.hpp"

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

using basefold_bench_common::BuildZZpX;
using basefold_bench_common::ComputeStats;
using basefold_bench_common::DeduceBasePrimeAndExponent;
using basefold_bench_common::MakeDeterministicCoefficients;
using basefold_bench_common::MsSince;
using basefold_bench_common::NormalizeMod;
using basefold_bench_common::ParseCoeffList;
using basefold_bench_common::ParseInt;
using basefold_bench_common::ParseLong;
using basefold_bench_common::ParseU64OrDie;
using basefold_bench_common::ParseZZ;
using basefold_bench_common::Pow2Checked;
using basefold_bench_common::Stats;
using basefold_bench_common::ValidateMonic;

struct ContextSpec {
  std::string label;
  ZZ scalar_modulus = ZZ(0);  // ZZ_p modulus (p for fields, p^s for rings)
  ZZ base_prime = ZZ(0);      // optional: the prime p (only used by checked paths)
  std::vector<ZZ> F_coeffs;     // extension modulus polynomial coefficients
  std::vector<ZZ> zeta_coeffs;  // ζ element coefficients
};

basefold::FoldableCodeParams BuildParams_k0_1(long c, long d, const ZZ &prime_p,
                                             const ZZ_pE &zeta) {
  return basefold_bench_common::BuildFoldableParamsK0Eq1(
      c, d, prime_p, zeta, "BuildParams_k0_1");
}

struct BenchResult {
  Stats encode;
  Stats top_merkle_build;
  Stats commit;
  std::uint64_t anti_opt_checksum = 0;
};

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
  std::vector<double> top_merkle_build_ms;
  std::vector<double> commit_ms;
  encode_ms.reserve(static_cast<std::size_t>(reps));
  top_merkle_build_ms.reserve(static_cast<std::size_t>(reps));
  commit_ms.reserve(static_cast<std::size_t>(reps));

  std::uint64_t anti_opt_checksum = 0;
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
      anti_opt_checksum ^= static_cast<std::uint64_t>(pi_d[0] != ZZ_pE(0));
    }
    if (!root_d.empty()) {
      anti_opt_checksum ^= static_cast<std::uint64_t>(root_d[0]);
    }

    if (iter >= 0) {
      encode_ms.push_back(MsSince(t0, t1));
      top_merkle_build_ms.push_back(MsSince(t1, t2));
      commit_ms.push_back(MsSince(t0, t2));
    }
  }

  BenchResult out;
  out.encode = ComputeStats(encode_ms);
  out.top_merkle_build = ComputeStats(top_merkle_build_ms);
  out.commit = ComputeStats(commit_ms);
  out.anti_opt_checksum = anti_opt_checksum;
  return out;
}

void PrintResult(const std::string &label, const ZZ &scalar_modulus, long c,
                 long k0, long d, long k_d, long n_d, int warmup, int reps,
                 const BenchResult &r) {
  std::cout << "\n[" << label << "] c=" << c << " k0=" << k0 << " d=" << d
            << "  mod=" << scalar_modulus << "  k_d=" << k_d
            << "  n_d=" << n_d
            << "  warmup=" << warmup << " reps=" << reps << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  hash backend " << basefold::SelectedHashBackendName() << "\n";
  std::cout << "  encode-only mean " << r.encode.mean_ms << " ms  (min "
            << r.encode.min_ms << ", max " << r.encode.max_ms << ")\n";
  std::cout << "  top-merkle-build mean " << r.top_merkle_build.mean_ms
            << " ms  (min " << r.top_merkle_build.min_ms << ", max "
            << r.top_merkle_build.max_ms << ")\n";
  std::cout << "  commit     mean " << r.commit.mean_ms << " ms  (min "
            << r.commit.min_ms << ", max " << r.commit.max_ms << ")\n";
  std::cout << "  anti-opt checksum " << r.anti_opt_checksum << "\n";
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
      << "  Auxiliary lines print encode-only and top-merkle-build-only splits.\n\n"
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

void RunOneContext(const ContextSpec &spec, long c, long k0, long d, int warmup,
                   int reps, bool auto_zeta_teich, std::uint64_t seed) {
  if (spec.scalar_modulus <= 1) {
    LogicError("RunOneContext: modulus must be > 1");
  }

  const ZZ modulus = spec.scalar_modulus;
  ZZ_pPush mod_push(modulus);

  ValidateMonic(spec.F_coeffs, spec.scalar_modulus, "F");
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
          c, d, (spec.base_prime > 1) ? spec.base_prime : spec.scalar_modulus,
          zeta);
    }
    return BuildParams_k0_gt1(
        c, k0, d,
        (spec.base_prime > 1) ? spec.base_prime : spec.scalar_modulus, zeta,
        seed);
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
  PrintResult(spec.label, spec.scalar_modulus, c, k0, d, k_d, n_d, warmup,
              reps, r);
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
  field.scalar_modulus = to_ZZ(2);
  field.base_prime = ZZ(0);
  field.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};      // x^2 + x + 1
  field.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};      // x

  ContextSpec ring;
  ring.label = "Ring";
  ring.scalar_modulus = to_ZZ(4);
  ring.base_prime = to_ZZ(2);
  ring.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};       // x^2 + x + 1
  ring.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};       // x

  for (int i = 1; i < argc; ++i) {
    const std::string arg(argv[i]);
    auto RequireNextArgValue = [&](const char *flag) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << flag << "\n";
        std::exit(2);
      }
      return argv[++i];
    };

    if (arg == "--mode") {
      const std::string m = RequireNextArgValue("--mode");
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
      if (!ParseLong(RequireNextArgValue("--d"), d) || d < 0) {
        std::cerr << "Invalid --d\n";
        return 2;
      }
    } else if (arg == "--c") {
      if (!ParseLong(RequireNextArgValue("--c"), c) || c <= 0) {
        std::cerr << "Invalid --c\n";
        return 2;
      }
    } else if (arg == "--k0") {
      if (!ParseLong(RequireNextArgValue("--k0"), k0) || k0 <= 0) {
        std::cerr << "Invalid --k0\n";
        return 2;
      }
    } else if (arg == "--warmup") {
      if (!ParseInt(RequireNextArgValue("--warmup"), warmup) || warmup < 0) {
        std::cerr << "Invalid --warmup\n";
        return 2;
      }
    } else if (arg == "--reps") {
      if (!ParseInt(RequireNextArgValue("--reps"), reps) || reps <= 0) {
        std::cerr << "Invalid --reps\n";
        return 2;
      }
    } else if (arg == "--seed") {
      seed = ParseU64OrDie(RequireNextArgValue("--seed"), "--seed");
    } else if (arg == "--auto-zeta") {
      const std::string mode = RequireNextArgValue("--auto-zeta");
      if (mode == "teich") {
        auto_zeta_teich = true;
      } else {
        std::cerr << "Invalid --auto-zeta (expected teich)\n";
        return 2;
      }
    } else if (arg == "--field-mod") {
      if (!ParseZZ(RequireNextArgValue("--field-mod"), field.scalar_modulus) ||
          field.scalar_modulus <= 1) {
        std::cerr << "Invalid --field-mod\n";
        return 2;
      }
    } else if (arg == "--field-F") {
      field.F_coeffs = ParseCoeffList(RequireNextArgValue("--field-F"));
    } else if (arg == "--field-zeta") {
      field.zeta_coeffs = ParseCoeffList(RequireNextArgValue("--field-zeta"));
    } else if (arg == "--ring-mod") {
      if (!ParseZZ(RequireNextArgValue("--ring-mod"), ring.scalar_modulus) ||
          ring.scalar_modulus <= 1) {
        std::cerr << "Invalid --ring-mod\n";
        return 2;
      }
    } else if (arg == "--ring-p") {
      if (!ParseZZ(RequireNextArgValue("--ring-p"), ring.base_prime) ||
          ring.base_prime <= 1) {
        std::cerr << "Invalid --ring-p\n";
        return 2;
      }
    } else if (arg == "--ring-F") {
      ring.F_coeffs = ParseCoeffList(RequireNextArgValue("--ring-F"));
    } else if (arg == "--ring-zeta") {
      ring.zeta_coeffs = ParseCoeffList(RequireNextArgValue("--ring-zeta"));
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
