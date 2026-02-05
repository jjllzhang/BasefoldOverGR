#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

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

long Pow2Checked(long e) {
  if (e < 0) LogicError("Pow2Checked: negative exponent");
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
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

long PolyDegree(const std::vector<long> &coeffs, long mod) {
  if (coeffs.empty())
    LogicError("PolyDegree: empty polynomial");
  long last = static_cast<long>(coeffs.size()) - 1;
  while (last > 0 &&
         NormalizeMod(coeffs[static_cast<std::size_t>(last)], mod) == 0) {
    --last;
  }
  return last;
}

std::uint64_t FixedCoeffByteWidth(long mod) {
  if (mod <= 1)
    LogicError("FixedCoeffByteWidth: mod must be > 1");
  std::uint64_t x = static_cast<std::uint64_t>(mod - 1);
  int bits = 0;
  while (x > 0) {
    ++bits;
    x >>= 1;
  }
  return static_cast<std::uint64_t>((bits + 7) / 8);
}

std::uint64_t FixedFieldElementBytes(long mod, long ext_degree) {
  if (ext_degree <= 0)
    LogicError("FixedFieldElementBytes: ext_degree must be > 0");
  const std::uint64_t coeff_bytes = FixedCoeffByteWidth(mod);
  return 8ULL +
         static_cast<std::uint64_t>(ext_degree) * (8ULL + coeff_bytes);
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
  total += 8;  // index
  total += field_elem_bytes;
  total += static_cast<unsigned __int128>(height) * kHashBytes;  // sibling hashes

  if (total > std::numeric_limits<std::uint64_t>::max())
    LogicError("MerkleOpeningBytes: overflow");
  return static_cast<std::uint64_t>(total);
}

std::uint64_t EstimateEvalProofSizeFormulaBytes(long mod, long c, long d,
                                                long num_queries,
                                                long ext_degree) {
  static constexpr std::uint64_t kHashBytes = 32;
  if (c <= 0)
    LogicError("EstimateEvalProofSizeFormulaBytes: c must be > 0");
  if (d < 0)
    LogicError("EstimateEvalProofSizeFormulaBytes: d must be >= 0");
  if (num_queries < 0)
    LogicError("EstimateEvalProofSizeFormulaBytes: queries must be >= 0");

  const std::uint64_t fe_bytes = FixedFieldElementBytes(mod, ext_degree);

  unsigned __int128 total = 0;
  total += static_cast<unsigned __int128>(d + 1) * kHashBytes;          // roots
  total += static_cast<unsigned __int128>(d) * 3ULL * fe_bytes;         // sumcheck polys
  total += static_cast<unsigned __int128>(c) * fe_bytes;                // pi0_full

  std::vector<std::uint64_t> open_bytes;
  open_bytes.resize(static_cast<std::size_t>(d + 1));

  std::uint64_t leaf_count = static_cast<std::uint64_t>(c);
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
      << "  bench_pcs_proof_size [--mode field|ring|both] [--c <int>] [--d <int>]\n"
      << "                     [--queries <int>] [--seed <u64>] [--formula]\n"
      << "                     [--field-mod <int>] [--field-F <a0,a1,...>] [--field-zeta <b0,b1,...>]\n"
      << "                     [--ring-mod <int>]  [--ring-p <int>] [--ring-F <a0,a1,...>] [--ring-zeta <b0,b1,...>]\n\n"
      << "Notes:\n"
      << "  KB is KiB (1024 bytes).\n"
      << "  By default, prover uses BaseFoldPCSProveEval (includes parameter/length checks and claimed_y == f(z)).\n\n"
      << "  With --formula, it does NOT run the prover. It estimates proof size from (c,d,queries)\n"
      << "  assuming fixed-size field elements and sha256-based Merkle hashing.\n\n"
      << "Examples:\n"
      << "  # GF(2^2) with F(x)=x^2+x+1 and zeta=x\n"
      << "  bench_pcs_proof_size --mode field --field-mod 2 --field-F 1,1,1 --field-zeta 0,1 --d 16 --queries 4\n"
      << "  # GR(4,2) with the same extension polynomial and zeta=x\n"
      << "  bench_pcs_proof_size --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --ring-zeta 0,1 --d 16 --queries 4\n";
}

void RunOneContext(const ContextSpec &spec, long c, long d, long num_queries,
                   bool formula_only,
                   std::uint64_t seed) {
  if (spec.mod <= 1) LogicError("RunOneContext: modulus must be > 1");
  if (c <= 0) LogicError("RunOneContext: c must be > 0");
  if (d < 0) LogicError("RunOneContext: d must be >= 0");
  if (num_queries < 0) LogicError("RunOneContext: queries must be >= 0");

  if (formula_only) {
    ValidateMonic(spec.F_coeffs, spec.mod, "F");
    const long r = PolyDegree(spec.F_coeffs, spec.mod);

    const std::uint64_t bytes =
        EstimateEvalProofSizeFormulaBytes(spec.mod, c, d, num_queries, r);
    const double kb = static_cast<double>(bytes) / 1024.0;

    std::cout << "\n[" << spec.label << "] c=" << c << " k0=1 d=" << d
              << "  mod=" << spec.mod << "  queries=" << num_queries
              << "  (formula)\n";
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "  proof size  " << kb << " KB  (" << bytes << " B)\n";
    return;
  }

  const ZZ modulus = to_ZZ(spec.mod);
  ZZ_pPush mod_push(modulus);

  ValidateMonic(spec.F_coeffs, spec.mod, "F");
  const ZZ_pX F = BuildZZpX(spec.F_coeffs);
  ZZ_pEPush e_push(F);

  const ZZ_pX zpoly = BuildZZpX(spec.zeta_coeffs);
  ZZ_pE zeta;
  conv(zeta, zpoly);

  const basefold::FoldableCodeParams params =
      BuildParams_k0_1(c, d, to_ZZ(spec.prime_p), zeta);

  const long k_d = Pow2Checked(d);  // k0==1
  const vec_ZZ_pE f_coeffs = MakeDeterministicCoefficients(k_d, seed);
  const std::vector<ZZ_pE> z = MakeDeterministicPoint(d, seed ^ 0xdeadbeefULL);

  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);

  const std::uint64_t bytes = basefold::BaseFoldPCSEvalProofSizeBytes(proof);
  const double kb = static_cast<double>(bytes) / 1024.0;

  std::cout << "\n[" << spec.label << "] c=" << c << " k0=1 d=" << d
            << "  mod=" << spec.mod << "  queries=" << num_queries << "\n";
  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  proof size  " << kb << " KB  (" << bytes << " B)\n";
}

}  // namespace

int main(int argc, char **argv) {
  long d = 16;
  long c = 2;
  long num_queries = 4;
  bool formula_only = false;
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
    } else if (arg == "--queries") {
      if (!ParseLong(NeedValue("--queries"), num_queries) || num_queries < 0) {
        std::cerr << "Invalid --queries\n";
        return 2;
      }
    } else if (arg == "--formula") {
      formula_only = true;
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
      RunOneContext(field, c, d, num_queries, formula_only, seed);
    if (do_ring)
      RunOneContext(ring, c, d, num_queries, formula_only, seed);
  } catch (const std::exception &e) {
    std::cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    std::cerr << "Unhandled non-std exception\n";
    return 2;
  }

  return 0;
}
