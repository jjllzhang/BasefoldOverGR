#include <NTL/ZZ.h>

#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using NTL::LogicError;
using NTL::to_ZZ;
using NTL::ZZ;

namespace {

struct ContextSpec {
  std::string label;
  ZZ mod = ZZ(0);      // ZZ_p modulus (p for fields, p^s for rings)
  ZZ prime_p = ZZ(0);  // optional metadata (not needed for byte-size estimation)
  std::vector<ZZ> F_coeffs;
};

bool IsPowerOfTwoLong(long n) { return n > 0 && (n & (n - 1)) == 0; }

long Log2ExactPowerOfTwoLong(long n) {
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

bool ParseZZ(const char *s, ZZ &out) {
  if (!s) return false;
  return ParseZZString(std::string(s), out);
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

    ZZ v;
    if (!ParseZZString(token, v)) LogicError("ParseCoeffList: bad integer token");

    out.push_back(v);
    pos = (comma == std::string::npos) ? s.size() : (comma + 1);
  }
  if (out.empty()) {
    LogicError("ParseCoeffList: empty list");
  }
  return out;
}

ZZ NormalizeMod(const ZZ &x, const ZZ &mod) {
  if (mod <= 0) LogicError("NormalizeMod: mod must be positive");
  ZZ r = x % mod;
  if (r < 0) r += mod;
  return r;
}

void ValidateMonic(const std::vector<ZZ> &coeffs, const ZZ &mod, const char *what) {
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

long PolyDegree(const std::vector<ZZ> &coeffs, const ZZ &mod) {
  if (coeffs.empty()) {
    LogicError("PolyDegree: empty polynomial");
  }
  long last = static_cast<long>(coeffs.size()) - 1;
  while (last > 0 &&
         NormalizeMod(coeffs[static_cast<std::size_t>(last)], mod) == 0) {
    --last;
  }
  return last;
}

std::uint64_t CheckedToU64(unsigned __int128 value, const char *where) {
  if (value > std::numeric_limits<std::uint64_t>::max()) {
    const std::string msg = std::string(where) + ": overflow";
    LogicError(msg.c_str());
  }
  return static_cast<std::uint64_t>(value);
}

std::uint64_t FixedCoeffByteWidth(const ZZ &mod) {
  if (mod <= 1) {
    LogicError("FixedCoeffByteWidth: mod must be > 1");
  }
  const long bits = NTL::NumBits(mod - 1);
  return static_cast<std::uint64_t>((bits + 7) / 8);
}

std::uint64_t FixedFieldElementBytes(const ZZ &mod, long ext_degree) {
  if (ext_degree <= 0) {
    LogicError("FixedFieldElementBytes: ext_degree must be > 0");
  }
  const std::uint64_t coeff_bytes = FixedCoeffByteWidth(mod);
  return 8ULL +
         static_cast<std::uint64_t>(ext_degree) * (8ULL + coeff_bytes);
}

std::size_t MerkleHeight(std::uint64_t leaf_count) {
  if (leaf_count <= 1) {
    return 0;
  }
  std::size_t height = 0;
  std::uint64_t n = leaf_count;
  while (n > 1) {
    if (n & 1ULL) {
      n += 1;
    }
    n /= 2;
    ++height;
  }
  return height;
}

std::uint64_t MerkleOpeningBytes(std::uint64_t leaf_count,
                                 std::uint64_t field_elem_bytes) {
  static constexpr std::uint64_t kHashBytes = 32;
  unsigned __int128 total = 0;
  total += 8;  // leaf index
  total += field_elem_bytes;
  total += static_cast<unsigned __int128>(MerkleHeight(leaf_count)) * kHashBytes;
  return CheckedToU64(total, "MerkleOpeningBytes");
}

std::uint64_t QueryIndexBytes(std::uint64_t upper_bound) {
  if (upper_bound <= 1) {
    return 0;
  }
  std::uint64_t t = upper_bound - 1;
  int bits = 0;
  while (t > 0) {
    ++bits;
    t >>= 1;
  }
  return static_cast<std::uint64_t>((bits + 7) / 8);
}

struct CommunicationEstimate {
  std::uint64_t point_dim = 0;
  std::uint64_t n0 = 0;
  std::uint64_t n_last = 0;
  std::uint64_t field_elem_bytes = 0;
  std::uint64_t query_index_bytes = 0;

  std::uint64_t roots_p2v = 0;
  std::uint64_t sumcheck_p2v = 0;
  std::uint64_t pi0_full_p2v = 0;
  std::uint64_t openings_p2v = 0;
  std::uint64_t total_p2v = 0;

  std::uint64_t challenges_v2p = 0;
  std::uint64_t query_indices_v2p = 0;
  std::uint64_t total_v2p_interactive = 0;

  std::uint64_t total_interactive = 0;
  std::uint64_t total_fs = 0;
};

CommunicationEstimate EstimateCommunicationBytes(const ZZ &mod, long c, long k0, long d,
                                                 long num_queries, long ext_degree) {
  static constexpr std::uint64_t kHashBytes = 32;
  if (mod <= 1) LogicError("EstimateCommunicationBytes: mod must be > 1");
  if (c <= 0) LogicError("EstimateCommunicationBytes: c must be > 0");
  if (k0 <= 0) LogicError("EstimateCommunicationBytes: k0 must be > 0");
  if (!IsPowerOfTwoLong(k0)) {
    LogicError("EstimateCommunicationBytes: k0 must be a power of two");
  }
  if (d < 0) LogicError("EstimateCommunicationBytes: d must be >= 0");
  if (num_queries < 0) {
    LogicError("EstimateCommunicationBytes: queries must be >= 0");
  }

  CommunicationEstimate out;
  out.point_dim = static_cast<std::uint64_t>(d + Log2ExactPowerOfTwoLong(k0));
  out.field_elem_bytes = FixedFieldElementBytes(mod, ext_degree);

  const std::uint64_t c_u64 = static_cast<std::uint64_t>(c);
  const std::uint64_t k0_u64 = static_cast<std::uint64_t>(k0);
  if (k0_u64 > std::numeric_limits<std::uint64_t>::max() / c_u64) {
    LogicError("EstimateCommunicationBytes: overflow in n0");
  }
  out.n0 = c_u64 * k0_u64;

  out.roots_p2v = CheckedToU64(
      static_cast<unsigned __int128>(d + 1) * kHashBytes, "roots_p2v");
  out.sumcheck_p2v = CheckedToU64(
      static_cast<unsigned __int128>(d) * 3ULL * out.field_elem_bytes,
      "sumcheck_p2v");
  out.pi0_full_p2v =
      CheckedToU64(static_cast<unsigned __int128>(out.n0) * out.field_elem_bytes,
                   "pi0_full_p2v");

  std::vector<std::uint64_t> open_bytes_by_level;
  open_bytes_by_level.resize(static_cast<std::size_t>(d + 1));
  std::vector<std::uint64_t> leaf_count_by_level;
  leaf_count_by_level.resize(static_cast<std::size_t>(d + 1));

  std::uint64_t leaf_count = out.n0;
  for (long level = 0; level <= d; ++level) {
    leaf_count_by_level[static_cast<std::size_t>(level)] = leaf_count;
    open_bytes_by_level[static_cast<std::size_t>(level)] =
        MerkleOpeningBytes(leaf_count, out.field_elem_bytes);
    if (level < d) {
      if (leaf_count > std::numeric_limits<std::uint64_t>::max() / 2ULL) {
        LogicError("EstimateCommunicationBytes: overflow in n_i");
      }
      leaf_count *= 2ULL;
    }
  }

  if (d > 0) {
    out.n_last = leaf_count_by_level[static_cast<std::size_t>(d - 1)];
  } else {
    out.n_last = out.n0;
  }

  unsigned __int128 per_query_p2v = 0;
  for (long i = 0; i < d; ++i) {
    per_query_p2v +=
        open_bytes_by_level[static_cast<std::size_t>(i)];  // folded at level i
    per_query_p2v +=
        2ULL * open_bytes_by_level[static_cast<std::size_t>(i + 1)];  // left/right
  }
  out.openings_p2v = CheckedToU64(
      static_cast<unsigned __int128>(num_queries) * per_query_p2v, "openings_p2v");

  out.total_p2v = CheckedToU64(
      static_cast<unsigned __int128>(out.roots_p2v) + out.sumcheck_p2v +
          out.pi0_full_p2v + out.openings_p2v,
      "total_p2v");

  out.challenges_v2p = CheckedToU64(
      static_cast<unsigned __int128>(d) * out.field_elem_bytes, "challenges_v2p");

  out.query_index_bytes = 0;
  out.query_indices_v2p = 0;
  if (d > 0) {
    out.query_index_bytes = QueryIndexBytes(out.n_last);
    out.query_indices_v2p = CheckedToU64(
        static_cast<unsigned __int128>(num_queries) * out.query_index_bytes,
        "query_indices_v2p");
  }

  out.total_v2p_interactive = CheckedToU64(
      static_cast<unsigned __int128>(out.challenges_v2p) + out.query_indices_v2p,
      "total_v2p_interactive");
  out.total_interactive = CheckedToU64(
      static_cast<unsigned __int128>(out.total_p2v) + out.total_v2p_interactive,
      "total_interactive");
  out.total_fs = out.total_p2v;
  return out;
}

void PrintHelp() {
  std::cout
      << "bench_pcs_communication (estimate PCS prover/verifier communication)\n\n"
      << "Usage:\n"
      << "  bench_pcs_communication [--mode field|ring|both] [--c <int>] [--k0 <int>] [--d <int>]\n"
      << "                          [--queries <int>] [--field-mod <decimal-int>] [--field-F <a0,a1,...>]\n"
      << "                          [--ring-mod <decimal-int>] [--ring-p <decimal-int>] [--ring-F <a0,a1,...>]\n\n"
      << "Estimation model:\n"
      << "  - No prover execution; formula-only from parameters.\n"
      << "  - Prover->Verifier follows BaseFoldPCSEvalProof payload layout:\n"
      << "      roots + sumcheck polys + pi0_full + query openings.\n"
      << "  - Verifier->Prover (interactive-equivalent):\n"
      << "      d field challenges r_i + query indices mu.\n"
      << "  - Fiat-Shamir (current implementation): V->P is 0, total equals proof payload.\n"
      << "  - Not counting external/public-input transport (commitment C, point z, claimed y).\n\n"
      << "Notes:\n"
      << "  - Byte model assumes fixed-size ZZ_pE serialization used in bench formula mode.\n"
      << "  - KB is KiB (1024 bytes).\n\n"
      << "Examples:\n"
      << "  bench_pcs_communication --mode field --field-mod 2 --field-F 1,1,1 --d 16 --queries 4\n"
      << "  bench_pcs_communication --mode ring  --ring-mod 4 --ring-p 2 --ring-F 1,1,1 --d 16 --queries 4\n";
}

void PrintEstimate(const ContextSpec &spec, long c, long k0, long d,
                   long num_queries) {
  ValidateMonic(spec.F_coeffs, spec.mod, "F");
  const long ext_degree = PolyDegree(spec.F_coeffs, spec.mod);

  const CommunicationEstimate est =
      EstimateCommunicationBytes(spec.mod, c, k0, d, num_queries, ext_degree);

  const auto ToKB = [](std::uint64_t bytes) -> double {
    return static_cast<double>(bytes) / 1024.0;
  };

  std::cout << "\n[" << spec.label << "] c=" << c << " k0=" << k0 << " d=" << d
            << "  mod=" << spec.mod << "  queries=" << num_queries << "\n";
  std::cout << "  point dim (d+log2(k0)) = " << est.point_dim << "\n";
  std::cout << "  n0 = c*k0 = " << est.n0;
  if (d > 0) {
    std::cout << " , n_(d-1) = " << est.n_last;
  }
  std::cout << "\n";
  std::cout << "  field elem bytes = " << est.field_elem_bytes << "\n";
  if (d > 0) {
    std::cout << "  mu index bytes   = " << est.query_index_bytes
              << " (interactive-equivalent)\n";
  }

  std::cout << std::fixed << std::setprecision(3);
  std::cout << "  prover -> verifier:\n";
  std::cout << "    roots     " << ToKB(est.roots_p2v) << " KB  ("
            << est.roots_p2v << " B)\n";
  std::cout << "    sumcheck  " << ToKB(est.sumcheck_p2v) << " KB  ("
            << est.sumcheck_p2v << " B)\n";
  std::cout << "    pi0_full  " << ToKB(est.pi0_full_p2v) << " KB  ("
            << est.pi0_full_p2v << " B)\n";
  std::cout << "    openings  " << ToKB(est.openings_p2v) << " KB  ("
            << est.openings_p2v << " B)\n";
  std::cout << "    total     " << ToKB(est.total_p2v) << " KB  ("
            << est.total_p2v << " B)\n";

  std::cout << "  verifier -> prover (interactive-equivalent):\n";
  std::cout << "    r_i       " << ToKB(est.challenges_v2p) << " KB  ("
            << est.challenges_v2p << " B)\n";
  std::cout << "    mu        " << ToKB(est.query_indices_v2p) << " KB  ("
            << est.query_indices_v2p << " B)\n";
  std::cout << "    total     " << ToKB(est.total_v2p_interactive) << " KB  ("
            << est.total_v2p_interactive << " B)\n";

  std::cout << "  total interactive communication:\n";
  std::cout << "    total     " << ToKB(est.total_interactive) << " KB  ("
            << est.total_interactive << " B)\n";

  std::cout << "  total Fiat-Shamir communication (current code path):\n";
  std::cout << "    total     " << ToKB(est.total_fs) << " KB  ("
            << est.total_fs << " B)\n";
}

}  // namespace

int main(int argc, char **argv) {
  long d = 16;
  long c = 2;
  long k0 = 1;
  long num_queries = 4;

  bool do_field = true;
  bool do_ring = true;

  ContextSpec field;
  field.label = "Field";
  field.mod = to_ZZ(2);
  field.prime_p = ZZ(0);
  field.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};  // x^2 + x + 1

  ContextSpec ring;
  ring.label = "Ring";
  ring.mod = to_ZZ(4);
  ring.prime_p = to_ZZ(2);
  ring.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};  // x^2 + x + 1

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
    } else if (arg == "--field-mod") {
      if (!ParseZZ(NeedValue("--field-mod"), field.mod) || field.mod <= 1) {
        std::cerr << "Invalid --field-mod\n";
        return 2;
      }
    } else if (arg == "--field-F") {
      field.F_coeffs = ParseCoeffList(NeedValue("--field-F"));
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
    if (do_field) {
      PrintEstimate(field, c, k0, d, num_queries);
    }
    if (do_ring) {
      PrintEstimate(ring, c, k0, d, num_queries);
    }
  } catch (const std::exception &e) {
    std::cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    std::cerr << "Unhandled non-std exception\n";
    return 2;
  }

  return 0;
}
