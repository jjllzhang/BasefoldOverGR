#include "BaseFold/IOPP.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>

#include <openssl/evp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "GaloisRing/Inverse.hpp"

using NTL::BytesFromZZ;
using NTL::coeff;
using NTL::LogicError;
using NTL::mat_ZZ_pE;
using NTL::mul;
using NTL::NumBytes;
using NTL::rep;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pX;

namespace basefold {
namespace {

long Pow2Checked(long e) {
  if (e < 0)
    LogicError("Pow2Checked: negative exponent");
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
}

void ValidateParamsOrThrow(const FoldableCodeParams &params) {
  (void)MessageLength(params);
}

void AppendU64(Bytes &out, std::uint64_t v) {
  for (int i = 7; i >= 0; --i) {
    out.push_back(static_cast<Byte>((v >> (8 * i)) & 0xff));
  }
}

Bytes ZZToBytes(const ZZ &x) {
  const long n = NumBytes(x);
  Bytes out;
  out.resize(static_cast<std::size_t>(n));
  if (n > 0) {
    BytesFromZZ(reinterpret_cast<unsigned char *>(out.data()), x, n);
  }
  return out;
}

Bytes SerializeFieldElement(const FieldElement &x) {
  const long r = ZZ_pE::degree();
  if (r <= 0)
    LogicError("SerializeFieldElement: invalid extension degree");

  const ZZ_pX poly = rep(x);
  Bytes out;
  AppendU64(out, static_cast<std::uint64_t>(r));
  for (long i = 0; i < r; ++i) {
    const ZZ c = rep(coeff(poly, i));
    const Bytes c_bytes = ZZToBytes(c);
    AppendU64(out, static_cast<std::uint64_t>(c_bytes.size()));
    out.insert(out.end(), c_bytes.begin(), c_bytes.end());
  }
  return out;
}

bool IsUnitInBaseRing(const NTL::ZZ_p &a) {
  if (a == 0) return false;
  NTL::ZZ g;
  NTL::GCD(g, NTL::rep(a), NTL::ZZ_p::modulus());
  return g == 1;
}

bool IsUnit(const ZZ_pE &a) {
  if (a == 0) return false;
  const NTL::ZZ_pX poly = rep(a);
  const long r = ZZ_pE::degree();
  for (long i = 0; i < r; ++i) {
    if (IsUnitInBaseRing(coeff(poly, i))) return true;
  }
  return false;
}

bool TryInvertUnit(ZZ_pE &inv_out, const ZZ_pE &a) {
  if (!IsUnit(a)) return false;

  const long r = ZZ_pE::degree();
  if (r <= 0) LogicError("TryInvertUnit: invalid extension degree");

  if (r == 1) {
    const NTL::ZZ a_rep = NTL::rep(coeff(rep(a), 0));
    NTL::ZZ inv_rep;
    if (NTL::InvModStatus(inv_rep, a_rep, NTL::ZZ_p::modulus()) != 0) return false;
    NTL::ZZ_p inv_base;
    NTL::conv(inv_base, inv_rep);
    NTL::ZZ_pX poly;
    NTL::clear(poly);
    NTL::SetCoeff(poly, 0, inv_base);
    NTL::conv(inv_out, poly);
    return true;
  }

  inv_out = Inv(a, r);
  if (inv_out == 0) return false;
  ZZ_pE one;
  NTL::set(one);
  return a * inv_out == one;
}

Bytes Sha256(const Byte *data, std::size_t len) {
  Bytes digest;
  digest.resize(32);
  unsigned int out_len = 0;
  const int ok = EVP_Digest(static_cast<const void *>(data), len,
                            reinterpret_cast<unsigned char *>(digest.data()),
                            &out_len, EVP_sha256(), nullptr);
  if (ok != 1) {
    LogicError("Sha256: EVP_Digest failed");
  }
  if (out_len != digest.size()) {
    LogicError("Sha256: unexpected digest size");
  }
  return digest;
}

Bytes Sha256(const Bytes &data) { return Sha256(data.data(), data.size()); }

Bytes HashWithPrefix(Byte prefix, const Bytes &payload) {
  Bytes in;
  in.reserve(1 + payload.size());
  in.push_back(prefix);
  in.insert(in.end(), payload.begin(), payload.end());
  return Sha256(in);
}

Bytes HashNode(const Bytes &left, const Bytes &right) {
  Bytes in;
  in.reserve(1 + left.size() + right.size());
  in.push_back(static_cast<Byte>(0x01));
  in.insert(in.end(), left.begin(), left.end());
  in.insert(in.end(), right.begin(), right.end());
  return Sha256(in);
}

Bytes HashLeaf(long index, const FieldElement &value) {
  Bytes in;
  in.reserve(1 + 8 + 64);
  in.push_back(static_cast<Byte>(0x00));
  AppendU64(in, static_cast<std::uint64_t>(index));
  const Bytes enc = SerializeFieldElement(value);
  in.insert(in.end(), enc.begin(), enc.end());
  return Sha256(in);
}

Bytes HashRootWithCount(long leaf_count, const Bytes &raw_root) {
  Bytes in;
  in.reserve(1 + 8 + raw_root.size());
  in.push_back(static_cast<Byte>(0x03));
  AppendU64(in, static_cast<std::uint64_t>(leaf_count));
  in.insert(in.end(), raw_root.begin(), raw_root.end());
  return Sha256(in);
}

std::size_t ExpectedMerkleHeight(long leaf_count) {
  if (leaf_count <= 0)
    return 0;
  std::size_t height = 0;
  long n = leaf_count;
  while (n > 1) {
    if (n & 1L)
      n += 1;
    n /= 2;
    ++height;
  }
  return height;
}

Bytes MerkleRootRaw(const std::vector<Bytes> &leaf_hashes) {
  if (leaf_hashes.empty())
    return HashWithPrefix(static_cast<Byte>(0x04), {});

  std::vector<Bytes> level = leaf_hashes;
  while (level.size() > 1) {
    if (level.size() % 2 == 1)
      level.push_back(level.back());
    std::vector<Bytes> next;
    next.reserve(level.size() / 2);
    for (std::size_t i = 0; i < level.size(); i += 2) {
      next.push_back(HashNode(level[i], level[i + 1]));
    }
    level = std::move(next);
  }
  return level[0];
}

bool SolveLinearSystemRref(vec_ZZ_pE &x_out, mat_ZZ_pE &aug) {
  const long m = aug.NumRows();
  const long n_plus_1 = aug.NumCols();
  if (n_plus_1 <= 0)
    LogicError("SolveLinearSystemRref: empty system");
  const long n = n_plus_1 - 1;

  x_out.SetLength(n);
  for (long i = 0; i < n; ++i)
    x_out[i] = ZZ_pE(0);

  std::vector<long> pivot_col_orig;
  pivot_col_orig.resize(static_cast<std::size_t>(m), -1);

  std::vector<long> col_perm;
  col_perm.resize(static_cast<std::size_t>(n));
  for (long c = 0; c < n; ++c)
    col_perm[static_cast<std::size_t>(c)] = c;

  long row = 0;
  for (long col = 0; col < n && row < m; ++col) {
    long pivot_row = -1;
    long pivot_col = -1;
    ZZ_pE inv_pivot;

    for (long c = col; c < n; ++c) {
      for (long r = row; r < m; ++r) {
        if (TryInvertUnit(inv_pivot, aug[r][c])) {
          pivot_row = r;
          pivot_col = c;
          break;
        }
      }
      if (pivot_row >= 0)
        break;
    }

    if (pivot_row < 0) {
      break;
    }

    if (pivot_col != col) {
      for (long r = 0; r < m; ++r) {
        std::swap(aug[r][col], aug[r][pivot_col]);
      }
      std::swap(col_perm[static_cast<std::size_t>(col)],
                col_perm[static_cast<std::size_t>(pivot_col)]);
    }

    if (pivot_row != row) {
      for (long c = 0; c < n_plus_1; ++c) {
        std::swap(aug[row][c], aug[pivot_row][c]);
      }
    }

    for (long c = col; c < n_plus_1; ++c) {
      aug[row][c] *= inv_pivot;
    }

    for (long r = 0; r < m; ++r) {
      if (r == row)
        continue;
      if (aug[r][col] == 0)
        continue;
      const ZZ_pE factor = aug[r][col];
      for (long c = col; c < n_plus_1; ++c) {
        aug[r][c] -= factor * aug[row][c];
      }
    }

    pivot_col_orig[static_cast<std::size_t>(row)] =
        col_perm[static_cast<std::size_t>(col)];
    ++row;
  }

  for (long r = 0; r < m; ++r) {
    bool all_zero = true;
    for (long c = 0; c < n; ++c) {
      if (aug[r][c] != 0) {
        all_zero = false;
        break;
      }
    }
    if (all_zero && aug[r][n] != 0)
      return false;
  }

  for (long r = 0; r < m; ++r) {
    const long pc_orig = pivot_col_orig[static_cast<std::size_t>(r)];
    if (pc_orig >= 0)
      x_out[pc_orig] = aug[r][n];
  }

  return true;
}

} // namespace

long MessageLengthAtLevel(const FoldableCodeParams &params, long level) {
  ValidateParamsOrThrow(params);
  if (level < 0 || level > params.d) {
    LogicError("MessageLengthAtLevel: level out of range");
  }
  const long pow2 = Pow2Checked(level);
  if (params.k0 > std::numeric_limits<long>::max() / pow2) {
    LogicError("MessageLengthAtLevel: overflow");
  }
  return params.k0 * pow2;
}

long CodewordLengthAtLevel(const FoldableCodeParams &params, long level) {
  ValidateParamsOrThrow(params);
  const long k = MessageLengthAtLevel(params, level);
  if (params.c > std::numeric_limits<long>::max() / k) {
    LogicError("CodewordLengthAtLevel: overflow");
  }
  return params.c * k;
}

void FoldingPoints(FieldElement &x_left, FieldElement &x_right,
                   const FoldableCodeParams &params, long level_i, long j) {
  ValidateParamsOrThrow(params);
  if (level_i < 0 || level_i >= params.d) {
    LogicError("FoldingPoints: level out of range");
  }
  const long n_i = CodewordLengthAtLevel(params, level_i);
  if (j < 0 || j >= n_i)
    LogicError("FoldingPoints: index out of range");

  const FieldElement &t = params.diag_T[static_cast<std::size_t>(level_i)][j];
  x_left = t;
  x_right = params.zeta * t;
}

FieldElement EvalLineAt(const FieldElement &x, const FieldElement &x1,
                        const FieldElement &y1, const FieldElement &x2,
                        const FieldElement &y2) {
  const FieldElement denom = x2 - x1;
  if (denom == 0)
    LogicError("EvalLineAt: x1 must not equal x2");
  ZZ_pE inv_denom;
  if (!TryInvertUnit(inv_denom, denom)) {
    LogicError("EvalLineAt: x2-x1 must be a unit");
  }
  return y1 + (x - x1) * (y2 - y1) * inv_denom;
}

void ProverCommitRound(Oracle &pi_i, const Oracle &pi_ip1,
                       const FieldElement &alpha_i, long level_i,
                       const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (level_i < 0 || level_i >= params.d) {
    LogicError("ProverCommitRound: level out of range");
  }

  const long n_i = CodewordLengthAtLevel(params, level_i);
  if (pi_ip1.length() != 2 * n_i) {
    LogicError("ProverCommitRound: pi_{i+1} has wrong length");
  }

  pi_i.SetLength(n_i);
  for (long j = 0; j < n_i; ++j) {
    FieldElement x1, x2;
    FoldingPoints(x1, x2, params, level_i, j);
    const FieldElement &y1 = pi_ip1[j];
    const FieldElement &y2 = pi_ip1[j + n_i];
    pi_i[j] = EvalLineAt(alpha_i, x1, y1, x2, y2);
  }
}

void ProverCommitAll(IOPPOracles &oracles, const Oracle &pi_d,
                     const IOPPChallenges &challenges,
                     const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (static_cast<long>(challenges.alphas.size()) != params.d) {
    LogicError("ProverCommitAll: challenges.alphas has wrong size");
  }

  const long n_d = CodewordLengthAtLevel(params, params.d);
  if (pi_d.length() != n_d)
    LogicError("ProverCommitAll: pi_d has wrong length");

  oracles.pi.resize(static_cast<std::size_t>(params.d + 1));
  oracles.pi[static_cast<std::size_t>(params.d)] = pi_d;

  for (long i = params.d; i-- > 0;) {
    ProverCommitRound(oracles.pi[static_cast<std::size_t>(i)],
                      oracles.pi[static_cast<std::size_t>(i + 1)],
                      challenges.alphas[static_cast<std::size_t>(i)], i,
                      params);
  }
}

IOPPQueryPlan MakeQueryPlan(long initial_mu, const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  IOPPQueryPlan plan;
  plan.initial_mu = initial_mu;
  plan.mu_by_level.resize(static_cast<std::size_t>(params.d));

  if (params.d == 0)
    return plan;

  const long n_last = CodewordLengthAtLevel(params, params.d - 1);
  if (initial_mu < 0 || initial_mu >= n_last) {
    LogicError("MakeQueryPlan: initial_mu out of range");
  }

  long mu = initial_mu;
  for (long i = params.d; i-- > 0;) {
    plan.mu_by_level[static_cast<std::size_t>(i)] = mu;
    if (i > 0) {
      const long n_prev = CodewordLengthAtLevel(params, i - 1);
      if (mu >= n_prev)
        mu -= n_prev;
    }
  }

  return plan;
}

bool VerifyQueryFromOpenings(const IOPPQueryPlan &plan,
                             const IOPPChallenges &challenges,
                             const IOPPQueryOpenings &openings,
                             const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (static_cast<long>(challenges.alphas.size()) != params.d)
    return false;
  if (static_cast<long>(plan.mu_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(openings.left.size()) != params.d)
    return false;
  if (static_cast<long>(openings.right.size()) != params.d)
    return false;
  if (static_cast<long>(openings.folded.size()) != params.d)
    return false;

  for (long i = params.d; i-- > 0;) {
    const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
    const long n_i = CodewordLengthAtLevel(params, i);
    if (mu < 0 || mu >= n_i)
      return false;

    FieldElement x1, x2;
    FoldingPoints(x1, x2, params, i, mu);
    const FieldElement &y1 = openings.left[static_cast<std::size_t>(i)];
    const FieldElement &y2 = openings.right[static_cast<std::size_t>(i)];
    const FieldElement expected = EvalLineAt(
        challenges.alphas[static_cast<std::size_t>(i)], x1, y1, x2, y2);
    if (expected != openings.folded[static_cast<std::size_t>(i)])
      return false;
  }

  return IsCodewordC0(openings.pi0_full, params);
}

bool VerifyQueryFromOracles(const IOPPQueryPlan &plan,
                            const IOPPChallenges &challenges,
                            const IOPPOracles &oracles,
                            const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (static_cast<long>(challenges.alphas.size()) != params.d)
    return false;
  if (static_cast<long>(plan.mu_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(oracles.pi.size()) != params.d + 1)
    return false;

  for (long i = 0; i <= params.d; ++i) {
    if (oracles.pi[static_cast<std::size_t>(i)].length() !=
        CodewordLengthAtLevel(params, i)) {
      return false;
    }
  }

  for (long i = params.d; i-- > 0;) {
    const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
    const long n_i = CodewordLengthAtLevel(params, i);
    if (mu < 0 || mu >= n_i)
      return false;

    FieldElement x1, x2;
    FoldingPoints(x1, x2, params, i, mu);
    const FieldElement &y1 = oracles.pi[static_cast<std::size_t>(i + 1)][mu];
    const FieldElement &y2 =
        oracles.pi[static_cast<std::size_t>(i + 1)][mu + n_i];
    const FieldElement expected = EvalLineAt(
        challenges.alphas[static_cast<std::size_t>(i)], x1, y1, x2, y2);
    if (expected != oracles.pi[static_cast<std::size_t>(i)][mu])
      return false;
  }

  return IsCodewordC0(oracles.pi[0], params);
}

MerkleRoot MerkleCommitOracle(const Oracle &oracle) {
  const long leaf_count = oracle.length();
  if (leaf_count < 0)
    LogicError("MerkleCommitOracle: invalid leaf count");
  if (leaf_count == 0)
    return HashRootWithCount(0, MerkleRootRaw({}));

  std::vector<Bytes> leaf_hashes;
  leaf_hashes.reserve(static_cast<std::size_t>(leaf_count));
  for (long i = 0; i < leaf_count; ++i) {
    leaf_hashes.push_back(HashLeaf(i, oracle[i]));
  }
  const Bytes raw_root = MerkleRootRaw(leaf_hashes);
  return HashRootWithCount(leaf_count, raw_root);
}

MerkleOpening MerkleOpenOracle(const Oracle &oracle, long index) {
  const long leaf_count = oracle.length();
  if (index < 0 || index >= leaf_count) {
    LogicError("MerkleOpenOracle: index out of range");
  }

  std::vector<Bytes> level;
  level.reserve(static_cast<std::size_t>(leaf_count));
  for (long i = 0; i < leaf_count; ++i) {
    level.push_back(HashLeaf(i, oracle[i]));
  }

  MerkleAuthPath path;
  const std::size_t height = ExpectedMerkleHeight(leaf_count);
  path.sibling_hashes.reserve(height);
  path.sibling_is_left.reserve(height);

  long idx = index;
  while (level.size() > 1) {
    if (level.size() % 2 == 1)
      level.push_back(level.back());

    const long sibling = (idx % 2 == 0) ? (idx + 1) : (idx - 1);
    path.sibling_hashes.push_back(level[static_cast<std::size_t>(sibling)]);
    path.sibling_is_left.push_back(idx % 2 == 1);

    std::vector<Bytes> next;
    next.reserve(level.size() / 2);
    for (std::size_t i = 0; i < level.size(); i += 2) {
      next.push_back(HashNode(level[i], level[i + 1]));
    }
    idx /= 2;
    level = std::move(next);
  }

  MerkleOpening opening;
  opening.index = index;
  opening.value = oracle[index];
  opening.auth_path = std::move(path);
  return opening;
}

bool MerkleVerifyOpening(const MerkleRoot &root, long leaf_count,
                         const MerkleOpening &opening) {
  if (leaf_count < 0)
    return false;
  if (opening.index < 0 || opening.index >= leaf_count)
    return false;

  const std::size_t expected_height = ExpectedMerkleHeight(leaf_count);
  if (opening.auth_path.sibling_hashes.size() != expected_height)
    return false;
  if (opening.auth_path.sibling_is_left.size() != expected_height)
    return false;

  Bytes cur = HashLeaf(opening.index, opening.value);
  for (std::size_t i = 0; i < expected_height; ++i) {
    const Bytes &sib = opening.auth_path.sibling_hashes[i];
    if (opening.auth_path.sibling_is_left[i]) {
      cur = HashNode(sib, cur);
    } else {
      cur = HashNode(cur, sib);
    }
  }

  const Bytes expected_root = HashRootWithCount(leaf_count, cur);
  return expected_root == root;
}

MerkleTree MerkleTree::Build(const Oracle &oracle) {
  const long leaf_count = oracle.length();
  if (leaf_count < 0)
    LogicError("MerkleTree::Build: invalid leaf count");

  MerkleTree t;
  t.leaf_count_ = leaf_count;
  t.levels_.clear();
  t.raw_root_.clear();

  if (leaf_count == 0) {
    t.raw_root_ = HashWithPrefix(static_cast<Byte>(0x04), {});
    return t;
  }

  std::vector<Bytes> level;
  level.reserve(static_cast<std::size_t>(leaf_count));
  for (long i = 0; i < leaf_count; ++i) {
    level.push_back(HashLeaf(i, oracle[i]));
  }

  while (level.size() > 1) {
    if (level.size() % 2 == 1)
      level.push_back(level.back());
    t.levels_.push_back(level);

    std::vector<Bytes> next;
    next.reserve(level.size() / 2);
    for (std::size_t i = 0; i < level.size(); i += 2) {
      next.push_back(HashNode(level[i], level[i + 1]));
    }
    level = std::move(next);
  }

  t.raw_root_ = level[0];
  return t;
}

MerkleRoot MerkleTree::Root() const {
  return HashRootWithCount(leaf_count_, raw_root_);
}

MerkleOpening MerkleTree::Open(const Oracle &oracle, long index) const {
  if (oracle.length() != leaf_count_) {
    LogicError("MerkleTree::Open: oracle length mismatch");
  }
  if (index < 0 || index >= leaf_count_) {
    LogicError("MerkleTree::Open: index out of range");
  }

  MerkleAuthPath path;
  const std::size_t height = ExpectedMerkleHeight(leaf_count_);
  path.sibling_hashes.reserve(height);
  path.sibling_is_left.reserve(height);

  long idx = index;
  for (std::size_t h = 0; h < height; ++h) {
    const std::vector<Bytes> &level = levels_[h];
    const long sibling = (idx % 2 == 0) ? (idx + 1) : (idx - 1);
    path.sibling_hashes.push_back(level[static_cast<std::size_t>(sibling)]);
    path.sibling_is_left.push_back(idx % 2 == 1);
    idx /= 2;
  }

  MerkleOpening opening;
  opening.index = index;
  opening.value = oracle[index];
  opening.auth_path = std::move(path);
  return opening;
}

IOPPChallenges
FiatShamirDeriveChallenges(FiatShamirTranscript &transcript,
                           const IOPPMerkleCommitments &commitments,
                           const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (static_cast<long>(commitments.roots_by_level.size()) != params.d + 1) {
    LogicError("FiatShamirDeriveChallenges: roots_by_level has wrong size");
  }

  IOPPChallenges out;
  out.alphas.resize(static_cast<std::size_t>(params.d));
  if (params.d == 0)
    return out;

  transcript.AbsorbBytes(
      commitments.roots_by_level[static_cast<std::size_t>(params.d)]);
  for (long i = params.d; i-- > 0;) {
    out.alphas[static_cast<std::size_t>(i)] =
        transcript.ChallengeFieldElement("alpha/" + std::to_string(i));
    transcript.AbsorbBytes(
        commitments.roots_by_level[static_cast<std::size_t>(i)]);
  }

  return out;
}

std::vector<IOPPQueryPlan>
FiatShamirDeriveQueryPlans(FiatShamirTranscript &transcript, long num_queries,
                           const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (num_queries < 0)
    LogicError("FiatShamirDeriveQueryPlans: negative count");

  std::vector<IOPPQueryPlan> plans;
  plans.reserve(static_cast<std::size_t>(num_queries));

  if (params.d == 0) {
    for (long q = 0; q < num_queries; ++q)
      plans.push_back(MakeQueryPlan(0, params));
    return plans;
  }

  const long n_last = CodewordLengthAtLevel(params, params.d - 1);
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    plans.push_back(MakeQueryPlan(mu, params));
  }

  return plans;
}

bool VerifyQueryFromMerkleOpenings(const IOPPQueryPlan &plan,
                                   const IOPPChallenges &challenges,
                                   const IOPPQueryMerkleOpenings &openings,
                                   const IOPPMerkleCommitments &commitments,
                                   const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (static_cast<long>(commitments.roots_by_level.size()) != params.d + 1) {
    return false;
  }
  if (static_cast<long>(challenges.alphas.size()) != params.d)
    return false;
  if (static_cast<long>(plan.mu_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(openings.left.size()) != params.d)
    return false;
  if (static_cast<long>(openings.right.size()) != params.d)
    return false;
  if (static_cast<long>(openings.folded.size()) != params.d)
    return false;

  if (MerkleCommitOracle(openings.pi0_full) != commitments.roots_by_level[0]) {
    return false;
  }

  for (long i = params.d; i-- > 0;) {
    const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
    const long n_i = CodewordLengthAtLevel(params, i);
    if (mu < 0 || mu >= n_i)
      return false;

    const MerkleOpening &left = openings.left[static_cast<std::size_t>(i)];
    const MerkleOpening &right = openings.right[static_cast<std::size_t>(i)];
    const MerkleOpening &folded = openings.folded[static_cast<std::size_t>(i)];

    if (left.index != mu)
      return false;
    if (right.index != mu + n_i)
      return false;
    if (folded.index != mu)
      return false;

    const long n_ip1 = 2 * n_i;
    if (!MerkleVerifyOpening(
            commitments.roots_by_level[static_cast<std::size_t>(i + 1)], n_ip1,
            left)) {
      return false;
    }
    if (!MerkleVerifyOpening(
            commitments.roots_by_level[static_cast<std::size_t>(i + 1)], n_ip1,
            right)) {
      return false;
    }
    if (!MerkleVerifyOpening(
            commitments.roots_by_level[static_cast<std::size_t>(i)], n_i,
            folded)) {
      return false;
    }

    FieldElement x1, x2;
    FoldingPoints(x1, x2, params, i, mu);
    const FieldElement expected =
        EvalLineAt(challenges.alphas[static_cast<std::size_t>(i)], x1,
                   left.value, x2, right.value);
    if (expected != folded.value)
      return false;
  }

  return IsCodewordC0(openings.pi0_full, params);
}

bool IsCodewordC0(const Oracle &pi0, const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  const long n0 = CodewordLengthAtLevel(params, 0);
  if (pi0.length() != n0) return false;

  vec_ZZ_pE msg0;
  return DecodeC0(msg0, pi0, params);
}

bool DecodeC0(vec_ZZ_pE &msg0_out, const Oracle &pi0,
              const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  const long n0 = CodewordLengthAtLevel(params, 0);
  if (pi0.length() != n0) return false;

  if (params.k0 == 1) {
    ZZ_pE inv_g;
    long pivot = -1;
    for (long j = 0; j < n0; ++j) {
      if (TryInvertUnit(inv_g, params.G0[0][j])) {
        pivot = j;
        break;
      }
    }
    if (pivot < 0) return false;

    msg0_out.SetLength(1);
    msg0_out[0] = pi0[pivot] * inv_g;

    vec_ZZ_pE rec;
    mul(rec, msg0_out, params.G0);
    return rec == pi0;
  }

  mat_ZZ_pE a;
  a.SetDims(n0, params.k0 + 1);

  for (long row = 0; row < n0; ++row) {
    for (long col = 0; col < params.k0; ++col) {
      a[row][col] = params.G0[col][row];
    }
    a[row][params.k0] = pi0[row];
  }

  if (!SolveLinearSystemRref(msg0_out, a)) return false;

  vec_ZZ_pE rec;
  mul(rec, msg0_out, params.G0);
  return rec == pi0;
}

} // namespace basefold
