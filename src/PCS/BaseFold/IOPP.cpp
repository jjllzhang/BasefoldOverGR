#include "PCS/BaseFold/IOPP.hpp"
#include "PCS/Common/Profile.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "GaloisRing/Inverse.hpp"

using NTL::coeff;
using NTL::LogicError;
using NTL::mat_ZZ_pE;
using NTL::mul;
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

bool IsUnitInBaseRing(const NTL::ZZ_p &a) {
  if (a == 0) return false;
  NTL::ZZ g;
  NTL::GCD(g, NTL::rep(a), NTL::ZZ_p::modulus());
  return g == 1;
}

bool BaseModulusIsPrime() {
  static thread_local ZZ cached_modulus;
  static thread_local bool cached_is_prime = false;
  static thread_local bool has_cache = false;

  const ZZ &modulus = NTL::ZZ_p::modulus();
  if (!has_cache || modulus != cached_modulus) {
    cached_modulus = modulus;
    cached_is_prime = NTL::ProbPrime(modulus);
    has_cache = true;
  }
  return cached_is_prime;
}

bool IsUnitInCurrentZZpEContext(const ZZ_pE &a) {
  if (a == 0) return false;
  const long r = ZZ_pE::degree();
  if (r <= 0) LogicError("IsUnitInCurrentZZpEContext: invalid extension degree");
  if (r == 1) {
    return IsUnitInBaseRing(coeff(rep(a), 0));
  }
  const NTL::ZZ_pX &poly = rep(a);
  for (long i = 0; i < r; ++i) {
    if (IsUnitInBaseRing(coeff(poly, i))) return true;
  }
  return false;
}

bool TryInvertUnitInCurrentZZpEContext(ZZ_pE &inv_out, const ZZ_pE &a) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->try_invert_unit_ns : nullptr,
                    prof ? &prof->try_invert_unit_calls : nullptr);

  if (a == 0) return false;
  const long r = ZZ_pE::degree();
  if (r <= 0) {
    LogicError("TryInvertUnitInCurrentZZpEContext: invalid extension degree");
  }

  // Context-aware single-entry cache: inversions are often repeated (e.g. when
  // diag_T is constant, all folding denominators match).
  static thread_local bool cache_valid = false;
  static thread_local ZZ cache_modulus;
  static thread_local long cache_r = 0;
  static thread_local ZZ_pE cache_a;
  static thread_local ZZ_pE cache_inv;

  const ZZ modulus = NTL::ZZ_p::modulus();
  if (cache_valid && cache_r == r && modulus == cache_modulus && a == cache_a) {
    if (prof != nullptr) ++prof->try_invert_unit_cache_hits;
    inv_out = cache_inv;
    return true;
  }

  ZZ_pE inv_candidate;
  bool ok = false;

  if (BaseModulusIsPrime()) {
    try {
      inv_candidate = NTL::inv(a);
      ok = true;
    } catch (...) {
      ok = false;
    }
    if (!ok) return false;
  } else {
    bool is_unit = false;
    if (prof != nullptr) {
      ScopedTimer is_unit_timer(&prof->is_unit_ns, &prof->is_unit_calls);
      is_unit = IsUnitInCurrentZZpEContext(a);
    } else {
      is_unit = IsUnitInCurrentZZpEContext(a);
    }
    if (!is_unit) return false;

    if (r == 1) {
      if (prof != nullptr) {
        ScopedTimer inv_timer(&prof->inv_fallback_ns,
                              &prof->inv_fallback_calls);
        const NTL::ZZ a_rep = NTL::rep(coeff(rep(a), 0));
        NTL::ZZ inv_rep;
        if (NTL::InvModStatus(inv_rep, a_rep, NTL::ZZ_p::modulus()) != 0)
          return false;
        NTL::ZZ_p inv_base;
        NTL::conv(inv_base, inv_rep);
        NTL::ZZ_pX poly;
        NTL::clear(poly);
        NTL::SetCoeff(poly, 0, inv_base);
        NTL::conv(inv_candidate, poly);
        ok = true;
      } else {
        const NTL::ZZ a_rep = NTL::rep(coeff(rep(a), 0));
        NTL::ZZ inv_rep;
        if (NTL::InvModStatus(inv_rep, a_rep, NTL::ZZ_p::modulus()) != 0)
          return false;
        NTL::ZZ_p inv_base;
        NTL::conv(inv_base, inv_rep);
        NTL::ZZ_pX poly;
        NTL::clear(poly);
        NTL::SetCoeff(poly, 0, inv_base);
        NTL::conv(inv_candidate, poly);
        ok = true;
      }
    } else {
      if (prof != nullptr) {
        ScopedTimer inv_timer(&prof->inv_fallback_ns,
                              &prof->inv_fallback_calls);
        inv_candidate = Inv(a, r);
      } else {
        inv_candidate = Inv(a, r);
      }
      if (inv_candidate == 0) return false;
      ZZ_pE one;
      NTL::set(one);
      ok = (a * inv_candidate == one);
    }

    if (!ok) return false;
  }

  inv_out = inv_candidate;

  cache_valid = true;
  cache_modulus = modulus;
  cache_r = r;
  cache_a = a;
  cache_inv = inv_candidate;

  return true;
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
        if (TryInvertUnitInCurrentZZpEContext(inv_pivot, aug[r][c])) {
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
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->eval_line_at_ns : nullptr,
                    prof ? &prof->eval_line_at_calls : nullptr);

  const FieldElement denom = x2 - x1;
  if (denom == 0)
    LogicError("EvalLineAt: x1 must not equal x2");
  ZZ_pE inv_denom;
  if (!TryInvertUnitInCurrentZZpEContext(inv_denom, denom)) {
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

  oracles.oracles_by_level.resize(static_cast<std::size_t>(params.d + 1));
  oracles.oracles_by_level[static_cast<std::size_t>(params.d)] = pi_d;

  for (long i = params.d; i-- > 0;) {
    ProverCommitRound(oracles.oracles_by_level[static_cast<std::size_t>(i)],
                      oracles.oracles_by_level[static_cast<std::size_t>(i + 1)],
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
  if (static_cast<long>(openings.upper_left_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(openings.upper_right_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(openings.folded_by_level.size()) != params.d)
    return false;

  for (long i = params.d; i-- > 0;) {
    const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
    const long n_i = CodewordLengthAtLevel(params, i);
    if (mu < 0 || mu >= n_i)
      return false;

    FieldElement x1, x2;
    FoldingPoints(x1, x2, params, i, mu);
    const FieldElement &y1 =
        openings.upper_left_by_level[static_cast<std::size_t>(i)];
    const FieldElement &y2 =
        openings.upper_right_by_level[static_cast<std::size_t>(i)];
    const FieldElement expected = EvalLineAt(
        challenges.alphas[static_cast<std::size_t>(i)], x1, y1, x2, y2);
    if (expected != openings.folded_by_level[static_cast<std::size_t>(i)])
      return false;
  }

  return IsCodewordC0(openings.pi0_codeword, params);
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
  if (static_cast<long>(oracles.oracles_by_level.size()) != params.d + 1)
    return false;

  for (long i = 0; i <= params.d; ++i) {
    if (oracles.oracles_by_level[static_cast<std::size_t>(i)].length() !=
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
    const FieldElement &y1 =
        oracles.oracles_by_level[static_cast<std::size_t>(i + 1)][mu];
    const FieldElement &y2 =
        oracles.oracles_by_level[static_cast<std::size_t>(i + 1)][mu + n_i];
    const FieldElement expected = EvalLineAt(
        challenges.alphas[static_cast<std::size_t>(i)], x1, y1, x2, y2);
    if (expected != oracles.oracles_by_level[static_cast<std::size_t>(i)][mu])
      return false;
  }

  return IsCodewordC0(oracles.oracles_by_level[0], params);
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

  {
    const MerkleRoot &root_d =
        commitments.roots_by_level[static_cast<std::size_t>(params.d)];
    transcript.AbsorbBytes(root_d.data(), root_d.size());
  }
  for (long i = params.d; i-- > 0;) {
    out.alphas[static_cast<std::size_t>(i)] =
        transcript.ChallengeFieldElement("alpha/" + std::to_string(i));
    const MerkleRoot &root_i =
        commitments.roots_by_level[static_cast<std::size_t>(i)];
    transcript.AbsorbBytes(root_i.data(), root_i.size());
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
      if (TryInvertUnitInCurrentZZpEContext(inv_g, params.G0[0][j])) {
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
