#include "PCS/Common/Sumcheck.hpp"
#include "PCS/Common/NtlParallel.hpp"
#include "PCS/Common/Profile.hpp"

#include <NTL/ZZ.h>

#include <cstddef>
#include <vector>

#include "PCS/Common/Multilinear.hpp"

using NTL::LogicError;
using NTL::vec_ZZ_pE;
using NTL::ZZ_pE;

namespace basefold {
namespace {

constexpr long kProductSumcheckCurrentPolyParallelThreshold = 512;
constexpr long kProductSumcheckReceiveChallengeParallelThreshold = 512;

bool IsPowerOfTwoLong(long n) { return n > 0 && (n & (n - 1)) == 0; }

long Log2ExactPowerOfTwoLong(long n) {
  if (!IsPowerOfTwoLong(n))
    LogicError("Log2ExactPowerOfTwoLong: not a power of two");
  long d = 0;
  while (n > 1) {
    n >>= 1;
    ++d;
  }
  return d;
}

ZZ_pE One() {
  ZZ_pE one;
  NTL::set(one);
  return one;
}

QuadraticPoly ZeroQuadraticPoly() { return {ZZ_pE(0), ZZ_pE(0), ZZ_pE(0)}; }

// In-place subset-sum (zeta) transform:
//   out[mask] = Σ_{S ⊆ mask} in[S]
vec_ZZ_pE BooleanEvalTableFromMonomialCoeffs(const vec_ZZ_pE &coeffs, long k) {
  const long n = coeffs.length();
  if (k < 0)
    LogicError("BooleanEvalTableFromMonomialCoeffs: negative k");
  if (n != (1L << k))
    LogicError("BooleanEvalTableFromMonomialCoeffs: length mismatch");

  vec_ZZ_pE eval = coeffs;
  for (long bit = 0; bit < k; ++bit) {
    const long step = 1L << bit;
    for (long mask = 0; mask < n; ++mask) {
      if (mask & step) {
        eval[mask] += eval[mask ^ step];
      }
    }
  }
  return eval;
}

long ValidateEvalTableLengthOrThrow(const vec_ZZ_pE &table,
                                    const char *func_name) {
  const long n = table.length();
  if (!IsPowerOfTwoLong(n)) {
    LogicError(
        (std::string(func_name) + ": table length must be a power of two")
            .c_str());
  }
  return Log2ExactPowerOfTwoLong(n);
}

void ValidateSameEvalTableShapeOrThrow(const vec_ZZ_pE &f_table,
                                       const vec_ZZ_pE &g_table,
                                       const char *func_name) {
  const long f_dim = ValidateEvalTableLengthOrThrow(f_table, func_name);
  const long g_dim = ValidateEvalTableLengthOrThrow(g_table, func_name);
  if (f_dim != g_dim) {
    LogicError(
        (std::string(func_name) + ": f and g table dimensions must match")
            .c_str());
  }
}

bool CheckQuadraticSumcheckChain(const std::vector<QuadraticPoly> &h_by_level,
                                 const std::vector<FieldElement> &r,
                                 const FieldElement &initial_claim) {
  const long d = static_cast<long>(h_by_level.size());
  if (d < 0)
    return false;
  if (static_cast<long>(r.size()) != d)
    return false;
  if (d == 0)
    return true;

  const ZZ_pE one = One();
  const ZZ_pE zero = ZZ_pE(0);

  const QuadraticPoly &h_d = h_by_level[static_cast<std::size_t>(d - 1)];
  if (h_d.Eval(zero) + h_d.Eval(one) != initial_claim)
    return false;

  for (long k = 1; k < d; ++k) {
    const QuadraticPoly &h_k = h_by_level[static_cast<std::size_t>(k - 1)];
    const QuadraticPoly &h_kp1 = h_by_level[static_cast<std::size_t>(k)];
    const ZZ_pE lhs = h_k.Eval(zero) + h_k.Eval(one);
    const ZZ_pE rhs = h_kp1.Eval(r[static_cast<std::size_t>(k)]);
    if (lhs != rhs)
      return false;
  }

  return true;
}

std::vector<vec_ZZ_pE> PrecomputePrefixEqByVars(const std::vector<ZZ_pE> &z) {
  const long d = static_cast<long>(z.size());
  if (d < 0)
    LogicError("PrecomputePrefixEqByVars: negative dimension");

  std::vector<vec_ZZ_pE> out;
  if (d == 0)
    return out;

  out.resize(static_cast<std::size_t>(d));
  const ZZ_pE one = One();

  out[0].SetLength(1);
  out[0][0] = one;

  // out[t] uses z[0..t-1] and has length 2^t.
  for (long t = 1; t < d; ++t) {
    const ZZ_pE z_var = z[static_cast<std::size_t>(t - 1)];
    const ZZ_pE f0 = one - z_var; // factor(var, 0)
    const ZZ_pE f1 = z_var;       // factor(var, 1)

    const vec_ZZ_pE &prev = out[static_cast<std::size_t>(t - 1)];
    const long old = prev.length();
    out[static_cast<std::size_t>(t)].SetLength(2 * old);
    for (long mask = 0; mask < old; ++mask) {
      const ZZ_pE base = prev[mask];
      out[static_cast<std::size_t>(t)][mask] = base * f0;
      out[static_cast<std::size_t>(t)][mask + old] = base * f1;
    }
  }

  return out;
}

} // namespace

FieldElement QuadraticPoly::Eval(const FieldElement &x) const {
  return a0 + a1 * x + a2 * x * x;
}

ProductSumcheckProver::ProductSumcheckProver(const FieldVec &f_eval_table,
                                             const FieldVec &g_eval_table) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->product_sumcheck_init_ns : nullptr,
                    prof ? &prof->product_sumcheck_init_calls : nullptr);
  ValidateSameEvalTableShapeOrThrow(f_eval_table, g_eval_table,
                                    "ProductSumcheckProver");
  d_ = ValidateEvalTableLengthOrThrow(f_eval_table, "ProductSumcheckProver");
  cur_k_ = d_;
  f_eval_table_ = f_eval_table;
  g_eval_table_ = g_eval_table;
}

ProductSumcheckProver
ProductSumcheckProver::FromMonomialCoeffs(const FieldVec &f_coeffs,
                                          const FieldVec &g_coeffs) {
  if (f_coeffs.length() != g_coeffs.length()) {
    LogicError("ProductSumcheckProver::FromMonomialCoeffs: f and g lengths "
               "must match");
  }
  const long d = ValidateEvalTableLengthOrThrow(
      f_coeffs, "ProductSumcheckProver::FromMonomialCoeffs");
  return ProductSumcheckProver(BooleanEvalTableFromMonomialCoeffs(f_coeffs, d),
                               BooleanEvalTableFromMonomialCoeffs(g_coeffs, d));
}

SumcheckMonomialPrecomputation
BuildSumcheckMonomialPrecomputation(const FieldVec &f_coeffs) {
  const long d = ValidateEvalTableLengthOrThrow(
      f_coeffs, "BuildSumcheckMonomialPrecomputation");
  SumcheckMonomialPrecomputation out;
  out.valid = true;
  out.d = d;
  out.f_eval_table = BooleanEvalTableFromMonomialCoeffs(f_coeffs, d);
  return out;
}

QuadraticPoly ProductSumcheckProver::CurrentPolynomial() const {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->product_sumcheck_current_poly_ns : nullptr,
                    prof ? &prof->product_sumcheck_current_poly_calls
                         : nullptr);
  if (cur_k_ <= 0) {
    LogicError(
        "ProductSumcheckProver::CurrentPolynomial: no remaining variables");
  }

  const long k = cur_k_;
  const long n = f_eval_table_.length();
  if (n != (1L << k) || g_eval_table_.length() != n) {
    LogicError(
        "ProductSumcheckProver::CurrentPolynomial: internal length mismatch");
  }

  const long half = 1L << (k - 1);
  QuadraticPoly out = ZeroQuadraticPoly();
  bool accumulated_in_parallel = false;
#if defined(BASEFOLD_USE_OPENMP)
  const int threads_to_use = pcs_common_internal::ChooseElementParallelThreads(
      half, kProductSumcheckCurrentPolyParallelThreshold);
  if (threads_to_use >= 2) {
    const pcs_common_internal::NtlThreadContextSnapshot ntl_ctx =
        pcs_common_internal::CaptureNtlThreadContextSnapshot();
    std::vector<QuadraticPoly> partials(
        static_cast<std::size_t>(threads_to_use), ZeroQuadraticPoly());

#pragma omp parallel num_threads(threads_to_use)
    {
      pcs_common_internal::InitNtlThreadContext(ntl_ctx);
      const int tid = omp_get_thread_num();
      QuadraticPoly &local = partials[static_cast<std::size_t>(tid)];

#pragma omp for schedule(static)
      for (long mask = 0; mask < half; ++mask) {
        const ZZ_pE f0 = f_eval_table_[mask];
        const ZZ_pE f1 = f_eval_table_[mask + half];
        const ZZ_pE delta_f = f1 - f0;

        const ZZ_pE g0 = g_eval_table_[mask];
        const ZZ_pE g1 = g_eval_table_[mask + half];
        const ZZ_pE delta_g = g1 - g0;

        local.a0 += f0 * g0;
        local.a1 += f0 * delta_g + delta_f * g0;
        local.a2 += delta_f * delta_g;
      }
    }

    for (const QuadraticPoly &partial : partials) {
      out.a0 += partial.a0;
      out.a1 += partial.a1;
      out.a2 += partial.a2;
    }
    accumulated_in_parallel = true;
  }
#endif
  if (!accumulated_in_parallel) {
    for (long mask = 0; mask < half; ++mask) {
      const ZZ_pE f0 = f_eval_table_[mask];
      const ZZ_pE f1 = f_eval_table_[mask + half];
      const ZZ_pE delta_f = f1 - f0;

      const ZZ_pE g0 = g_eval_table_[mask];
      const ZZ_pE g1 = g_eval_table_[mask + half];
      const ZZ_pE delta_g = g1 - g0;

      out.a0 += f0 * g0;
      out.a1 += f0 * delta_g + delta_f * g0;
      out.a2 += delta_f * delta_g;
    }
  }

  return out;
}

void ProductSumcheckProver::ReceiveChallenge(const FieldElement &r_kminus1) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(
      prof ? &prof->product_sumcheck_receive_challenge_ns : nullptr,
      prof ? &prof->product_sumcheck_receive_challenge_calls : nullptr);
  if (cur_k_ <= 0) {
    LogicError(
        "ProductSumcheckProver::ReceiveChallenge: no remaining variables");
  }

  const long k = cur_k_;
  const long n = f_eval_table_.length();
  if (n != (1L << k) || g_eval_table_.length() != n) {
    LogicError(
        "ProductSumcheckProver::ReceiveChallenge: internal length mismatch");
  }

  const long half = n / 2;
  pcs_common_internal::ForEachIndexMaybeParallel(
      0, half, kProductSumcheckReceiveChallengeParallelThreshold, [&](long i) {
        const ZZ_pE f0 = f_eval_table_[i];
        const ZZ_pE f1 = f_eval_table_[i + half];
        f_eval_table_[i] = f0 + (f1 - f0) * r_kminus1;

        const ZZ_pE g0 = g_eval_table_[i];
        const ZZ_pE g1 = g_eval_table_[i + half];
        g_eval_table_[i] = g0 + (g1 - g0) * r_kminus1;
      });
  f_eval_table_.SetLength(half);
  g_eval_table_.SetLength(half);
  --cur_k_;
}

SumcheckProver::SumcheckProver(const FieldVec &f_coeffs,
                               const std::vector<FieldElement> &z)
    : SumcheckProver(BuildSumcheckMonomialPrecomputation(f_coeffs), z) {}

SumcheckProver::SumcheckProver(
    const SumcheckMonomialPrecomputation &precomputation,
    const std::vector<FieldElement> &z)
    : z_(z) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->sumcheck_init_ns : nullptr,
                    prof ? &prof->sumcheck_init_calls : nullptr);

  if (!precomputation.valid) {
    LogicError("SumcheckProver: precomputation must be valid");
  }

  d_ = ValidateEvalTableLengthOrThrow(precomputation.f_eval_table,
                                      "SumcheckProver");
  if (precomputation.d != d_) {
    LogicError("SumcheckProver: precomputation dimension mismatch");
  }
  if (z_.size() != static_cast<std::size_t>(d_))
    LogicError("SumcheckProver: z dimension mismatch");

  cur_k_ = d_;
  f_eval_table_ = precomputation.f_eval_table;
  prefix_eq_by_vars_ = PrecomputePrefixEqByVars(z_);
  suffix_eq_prod_ = One();
}

QuadraticPoly SumcheckProver::CurrentPolynomial() const {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->sumcheck_current_poly_ns : nullptr,
                    prof ? &prof->sumcheck_current_poly_calls : nullptr);

  if (cur_k_ <= 0)
    LogicError("SumcheckProver::CurrentPolynomial: no remaining variables");

  const long k = cur_k_;
  const long n = f_eval_table_.length();
  if (n != (1L << k))
    LogicError("SumcheckProver::CurrentPolynomial: internal length mismatch");

  const long half = 1L << (k - 1);

  if (static_cast<long>(prefix_eq_by_vars_.size()) != d_)
    LogicError("SumcheckProver::CurrentPolynomial: prefix table size mismatch");
  const vec_ZZ_pE &prefix = prefix_eq_by_vars_[static_cast<std::size_t>(k - 1)];
  if (prefix.length() != half)
    LogicError("SumcheckProver::CurrentPolynomial: prefix size mismatch");

  const ZZ_pE one = One();
  const ZZ_pE z_k = z_[static_cast<std::size_t>(k - 1)];
  const ZZ_pE factor0 = one - z_k; // eq-factor at X_k=0
  const ZZ_pE factor1 = z_k;       // eq-factor at X_k=1
  const ZZ_pE delta_factor = factor1 - factor0;

  QuadraticPoly out;
  out.a0 = ZZ_pE(0);
  out.a1 = ZZ_pE(0);
  out.a2 = ZZ_pE(0);

  for (long mask = 0; mask < half; ++mask) {
    const ZZ_pE common = prefix[mask] * suffix_eq_prod_;

    const ZZ_pE eq0 = common * factor0;
    const ZZ_pE delta_eq = common * delta_factor;

    const ZZ_pE f0 = f_eval_table_[mask];
    const ZZ_pE f1 = f_eval_table_[mask + half];
    const ZZ_pE delta_f = f1 - f0;

    out.a0 += f0 * eq0;
    out.a1 += f0 * delta_eq + delta_f * eq0;
    out.a2 += delta_f * delta_eq;
  }

  return out;
}

void SumcheckProver::ReceiveChallenge(const FieldElement &r_kminus1) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->sumcheck_receive_challenge_ns : nullptr,
                    prof ? &prof->sumcheck_receive_challenge_calls : nullptr);

  if (cur_k_ <= 0)
    LogicError("SumcheckProver::ReceiveChallenge: no remaining variables");

  const long k = cur_k_;
  const long n = f_eval_table_.length();
  if (n != (1L << k))
    LogicError("SumcheckProver::ReceiveChallenge: internal length mismatch");

  suffix_eq_prod_ *= EqFactor(z_[static_cast<std::size_t>(k - 1)], r_kminus1);

  const long half = n / 2;
  for (long i = 0; i < half; ++i) {
    const ZZ_pE f0 = f_eval_table_[i];
    const ZZ_pE f1 = f_eval_table_[i + half];
    f_eval_table_[i] = f0 + (f1 - f0) * r_kminus1;
  }
  f_eval_table_.SetLength(half);
  --cur_k_;
}

bool CheckSumcheckRelations(const std::vector<QuadraticPoly> &h_by_level,
                            const std::vector<FieldElement> &r,
                            const FieldElement &claimed_y) {
  return CheckQuadraticSumcheckChain(h_by_level, r, claimed_y);
}

bool CheckProductSumcheckChain(const FieldElement &initial_claim,
                               const std::vector<QuadraticPoly> &h_by_level,
                               const std::vector<FieldElement> &r) {
  return CheckQuadraticSumcheckChain(h_by_level, r, initial_claim);
}

} // namespace basefold
