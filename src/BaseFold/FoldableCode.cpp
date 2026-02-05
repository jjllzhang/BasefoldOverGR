#include "BaseFold/FoldableCode.hpp"

#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_pE.h>

using NTL::ZZ;
using NTL::LogicError;
using NTL::coeff;
using NTL::mul;
using NTL::rep;
using NTL::set;
using NTL::vec_ZZ_pE;
using NTL::ZZ_pE;
using NTL::ZZ_pX;

namespace basefold {
namespace {

long Pow2(long e) {
  if (e < 0) {
    LogicError("Pow2: negative exponent");
  }
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2: exponent too large for long");
  }
  return 1L << e;
}

bool IsUnitByReductionModP(const ZZ_pE &a, const ZZ &p) {
  if (p <= 1) {
    return a != 0;
  }
  if (a == 0) return false;

  const ZZ_pX poly = rep(a);
  const long r = ZZ_pE::degree();
  for (long i = 0; i < r; ++i) {
    ZZ c = rep(coeff(poly, i));
    c %= p;
    if (c != 0) return true;
  }
  return false;
}

void ValidateParams(const FoldableCodeParams &params) {
  static thread_local bool cache_valid = false;
  static thread_local const FoldableCodeParams *cache_params = nullptr;
  static thread_local ZZ cache_modulus;
  static thread_local long cache_degree = 0;
  static thread_local long cache_c = 0;
  static thread_local long cache_k0 = 0;
  static thread_local long cache_d = 0;
  static thread_local long cache_g0_rows = 0;
  static thread_local long cache_g0_cols = 0;
  static thread_local ZZ cache_p;

  if (cache_valid && cache_params == &params &&
      cache_modulus == NTL::ZZ_p::modulus() &&
      cache_degree == ZZ_pE::degree() && cache_c == params.c &&
      cache_k0 == params.k0 && cache_d == params.d &&
      cache_g0_rows == params.G0.NumRows() &&
      cache_g0_cols == params.G0.NumCols() && cache_p == params.p &&
      static_cast<long>(params.diag_T.size()) == params.d) {
    return;
  }

  if (params.c <= 0) LogicError("EncodeFoldable: c must be positive");
  if (params.k0 <= 0) LogicError("EncodeFoldable: k0 must be positive");
  if (params.d < 0) LogicError("EncodeFoldable: d must be non-negative");

  const long n0 = params.c * params.k0;
  if (params.G0.NumRows() != params.k0 || params.G0.NumCols() != n0) {
    LogicError("EncodeFoldable: G0 has wrong dimensions");
  }

  if (params.zeta == 0) LogicError("EncodeFoldable: zeta must be non-zero");
  if (params.zeta == 1) LogicError("EncodeFoldable: zeta must not equal 1");
  if (params.p > 1) {
    if (!IsUnitByReductionModP(params.zeta, params.p)) {
      LogicError("EncodeFoldable: zeta must be a unit");
    }
    ZZ_pE one;
    set(one);
    if (!IsUnitByReductionModP(one - params.zeta, params.p)) {
      LogicError("EncodeFoldable: 1-zeta must be a unit");
    }
  }

  if (static_cast<long>(params.diag_T.size()) != params.d) {
    LogicError("EncodeFoldable: diag_T must have size d");
  }

  for (long i = 0; i < params.d; ++i) {
    const long ni = params.c * params.k0 * Pow2(i);
    if (params.diag_T[static_cast<std::size_t>(i)].length() != ni) {
      LogicError("EncodeFoldable: diag_T[i] has wrong length");
    }
    for (long j = 0; j < ni; ++j) {
      const ZZ_pE &tij = params.diag_T[static_cast<std::size_t>(i)][j];
      if (tij == 0) {
        LogicError("EncodeFoldable: diag_T[i] entries must be non-zero");
      }
      if (params.p > 1 && !IsUnitByReductionModP(tij, params.p)) {
        LogicError("EncodeFoldable: diag_T[i] entries must be units");
      }
    }
  }

  cache_valid = true;
  cache_params = &params;
  cache_modulus = NTL::ZZ_p::modulus();
  cache_degree = ZZ_pE::degree();
  cache_c = params.c;
  cache_k0 = params.k0;
  cache_d = params.d;
  cache_g0_rows = params.G0.NumRows();
  cache_g0_cols = params.G0.NumCols();
  cache_p = params.p;
}

vec_ZZ_pE Encode0(const vec_ZZ_pE &msg, const FoldableCodeParams &params) {
  vec_ZZ_pE out;
  mul(out, msg, params.G0);  // row-vector * matrix
  return out;
}

void EncodeRec(vec_ZZ_pE &out, const vec_ZZ_pE &msg, long level,
               const FoldableCodeParams &params) {
  if (level == 0) {
    out = Encode0(msg, params);
    return;
  }

  const long half_k = params.k0 * Pow2(level - 1);
  if (msg.length() != 2 * half_k) {
    LogicError("EncodeFoldable: message length mismatch at recursion level");
  }

  vec_ZZ_pE ml, mr;
  ml.SetLength(half_k);
  mr.SetLength(half_k);
  for (long i = 0; i < half_k; ++i) {
    ml[i] = msg[i];
    mr[i] = msg[i + half_k];
  }

  vec_ZZ_pE l, r;
  EncodeRec(l, ml, level - 1, params);
  EncodeRec(r, mr, level - 1, params);

  const vec_ZZ_pE &t =
      params.diag_T[static_cast<std::size_t>(level - 1)];  // length n_{level-1}
  const long half_n = l.length();
  if (t.length() != half_n || r.length() != half_n) {
    LogicError("EncodeFoldable: internal codeword length mismatch");
  }

  out.SetLength(2 * half_n);
  for (long i = 0; i < half_n; ++i) {
    const NTL::ZZ_pE tr = t[i] * r[i];
    out[i] = l[i] + tr;
    out[i + half_n] = l[i] + params.zeta * tr;
  }
}

void EncodeRecUnchecked(vec_ZZ_pE &out, const vec_ZZ_pE &msg, long level,
                        const FoldableCodeParams &params) {
  if (level == 0) {
    mul(out, msg, params.G0);  // row-vector * matrix
    return;
  }

  const long half_k = msg.length() / 2;

  vec_ZZ_pE ml, mr;
  ml.SetLength(half_k);
  mr.SetLength(half_k);
  for (long i = 0; i < half_k; ++i) {
    ml[i] = msg[i];
    mr[i] = msg[i + half_k];
  }

  vec_ZZ_pE l, r;
  EncodeRecUnchecked(l, ml, level - 1, params);
  EncodeRecUnchecked(r, mr, level - 1, params);

  const vec_ZZ_pE &t =
      params.diag_T[static_cast<std::size_t>(level - 1)];  // length n_{level-1}
  const long half_n = l.length();

  out.SetLength(2 * half_n);
  for (long i = 0; i < half_n; ++i) {
    const NTL::ZZ_pE tr = t[i] * r[i];
    out[i] = l[i] + tr;
    out[i + half_n] = l[i] + params.zeta * tr;
  }
}

}  // namespace

long MessageLength(const FoldableCodeParams &params) {
  ValidateParams(params);
  return params.k0 * Pow2(params.d);
}

long CodewordLength(const FoldableCodeParams &params) {
  ValidateParams(params);
  return params.c * MessageLength(params);
}

void EncodeFoldable(vec_ZZ_pE &out, const vec_ZZ_pE &msg,
                    const FoldableCodeParams &params) {
  ValidateParams(params);

  const long kd = params.k0 * Pow2(params.d);
  if (msg.length() != kd) {
    LogicError("EncodeFoldable: msg has wrong length");
  }

  EncodeRec(out, msg, params.d, params);
}

void EncodeFoldableUnchecked(vec_ZZ_pE &out, const vec_ZZ_pE &msg,
                             const FoldableCodeParams &params) {
  EncodeRecUnchecked(out, msg, params.d, params);
}

}  // namespace basefold
