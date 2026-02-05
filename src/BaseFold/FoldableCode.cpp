#include "BaseFold/FoldableCode.hpp"

#include <NTL/mat_ZZ_pE.h>

using NTL::LogicError;
using NTL::mul;
using NTL::vec_ZZ_pE;

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

void ValidateParams(const FoldableCodeParams &params) {
  if (params.c <= 0) LogicError("EncodeFoldable: c must be positive");
  if (params.k0 <= 0) LogicError("EncodeFoldable: k0 must be positive");
  if (params.d < 0) LogicError("EncodeFoldable: d must be non-negative");

  const long n0 = params.c * params.k0;
  if (params.G0.NumRows() != params.k0 || params.G0.NumCols() != n0) {
    LogicError("EncodeFoldable: G0 has wrong dimensions");
  }

  if (params.zeta == 0) LogicError("EncodeFoldable: zeta must be non-zero");
  if (params.zeta == 1) LogicError("EncodeFoldable: zeta must not equal 1");

  if (static_cast<long>(params.diag_T.size()) != params.d) {
    LogicError("EncodeFoldable: diag_T must have size d");
  }

  for (long i = 0; i < params.d; ++i) {
    const long ni = params.c * params.k0 * Pow2(i);
    if (params.diag_T[static_cast<std::size_t>(i)].length() != ni) {
      LogicError("EncodeFoldable: diag_T[i] has wrong length");
    }
    for (long j = 0; j < ni; ++j) {
      if (params.diag_T[static_cast<std::size_t>(i)][j] == 0) {
        LogicError("EncodeFoldable: diag_T[i] entries must be non-zero");
      }
    }
  }
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

}  // namespace basefold

