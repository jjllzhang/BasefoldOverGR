#include "PCS/BaseFold/FoldableCode.hpp"
#include "PCS/Common/Profile.hpp"

#include <NTL/ZZ_p.h>
#include <NTL/mat_ZZ_pE.h>

#include <limits>

#if defined(BASEFOLD_USE_OPENMP)
#include <omp.h>
#endif

using NTL::ZZ;
using NTL::LogicError;
using NTL::clear;
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

template <typename Fn>
void ForEachIndexMaybeParallel(long begin, long end, long parallel_threshold,
                               const Fn &fn) {
  if (end <= begin) return;
#if defined(BASEFOLD_USE_OPENMP)
  const long work_items = end - begin;
  if (work_items >= parallel_threshold) {
    const int max_threads = omp_get_max_threads();
    int threads_to_use = static_cast<int>(work_items / parallel_threshold);
    if (threads_to_use > max_threads) threads_to_use = max_threads;
    if (threads_to_use >= 2) {
      const ZZ base_modulus = NTL::ZZ_p::modulus();
      const ZZ_pX extension_modulus = NTL::ZZ_pE::modulus().val();
#pragma omp parallel num_threads(threads_to_use)
      {
        NTL::ZZ_p::init(base_modulus);
        NTL::ZZ_pE::init(extension_modulus);
#pragma omp for schedule(static)
        for (long i = begin; i < end; ++i) {
          fn(i);
        }
      }
      return;
    }
  }
#endif
  for (long i = begin; i < end; ++i) {
    fn(i);
  }
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

vec_ZZ_pE EncodeLevel0WithG0(const vec_ZZ_pE &msg,
                             const FoldableCodeParams &params) {
  vec_ZZ_pE out;
  mul(out, msg, params.G0);  // row-vector * matrix
  return out;
}

bool IsAllOnes(const vec_ZZ_pE &v) {
  ZZ_pE one;
  set(one);
  for (long i = 0; i < v.length(); ++i) {
    if (v[i] != one) return false;
  }
  return true;
}

bool IsSystematicRepeatedIdentity(const FoldableCodeParams &params) {
  const long c = params.c;
  const long k0 = params.k0;
  if (c <= 0 || k0 <= 0) return false;
  if (params.G0.NumRows() != k0) return false;
  if (params.G0.NumCols() != c * k0) return false;

  ZZ_pE zero;
  clear(zero);
  ZZ_pE one;
  set(one);

  for (long block = 0; block < c; ++block) {
    const long base = block * k0;
    for (long col = 0; col < k0; ++col) {
      for (long row = 0; row < k0; ++row) {
        const ZZ_pE &entry = params.G0[row][base + col];
        if (row == col) {
          if (entry != one) return false;
        } else if (entry != zero) {
          return false;
        }
      }
    }
  }
  return true;
}

struct EncoderFastPathCache {
  bool valid = false;
  const FoldableCodeParams *params = nullptr;
  ZZ modulus;
  long degree = 0;
  long c = 0;
  long k0 = 0;
  long d = 0;
  long g0_rows = 0;
  long g0_cols = 0;
  std::vector<unsigned char> diag_t_all_ones;
  bool systematic_repeated_identity_g0 = false;
};

bool MatchesFastPathCache(const EncoderFastPathCache &cache,
                          const FoldableCodeParams &params) {
  return cache.valid && cache.params == &params &&
         cache.modulus == NTL::ZZ_p::modulus() &&
         cache.degree == ZZ_pE::degree() && cache.c == params.c &&
         cache.k0 == params.k0 && cache.d == params.d &&
         cache.g0_rows == params.G0.NumRows() &&
         cache.g0_cols == params.G0.NumCols();
}

const EncoderFastPathCache &GetEncoderFastPathCache(
    const FoldableCodeParams &params) {
  static thread_local EncoderFastPathCache cache;
  if (MatchesFastPathCache(cache, params)) {
    return cache;
  }

  cache.valid = true;
  cache.params = &params;
  cache.modulus = NTL::ZZ_p::modulus();
  cache.degree = ZZ_pE::degree();
  cache.c = params.c;
  cache.k0 = params.k0;
  cache.d = params.d;
  cache.g0_rows = params.G0.NumRows();
  cache.g0_cols = params.G0.NumCols();
  cache.diag_t_all_ones.assign(static_cast<std::size_t>(params.d), 0);
  for (long level = 0; level < params.d; ++level) {
    cache.diag_t_all_ones[static_cast<std::size_t>(level)] =
        IsAllOnes(params.diag_T[static_cast<std::size_t>(level)]) ? 1 : 0;
  }
  cache.systematic_repeated_identity_g0 = IsSystematicRepeatedIdentity(params);
  return cache;
}

void EncodeFoldable_k0_1_Iterative(vec_ZZ_pE &out, const vec_ZZ_pE &msg,
                                  const FoldableCodeParams &params) {
  const long c = params.c;
  const long d = params.d;
  const long kd = msg.length();

  if (c <= 0) LogicError("EncodeFoldable_k0_1_Iterative: c must be positive");
  if (d < 0) LogicError("EncodeFoldable_k0_1_Iterative: d must be non-negative");
  if (kd < 0) LogicError("EncodeFoldable_k0_1_Iterative: msg length invalid");

  if (kd > std::numeric_limits<long>::max() / c) {
    LogicError("EncodeFoldable_k0_1_Iterative: overflow in n_d");
  }
  const long nd = kd * c;

  out.SetLength(nd);

  const ZZ_pE &zeta = params.zeta;
  const long parallel_threshold = 1024;
  const EncoderFastPathCache &fast_path = GetEncoderFastPathCache(params);

  // Level 0: encode each scalar message symbol using the single-row G0.
  ForEachIndexMaybeParallel(0, kd, parallel_threshold, [&](long block) {
    const ZZ_pE &m = msg[block];
    const long base = block * c;
    for (long j = 0; j < c; ++j) {
      out[base + j] = m * params.G0[0][j];
    }
  });

  long block_len = c;  // n_0
  long blocks = kd;    // number of blocks at current level
  for (long level = 0; level < d; ++level) {
    const vec_ZZ_pE &t = params.diag_T[static_cast<std::size_t>(level)];
    const bool t_all_ones =
        fast_path.diag_t_all_ones[static_cast<std::size_t>(level)] != 0;
    const long new_block_len = 2 * block_len;

    const long half_blocks = blocks / 2;
    if (t_all_ones) {
      ForEachIndexMaybeParallel(0, half_blocks, parallel_threshold,
                                [&](long b) {
        const long base = b * new_block_len;
        const long left = base;
        const long right = base + block_len;
        for (long i = 0; i < block_len; ++i) {
          const ZZ_pE code_left = out[left + i];
          const ZZ_pE code_right = out[right + i];
          out[left + i] = code_left + code_right;
          out[right + i] = code_left + zeta * code_right;
        }
      });
    } else {
      ForEachIndexMaybeParallel(0, half_blocks, parallel_threshold,
                                [&](long b) {
        const long base = b * new_block_len;
        const long left = base;
        const long right = base + block_len;
        for (long i = 0; i < block_len; ++i) {
          const ZZ_pE code_left = out[left + i];
          const ZZ_pE twisted_right = t[i] * out[right + i];
          out[left + i] = code_left + twisted_right;
          out[right + i] = code_left + zeta * twisted_right;
        }
      });
    }
    block_len = new_block_len;
    blocks /= 2;
  }
}

void EncodeFoldable_k0_gt1_Iterative(vec_ZZ_pE &out, const vec_ZZ_pE &msg,
                                    const FoldableCodeParams &params) {
  const long c = params.c;
  const long k0 = params.k0;
  const long d = params.d;
  const long kd = msg.length();

  if (c <= 0) LogicError("EncodeFoldable_k0_gt1_Iterative: c must be positive");
  if (k0 <= 1) LogicError("EncodeFoldable_k0_gt1_Iterative: k0 must be > 1");
  if (d < 0)
    LogicError("EncodeFoldable_k0_gt1_Iterative: d must be non-negative");
  if (kd < 0)
    LogicError("EncodeFoldable_k0_gt1_Iterative: msg length invalid");

  if (kd > std::numeric_limits<long>::max() / c) {
    LogicError("EncodeFoldable_k0_gt1_Iterative: overflow in n_d");
  }
  const long nd = kd * c;

  const long n0 = params.G0.NumCols();
  if (n0 <= 0) LogicError("EncodeFoldable_k0_gt1_Iterative: invalid n0");
  if (kd % k0 != 0)
    LogicError("EncodeFoldable_k0_gt1_Iterative: msg length not multiple of k0");
  const long blocks0 = kd / k0;

  out.SetLength(nd);
  const long parallel_threshold = 1024;
  const EncoderFastPathCache &fast_path = GetEncoderFastPathCache(params);

  // Level 0: encode each length-k0 message block using the k0 x n0 generator.
  // out is laid out as consecutive codewords of length n0.
  if (fast_path.systematic_repeated_identity_g0) {
    ForEachIndexMaybeParallel(0, blocks0, parallel_threshold, [&](long block) {
      const long msg_base = block * k0;
      const long out_base = block * n0;
      for (long copy = 0; copy < c; ++copy) {
        const long copy_base = out_base + copy * k0;
        for (long r = 0; r < k0; ++r) {
          out[copy_base + r] = msg[msg_base + r];
        }
      }
    });
  } else {
    ForEachIndexMaybeParallel(0, blocks0, parallel_threshold, [&](long block) {
      const long msg_base = block * k0;
      const long out_base = block * n0;

      for (long j = 0; j < n0; ++j) {
        clear(out[out_base + j]);
      }
      for (long r = 0; r < k0; ++r) {
        const ZZ_pE &m = msg[msg_base + r];
        for (long j = 0; j < n0; ++j) {
          out[out_base + j] += m * params.G0[r][j];
        }
      }
    });
  }

  const ZZ_pE &zeta = params.zeta;

  long block_len = n0;     // n_0
  long blocks = blocks0;   // number of blocks at current level
  for (long level = 0; level < d; ++level) {
    const vec_ZZ_pE &t = params.diag_T[static_cast<std::size_t>(level)];
    const bool t_all_ones =
        fast_path.diag_t_all_ones[static_cast<std::size_t>(level)] != 0;
    const long new_block_len = 2 * block_len;

    const long half_blocks = blocks / 2;
    if (t_all_ones) {
      ForEachIndexMaybeParallel(0, half_blocks, parallel_threshold,
                                [&](long b) {
        const long base = b * new_block_len;
        const long left = base;
        const long right = base + block_len;
        for (long i = 0; i < block_len; ++i) {
          const ZZ_pE code_left = out[left + i];
          const ZZ_pE code_right = out[right + i];
          out[left + i] = code_left + code_right;
          out[right + i] = code_left + zeta * code_right;
        }
      });
    } else {
      ForEachIndexMaybeParallel(0, half_blocks, parallel_threshold,
                                [&](long b) {
        const long base = b * new_block_len;
        const long left = base;
        const long right = base + block_len;
        for (long i = 0; i < block_len; ++i) {
          const ZZ_pE code_left = out[left + i];
          const ZZ_pE twisted_right = t[i] * out[right + i];
          out[left + i] = code_left + twisted_right;
          out[right + i] = code_left + zeta * twisted_right;
        }
      });
    }

    block_len = new_block_len;
    blocks /= 2;
  }
}

void EncodeRec(vec_ZZ_pE &out, const vec_ZZ_pE &msg, long level,
               const FoldableCodeParams &params) {
  if (level == 0) {
    out = EncodeLevel0WithG0(msg, params);
    return;
  }

  const long half_k = params.k0 * Pow2(level - 1);
  if (msg.length() != 2 * half_k) {
    LogicError("EncodeFoldable: message length mismatch at recursion level");
  }

  vec_ZZ_pE msg_left, msg_right;
  msg_left.SetLength(half_k);
  msg_right.SetLength(half_k);
  for (long i = 0; i < half_k; ++i) {
    msg_left[i] = msg[i];
    msg_right[i] = msg[i + half_k];
  }

  vec_ZZ_pE code_left, code_right;
  EncodeRec(code_left, msg_left, level - 1, params);
  EncodeRec(code_right, msg_right, level - 1, params);

  const vec_ZZ_pE &t =
      params.diag_T[static_cast<std::size_t>(level - 1)];  // length n_{level-1}
  const long half_n = code_left.length();
  if (t.length() != half_n || code_right.length() != half_n) {
    LogicError("EncodeFoldable: internal codeword length mismatch");
  }

  out.SetLength(2 * half_n);
  for (long i = 0; i < half_n; ++i) {
    const NTL::ZZ_pE twisted_right = t[i] * code_right[i];
    out[i] = code_left[i] + twisted_right;
    out[i + half_n] = code_left[i] + params.zeta * twisted_right;
  }
}

void EncodeRecUnchecked(vec_ZZ_pE &out, const vec_ZZ_pE &msg, long level,
                        const FoldableCodeParams &params) {
  if (level == 0) {
    mul(out, msg, params.G0);  // row-vector * matrix
    return;
  }

  const long half_k = msg.length() / 2;

  vec_ZZ_pE msg_left, msg_right;
  msg_left.SetLength(half_k);
  msg_right.SetLength(half_k);
  for (long i = 0; i < half_k; ++i) {
    msg_left[i] = msg[i];
    msg_right[i] = msg[i + half_k];
  }

  vec_ZZ_pE code_left, code_right;
  EncodeRecUnchecked(code_left, msg_left, level - 1, params);
  EncodeRecUnchecked(code_right, msg_right, level - 1, params);

  const vec_ZZ_pE &t =
      params.diag_T[static_cast<std::size_t>(level - 1)];  // length n_{level-1}
  const long half_n = code_left.length();

  out.SetLength(2 * half_n);
  for (long i = 0; i < half_n; ++i) {
    const NTL::ZZ_pE twisted_right = t[i] * code_right[i];
    out[i] = code_left[i] + twisted_right;
    out[i + half_n] = code_left[i] + params.zeta * twisted_right;
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

  if (params.k0 == 1) {
    EncodeFoldable_k0_1_Iterative(out, msg, params);
    return;
  }
  EncodeFoldable_k0_gt1_Iterative(out, msg, params);
}

void EncodeFoldableUnchecked(vec_ZZ_pE &out, const vec_ZZ_pE &msg,
                             const FoldableCodeParams &params) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->encode_foldable_unchecked_ns : nullptr,
                    prof ? &prof->encode_foldable_unchecked_calls : nullptr);

  if (params.k0 == 1) {
    EncodeFoldable_k0_1_Iterative(out, msg, params);
    return;
  }
  EncodeFoldable_k0_gt1_Iterative(out, msg, params);
}

}  // namespace basefold
