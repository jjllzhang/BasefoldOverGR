#include <exception>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "BaseFold/FoldableCode.hpp"
#include "tests/test_common.hpp"

using NTL::conv;
using NTL::coeff;
using NTL::mat_ZZ_pE;
using NTL::mul;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::rep;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using std::cerr;
using std::cout;
using std::exception;
using std::mt19937;
using std::ostringstream;
using std::uniform_int_distribution;
using std::string;
using std::size_t;
using std::vector;

int g_failures = 0;

namespace {

ZZ_pE ElemFromIndexGF4(int idx, const ZZ_pE &alpha) {
  const ZZ_pE one = testutil::ConstZZpE(1);
  if (idx == 0) return ZZ_pE(0);
  if (idx == 1) return one;
  if (idx == 2) return alpha;
  return alpha + one;
}

ZZ_pE NonZeroElemFromIndexGF4(int idx, const ZZ_pE &alpha) {
  const ZZ_pE one = testutil::ConstZZpE(1);
  if (idx == 0) return one;
  if (idx == 1) return alpha;
  return alpha + one;
}

ZZ_pE ElemFromCoeffs(const vector<long> &coeffs) {
  ZZ_pX poly;
  for (size_t i = 0; i < coeffs.size(); ++i) {
    SetCoeff(poly, static_cast<long>(i), coeffs[i]);
  }
  ZZ_pE out;
  conv(out, poly);
  return out;
}

string ElemToCoeffString(const ZZ_pE &element) {
  const long s = ZZ_pE::degree();
  const ZZ_pX poly = rep(element);

  ostringstream out;
  out << "[";
  for (long i = 0; i < s; ++i) {
    long c = 0;
    conv(c, coeff(poly, i));
    out << c;
    if (i + 1 < s)
      out << ",";
  }
  out << "]";
  return out.str();
}

void PrintVec(const char *label, const vec_ZZ_pE &vec) {
  cout << "  " << label << " (len=" << vec.length() << "):\n";
  for (long i = 0; i < vec.length(); ++i) {
    cout << "    " << i << ": " << ElemToCoeffString(vec[i]) << "\n";
  }
}

mat_ZZ_pE BuildVandermondeG0(long k0, const vec_ZZ_pE &points) {
  mat_ZZ_pE G0;
  const long n0 = points.length();
  G0.SetDims(k0, n0);

  for (long j = 0; j < n0; ++j) {
    ZZ_pE cur = testutil::ConstZZpE(1);
    for (long r = 0; r < k0; ++r) {
      G0[r][j] = cur;
      cur *= points[j];
    }
  }

  return G0;
}

vector<vec_ZZ_pE> SampleDiagT(long c, long k0, long d, const ZZ_pE &alpha,
                              mt19937 &rng) {
  uniform_int_distribution<int> dist_nonzero(0, 2);
  vector<vec_ZZ_pE> diag_T;
  diag_T.resize(static_cast<size_t>(d));

  for (long level = 0; level < d; ++level) {
    const long ni = c * k0 * (1L << level);
    diag_T[static_cast<size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      const int choice = dist_nonzero(rng);
      diag_T[static_cast<size_t>(level)][i] =
          NonZeroElemFromIndexGF4(choice, alpha);
    }
  }

  return diag_T;
}

vec_ZZ_pE SampleMessage(long k0, long d, const ZZ_pE &alpha, mt19937 &rng) {
  uniform_int_distribution<int> dist_any(0, 3);
  vec_ZZ_pE msg;
  msg.SetLength(k0 * (1L << d));
  for (long i = 0; i < msg.length(); ++i) {
    msg[i] = ElemFromIndexGF4(dist_any(rng), alpha);
  }
  return msg;
}

ZZ_pE SampleElementGR42(mt19937 &rng) {
  uniform_int_distribution<int> dist_coeff(0, 3);
  const long c0 = dist_coeff(rng);
  const long c1 = dist_coeff(rng);
  return ElemFromCoeffs({c0, c1});
}

ZZ_pE SampleUnitGR42(mt19937 &rng) {
  uniform_int_distribution<int> dist_coeff(0, 3);
  while (true) {
    const long c0 = dist_coeff(rng);
    const long c1 = dist_coeff(rng);
    if ((c0 % 2 == 0) && (c1 % 2 == 0)) continue;
    return ElemFromCoeffs({c0, c1});
  }
}

vector<vec_ZZ_pE> SampleDiagT_GR42(long c, long k0, long d, mt19937 &rng) {
  vector<vec_ZZ_pE> diag_T;
  diag_T.resize(static_cast<size_t>(d));

  for (long level = 0; level < d; ++level) {
    const long ni = c * k0 * (1L << level);
    diag_T[static_cast<size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      diag_T[static_cast<size_t>(level)][i] = SampleUnitGR42(rng);
    }
  }

  return diag_T;
}

vec_ZZ_pE SampleMessageGR42(long k0, long d, mt19937 &rng) {
  vec_ZZ_pE msg;
  msg.SetLength(k0 * (1L << d));
  for (long i = 0; i < msg.length(); ++i) {
    msg[i] = SampleElementGR42(rng);
  }
  return msg;
}

mat_ZZ_pE BuildFoldableGeneratorMatrix(const mat_ZZ_pE &G0,
                                       const vector<vec_ZZ_pE> &diag_T,
                                       const ZZ_pE &zeta) {
  mat_ZZ_pE G = G0;

  for (size_t level = 0; level < diag_T.size(); ++level) {
    const vec_ZZ_pE &t = diag_T[level];
    const long k = G.NumRows();
    const long n = G.NumCols();
    CHECK_EQ(t.length(), n);

    mat_ZZ_pE next;
    next.SetDims(2 * k, 2 * n);

    for (long r = 0; r < k; ++r) {
      for (long c = 0; c < n; ++c) {
        const ZZ_pE &v = G[r][c];
        next[r][c] = v;
        next[r][c + n] = v;
        next[r + k][c] = v * t[c];
        next[r + k][c + n] = v * (zeta * t[c]);
      }
    }

    G = next;
  }

  return G;
}

void TestFoldableEncode_OverBinaryField() {
  testutil::PrintInfo(
      "Algorithm 1 (Encd): compares recursive EncodeFoldable vs msg * Gd");

  const ZZ p = to_ZZ(2);
  ZZ_pPush p_push(p);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  const long c = 2;
  const long k0 = 2;
  const long d = 2;
  const long n0 = c * k0;
  testutil::PrintInfo("Field: GF(2^2), F(x)=x^2+x+1");
  testutil::PrintInfo("Params: c=2, k0=2, d=2 => k_d=8, n_d=16");

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);
  const ZZ_pE zeta = alpha;

  vec_ZZ_pE points;
  points.SetLength(n0);
  points[0] = ZZ_pE(0);
  points[1] = testutil::ConstZZpE(1);
  points[2] = alpha;
  points[3] = alpha + testutil::ConstZZpE(1);

  mat_ZZ_pE G0;
  G0.SetDims(k0, n0);
  for (long j = 0; j < n0; ++j) {
    G0[0][j] = testutil::ConstZZpE(1);
    G0[1][j] = points[j];
  }

  vector<vec_ZZ_pE> diag_T;
  diag_T.resize(d);

  diag_T[0].SetLength(/*n0=*/4);
  diag_T[0][0] = testutil::ConstZZpE(1);
  diag_T[0][1] = alpha;
  diag_T[0][2] = alpha + testutil::ConstZZpE(1);
  diag_T[0][3] = alpha;

  diag_T[1].SetLength(/*n1=*/8);
  for (long i = 0; i < diag_T[1].length(); ++i) {
    if (i % 3 == 0)
      diag_T[1][i] = testutil::ConstZZpE(1);
    if (i % 3 == 1)
      diag_T[1][i] = alpha;
    if (i % 3 == 2)
      diag_T[1][i] = alpha + testutil::ConstZZpE(1);
  }

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = p;
  params.zeta = zeta;
  params.G0 = G0;
  params.diag_T = diag_T;

  vec_ZZ_pE msg;
  msg.SetLength(k0 * (1L << d));
  msg[0] = testutil::ConstZZpE(0);
  msg[1] = testutil::ConstZZpE(1);
  msg[2] = alpha;
  msg[3] = alpha + testutil::ConstZZpE(1);
  msg[4] = testutil::ConstZZpE(1);
  msg[5] = alpha;
  msg[6] = testutil::ConstZZpE(0);
  msg[7] = alpha + testutil::ConstZZpE(1);

  testutil::PrintInfo("Element format: [c0,c1] means c0 + c1*x in GF(2^2).");
  PrintVec("msg", msg);

  vec_ZZ_pE got;
  basefold::EncodeFoldable(got, msg, params);

  PrintVec("codeword", got);

  const mat_ZZ_pE Gd = BuildFoldableGeneratorMatrix(G0, diag_T, zeta);
  vec_ZZ_pE expected;
  mul(expected, msg, Gd);

  CHECK_EQ(got.length(), expected.length());
  for (long i = 0; i < got.length(); ++i) {
    CHECK_EQ(got[i], expected[i]);
  }
}

void TestFoldableEncode_OverGaloisRing() {
  testutil::PrintInfo(
      "Ring: GR(4,2), compares EncodeFoldable vs msg * Gd (with unit checks)");

  const ZZ p = to_ZZ(2);
  const ZZ mod = ZZ(4);
  ZZ_pPush mod_push(mod);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  const long c = 2;
  const long k0 = 2;
  const long d = 2;
  const long n0 = c * k0;
  testutil::PrintInfo("Ring: ZZ_p modulus = 4, extension degree = 2");
  testutil::PrintInfo("Params: c=2, k0=2, d=2 => k_d=8, n_d=16");

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);
  const ZZ_pE zeta = alpha;

  vec_ZZ_pE points;
  points.SetLength(n0);
  points[0] = ZZ_pE(0);
  points[1] = testutil::ConstZZpE(1);
  points[2] = alpha;
  points[3] = alpha + testutil::ConstZZpE(1);

  const mat_ZZ_pE G0 = BuildVandermondeG0(k0, points);

  vector<vec_ZZ_pE> diag_T;
  diag_T.resize(d);

  diag_T[0].SetLength(/*n0=*/4);
  diag_T[0][0] = testutil::ConstZZpE(1);
  diag_T[0][1] = alpha;
  diag_T[0][2] = alpha + testutil::ConstZZpE(1);
  diag_T[0][3] = alpha;

  diag_T[1].SetLength(/*n1=*/8);
  for (long i = 0; i < diag_T[1].length(); ++i) {
    if (i % 3 == 0)
      diag_T[1][i] = testutil::ConstZZpE(1);
    if (i % 3 == 1)
      diag_T[1][i] = alpha;
    if (i % 3 == 2)
      diag_T[1][i] = alpha + testutil::ConstZZpE(1);
  }

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = p;
  params.zeta = zeta;
  params.G0 = G0;
  params.diag_T = diag_T;

  const ZZ_pE two = testutil::ConstZZpE(2);
  vec_ZZ_pE msg;
  msg.SetLength(k0 * (1L << d));
  msg[0] = testutil::ConstZZpE(0);
  msg[1] = testutil::ConstZZpE(1);
  msg[2] = alpha;
  msg[3] = alpha + testutil::ConstZZpE(1);
  msg[4] = two;
  msg[5] = two * alpha;
  msg[6] = two * (alpha + testutil::ConstZZpE(1));
  msg[7] = testutil::ConstZZpE(3);

  testutil::PrintInfo("Element format: [c0,c1] means c0 + c1*x in GR(4,2).");
  PrintVec("msg", msg);

  vec_ZZ_pE got;
  basefold::EncodeFoldable(got, msg, params);
  PrintVec("codeword", got);

  const mat_ZZ_pE Gd = BuildFoldableGeneratorMatrix(G0, diag_T, zeta);
  vec_ZZ_pE expected;
  mul(expected, msg, Gd);

  CHECK_EQ(got.length(), expected.length());
  for (long i = 0; i < got.length(); ++i) {
    CHECK_EQ(got[i], expected[i]);
  }
}

void TestFoldableEncode_Randomized() {
  testutil::PrintInfo(
      "Randomized: multiple (c,k0,d), random diag_T and random messages");

  const ZZ p = to_ZZ(2);
  ZZ_pPush p_push(p);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const unsigned int seed = 1337U;
  testutil::PrintInfo("PRNG seed = 1337");
  mt19937 rng(seed);

  struct Case {
    long c;
    long k0;
    long d;
    int trials;
    int messages;
  };

  const vector<Case> cases = {
      {2, 2, 0, 4, 16}, {2, 2, 1, 4, 16}, {2, 2, 3, 3, 8},
      {1, 4, 2, 3, 8},  {4, 1, 4, 2, 8},  {1, 3, 2, 3, 8},
  };

  uniform_int_distribution<int> dist_zeta(0, 1);
  for (size_t case_id = 0; case_id < cases.size(); ++case_id) {
    const Case &cs = cases[case_id];
    const long n0 = cs.c * cs.k0;
    CHECK_LE(n0, 4);

    ostringstream hdr;
    hdr << "Case " << case_id << ": c=" << cs.c << ", k0=" << cs.k0
        << ", d=" << cs.d << " (n0=" << n0 << "), trials=" << cs.trials
        << ", messages=" << cs.messages;
    testutil::PrintInfo(hdr.str());

    vec_ZZ_pE points;
    points.SetLength(n0);
    for (long j = 0; j < n0; ++j) {
      points[j] = ElemFromIndexGF4(static_cast<int>(j), alpha);
    }

    const mat_ZZ_pE G0 = BuildVandermondeG0(cs.k0, points);

    for (int trial = 0; trial < cs.trials; ++trial) {
      const vector<vec_ZZ_pE> diag_T =
          SampleDiagT(cs.c, cs.k0, cs.d, alpha, rng);
      const ZZ_pE zeta =
          dist_zeta(rng) == 0 ? alpha : (alpha + testutil::ConstZZpE(1));

      basefold::FoldableCodeParams params;
      params.c = cs.c;
      params.k0 = cs.k0;
      params.d = cs.d;
      params.p = p;
      params.zeta = zeta;
      params.G0 = G0;
      params.diag_T = diag_T;

      const mat_ZZ_pE Gd = BuildFoldableGeneratorMatrix(G0, diag_T, zeta);

      for (int msg_id = 0; msg_id < cs.messages; ++msg_id) {
        const vec_ZZ_pE msg = SampleMessage(cs.k0, cs.d, alpha, rng);

        vec_ZZ_pE got;
        basefold::EncodeFoldable(got, msg, params);

        vec_ZZ_pE expected;
        mul(expected, msg, Gd);

        CHECK_EQ(got.length(), expected.length());
        for (long i = 0; i < got.length(); ++i) {
          if (got[i] != expected[i]) {
            ostringstream where;
            where << "Mismatch at idx=" << i << " (trial=" << trial
                  << ", msg_id=" << msg_id << ")";
            testutil::PrintInfo(where.str());
            PrintVec("msg", msg);
            PrintVec("got", got);
            PrintVec("expected", expected);
            CHECK_EQ(got[i], expected[i]);
            return;
          }
        }
      }
    }
  }
}

void TestFoldableEncode_Randomized_GR42() {
  testutil::PrintInfo(
      "Randomized ring: GR(4,2) with multiple (c,k0,d), random diag_T and random messages");

  const ZZ p = to_ZZ(2);
  const ZZ mod = ZZ(4);
  ZZ_pPush mod_push(mod);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const unsigned int seed = 424242U;
  testutil::PrintInfo("PRNG seed = 424242");
  mt19937 rng(seed);

  struct Case {
    long c;
    long k0;
    long d;
    int trials;
    int messages;
  };

  const vector<Case> cases = {
      {2, 2, 0, 4, 16}, {2, 2, 1, 4, 16}, {2, 2, 3, 3, 8},
      {1, 4, 2, 3, 8},  {4, 1, 4, 2, 8},  {1, 3, 2, 3, 8},
  };

  uniform_int_distribution<int> dist_zeta(0, 1);
  for (size_t case_id = 0; case_id < cases.size(); ++case_id) {
    const Case &cs = cases[case_id];
    const long n0 = cs.c * cs.k0;
    CHECK_LE(n0, 4);

    ostringstream hdr;
    hdr << "Case " << case_id << ": c=" << cs.c << ", k0=" << cs.k0
        << ", d=" << cs.d << " (n0=" << n0 << "), trials=" << cs.trials
        << ", messages=" << cs.messages;
    testutil::PrintInfo(hdr.str());

    vec_ZZ_pE points;
    points.SetLength(n0);
    for (long j = 0; j < n0; ++j) {
      points[j] = ElemFromIndexGF4(static_cast<int>(j), alpha);
    }

    const mat_ZZ_pE G0 = BuildVandermondeG0(cs.k0, points);

    for (int trial = 0; trial < cs.trials; ++trial) {
      const vector<vec_ZZ_pE> diag_T =
          SampleDiagT_GR42(cs.c, cs.k0, cs.d, rng);
      const ZZ_pE zeta =
          dist_zeta(rng) == 0 ? alpha : (alpha + testutil::ConstZZpE(1));

      basefold::FoldableCodeParams params;
      params.c = cs.c;
      params.k0 = cs.k0;
      params.d = cs.d;
      params.p = p;
      params.zeta = zeta;
      params.G0 = G0;
      params.diag_T = diag_T;

      const mat_ZZ_pE Gd = BuildFoldableGeneratorMatrix(G0, diag_T, zeta);

      for (int msg_id = 0; msg_id < cs.messages; ++msg_id) {
        const vec_ZZ_pE msg = SampleMessageGR42(cs.k0, cs.d, rng);

        vec_ZZ_pE got;
        basefold::EncodeFoldable(got, msg, params);

        vec_ZZ_pE expected;
        mul(expected, msg, Gd);

        CHECK_EQ(got.length(), expected.length());
        for (long i = 0; i < got.length(); ++i) {
          if (got[i] != expected[i]) {
            ostringstream where;
            where << "Mismatch at idx=" << i << " (trial=" << trial
                  << ", msg_id=" << msg_id << ")";
            testutil::PrintInfo(where.str());
            PrintVec("msg", msg);
            PrintVec("got", got);
            PrintVec("expected", expected);
            CHECK_EQ(got[i], expected[i]);
            return;
          }
        }
      }
    }
  }
}

} // namespace

int main() {
  try {
    RUN_TEST(TestFoldableEncode_OverBinaryField);
    RUN_TEST(TestFoldableEncode_OverGaloisRing);
    RUN_TEST(TestFoldableEncode_Randomized);
    RUN_TEST(TestFoldableEncode_Randomized_GR42);
  } catch (const exception &e) {
    cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    cerr << "Unhandled non-std exception\n";
    return 2;
  }

  if (g_failures == 0) {
    cout << "\nAll tests passed.\n";
    return 0;
  }

  cerr << "\n" << g_failures << " test(s) failed.\n";
  return 1;
}
