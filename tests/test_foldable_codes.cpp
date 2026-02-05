#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "BaseFold/FoldableCode.hpp"
#include "test_common.hpp"

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
using std::ostringstream;
using std::string;
using std::size_t;
using std::vector;

int g_failures = 0;

namespace {

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

} // namespace

int main() {
  try {
    RUN_TEST(TestFoldableEncode_OverBinaryField);
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
