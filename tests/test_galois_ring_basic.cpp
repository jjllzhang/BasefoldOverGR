#include "GaloisRing/HenselLift.hpp"
#include "GaloisRing/PrimitiveElement.hpp"
#include "GaloisRing/utils.hpp"

#include "test_common.hpp"

#include <algorithm>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <vector>

using namespace NTL;
using namespace std;

int g_failures = 0;

namespace {

vector<long> ToLongVec(const ZZ_pX& poly) {
  vector<long> out;
  ZZpX2long(poly, out);
  return out;
}

vector<long> ToLongVec(const ZZ_pE& element, long s) {
  vector<long> out;
  ZzpE2Veclong(element, out, s);
  return out;
}

vector<long> ToLongVec(const ZZ_pEX& poly, long s) {
  vector<long> out;
  ZZpEX2long(poly, out, s);
  return out;
}

void TestUtilsBasics() {
  testutil::PrintInfo(
      "Utilities without NTL contexts (factors, padding, power-of-two helpers)");

  CHECK_EQ(FindFactor(8), (vector<long>{1, 2, 4, 8}));

  CHECK(isPowerOfTwo(1));
  CHECK(isPowerOfTwo(2));
  CHECK(!isPowerOfTwo(3));
  CHECK(!isPowerOfTwo(0));
  CHECK(!isPowerOfTwo(-4));
  CHECK_EQ(nextPowerOf2(1), 1);
  CHECK_EQ(nextPowerOf2(5), 8);
  CHECK_EQ(nextPowerOf2(8), 8);
  CHECK_EQ(nextPowerOf2(0), 1);

  CHECK_EQ(Veclong2String({1, 0, 23}), string("1023"));

  CHECK_EQ(PadVectorToLength({1, 2, 3}, 5), (vector<long>{1, 2, 3, 0, 0}));
  CHECK_EQ(PadVectorToLength({1, 2, 3}, 2), (vector<long>{1, 2}));
  CHECK_EQ(TrimVector({1, 2, 0, 0}), (vector<long>{1, 2}));
  CHECK_EQ(TrimVector({0, 0, 0}), (vector<long>{}));
  CHECK_EQ(TrimVector({0}), (vector<long>{}));
  CHECK_EQ(TrimVector({}), (vector<long>{}));

  CHECK_EQ(SplitAndPadVector({1, 2, 3, 4, 5, 6, 7, 8}, /*segmentLength=*/6,
                             /*numSegments=*/2),
           (vector<long>{1, 2, 3, 4, 0, 0, 5, 6, 7, 8, 0, 0}));

  CHECK_EQ(splitVector({1, 2, 3, 4, 5, 6, 7}, 3),
           (vector<vector<long>>{{1, 2, 3}, {4, 5, 6}, {7, 0, 0}}));

  CHECK_EQ(nearestPerfectSquare(0), 1);
  CHECK_EQ(nearestPerfectSquare(1), 4);
  CHECK_EQ(nearestPerfectSquare(5), 9);

  vector<long> tmp = {1, 2, 3};
  print(tmp);
}

void TestConversionsAndPolyOps_Field() {
  const ZZ p = to_ZZ(7);
  const long s = 2;
  testutil::PrintInfo("Field conversions in GF(7^2), s=2, F(x)=x^2+1");

  ZZ_pPush p_push(p);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  {
    const vector<long> coeffs = {1, 2, 6};
    const ZZ_pX poly = long2ZZpX(coeffs);
    CHECK_EQ(ToLongVec(poly), coeffs);
  }

  {
    const vector<long> coeffs = {3, 5};  // 3 + 5*x
    const ZZ_pE a = long2ZZpE(coeffs);
    CHECK_EQ(ToLongVec(a, s), coeffs);
  }

  {
    vec_ZZ_pE v;
    v.SetLength(2);
    v[0] = long2ZZpE({1, 2});
    v[1] = long2ZZpE({3, 4});

    vector<long> packed;
    VeczzpE2Veclong(v, packed, s);
    CHECK_EQ(packed, (vector<long>{1, 2, 3, 4}));

    vector<string> packed_str;
    VeczzpE2Vecstring(v, packed_str, s);
    CHECK_EQ(packed_str, (vector<string>{"12", "34"}));

    CHECK(allNonZero(v));
    v[1] = ZZ_pE(0);
    CHECK(!allNonZero(v));
  }

  {
    const vector<long> packed = {1, 2, 3, 4, 5, 6};
    ZZ_pEX poly;
    clear(poly);
    Long2ZZpEX(packed, poly, s);
    CHECK_EQ(ToLongVec(poly, s), packed);

    ZZ_pEX poly2;
    clear(poly2);
    Long2ZZpEX2(packed, poly2, s, /*n=*/3);
    CHECK_EQ(ToLongVec(poly2, s), packed);
  }

  {
    ZZ_pX poly;
    SetCoeff(poly, 0, 1);
    SetCoeff(poly, 1, 2);
    SetCoeff(poly, 3, 4);

    vector<long> v1;
    ZZpX2long(poly, v1);

    vector<long> v2;
    fillIrred(poly, v2);
    CHECK_EQ(v1, v2);

    vec_ZZ_pE pts;
    pts.SetLength(2);
    pts[0] = long2ZZpE({1, 0});
    pts[1] = long2ZZpE({2, 0});

    vector<long> w1;
    VeczzpE2Veclong(pts, w1, s);

    vector<long> w2;
    fillInterpolation(pts, w2, s);

    CHECK_EQ(w1, w2);
  }
}

void TestInverse_Field() {
  const ZZ p = to_ZZ(7);
  const long s = 2;
  testutil::PrintInfo("Inverse in GF(7^2): checks Inv() and Inv2() vs NTL inv");

  ZZ_pPush p_push(p);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  ZZ_pE x = testutil::ConstZZpE(0);
  {
    ZZ_pX H;
    SetCoeff(H, 1, 1);
    conv(x, H);
  }

  const ZZ_pE a = x + testutil::ConstZZpE(3);
  CHECK_MSG(a != 0, "a must be non-zero");

  const ZZ_pE a_inv_custom = Inv(a, s);
  ZZ_pE one;
  set(one);
  CHECK_EQ(a * a_inv_custom, one);

  const ZZ modulus_before = ZZ_p::modulus();
  const ZZ_pE a_inv2 = Inv2(a, F, p, s, /*k=*/1);
  CHECK_EQ(ZZ_p::modulus(), modulus_before);
  CHECK_EQ(a * a_inv2, one);
  CHECK_EQ(a_inv2, inv(a));
}

void TestFindPrimitivePoly() {
  const ZZ p = to_ZZ(2);
  const long n = 3;
  testutil::PrintInfo("FindPrimitivePoly in GF(2): checks monic + irreducible");

  ZZ_pPush p_push(p);

  ZZ_pX g;
  clear(g);
  FindPrimitivePoly(g, p, n);

  CHECK_EQ(deg(g), n);
  CHECK_EQ(LeadCoeff(g), ZZ_p(1));
  CHECK_EQ(DetIrredTest(g), 1);
}

void TestInterpolateForGR() {
  testutil::PrintInfo("interpolate_for_GR: field case equals NTL; ring case evaluates correctly");

  {
    const ZZ p = to_ZZ(7);
    const long s = 2;
    testutil::PrintInfo("  Case A: field (l=1), p=7, s=2");

    ZZ_pPush p_push(p);

    ZZ_pX F;
    SetCoeff(F, 2, 1);
    SetCoeff(F, 0, 1);
    ZZ_pEPush e_push(F);

    vec_ZZ_pE a, b;
    a.SetLength(3);
    b.SetLength(3);
    a[0] = testutil::ConstZZpE(0);
    a[1] = testutil::ConstZZpE(1);
    a[2] = testutil::ConstZZpE(2);
    b[0] = testutil::ConstZZpE(1);
    b[1] = testutil::ConstZZpE(4);
    b[2] = testutil::ConstZZpE(2);

    ZZ_pEX f_ntl;
    interpolate(f_ntl, a, b);

    ZZ_pEX f_custom;
    interpolate_for_GR(f_custom, a, b, p, /*l=*/1, s);

    CHECK_EQ(f_custom, f_ntl);
  }

  {
    const ZZ p = to_ZZ(5);
    const long l = 2;
    const long s = 1;
    const ZZ mod = power(p, l);
    testutil::PrintInfo("  Case B: ring (l=2), modulus=25, s=1");

    ZZ_pPush p_push(mod);

    ZZ_pX F;
    SetCoeff(F, 1, 1);
    ZZ_pEPush e_push(F);

    vec_ZZ_pE a, b;
    a.SetLength(3);
    b.SetLength(3);

    a[0] = testutil::ConstZZpE(0);
    a[1] = testutil::ConstZZpE(1);
    a[2] = testutil::ConstZZpE(2);

    for (long i = 0; i < 3; i++) {
      b[i] = a[i] * a[i] + testutil::ConstZZpE(1);
    }

    ZZ_pEX f_custom;
    interpolate_for_GR(f_custom, a, b, p, l, s);

    for (long i = 0; i < 3; i++) {
      CHECK_EQ(eval(f_custom, a[i]), b[i]);
    }
  }
}

void TestHenselLift_Smoke() {
  const ZZ p = to_ZZ(3);
  const long n = 1;
  const ZZ mod = power(p, n + 1);
  testutil::PrintInfo("HenselLift smoke: f=x^2-1, g=x-1, p=3, lift to 9");

  ZZ_pPush p_push(p);
  const ZZ modulus_before = ZZ_p::modulus();

  ZZ_pX f_mod_p;
  SetCoeff(f_mod_p, 2, 1);
  SetCoeff(f_mod_p, 0, -1);

  ZZ_pX g_mod_p;
  SetCoeff(g_mod_p, 1, 1);
  SetCoeff(g_mod_p, 0, -1);

  const ZZ expected_c0_mod_p = rep(coeff(g_mod_p, 0));
  const ZZ expected_c1_mod_p = rep(coeff(g_mod_p, 1));

  ZZ_pX g_lift;
  clear(g_lift);
  HenselLift(g_lift, f_mod_p, g_mod_p, p, n);

  CHECK_EQ(ZZ_p::modulus(), modulus_before);

  {
    ZZ_pPush mod_push(mod);
    ZZ_pX f_lift;
    SetCoeff(f_lift, 2, 1);
    SetCoeff(f_lift, 0, -1);

    ZZ_pX q, r;
    DivRem(q, r, f_lift, g_lift);
    CHECK(IsZero(r));
  }

  auto ReduceCoeffModP = [&](const ZZ_pX& poly, long i) -> ZZ {
    ZZ c = rep(coeff(poly, i));
    c %= p;
    if (c < 0) c += p;
    return c;
  };

  CHECK_EQ(ReduceCoeffModP(g_lift, 0), expected_c0_mod_p % p);
  CHECK_EQ(ReduceCoeffModP(g_lift, 1), expected_c1_mod_p % p);
}

void TestPrimitiveElement_Smoke() {
  const ZZ p = to_ZZ(2);
  const long k = 1;
  const long s = 30;
  const ZZ mod = power(p, k);
  testutil::PrintInfo("Primitive element smoke: p=2, k=1, s=30 (hard-coded F in implementation)");

  ZZ_pPush p_push(mod);

  ZZ_pX F;
  SetCoeff(F, 30, 1);
  SetCoeff(F, 16, 1);
  SetCoeff(F, 15, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  const ZZ modulus_before = ZZ_p::modulus();
  const ZZ_pX modulus_poly_before = ZZ_pE::modulus().val();
  const long degree_before = ZZ_pE::degree();

  ZZ_pE pe = FindPrimitiveElement(p, k, s);

  CHECK_EQ(ZZ_p::modulus(), modulus_before);
  CHECK(ZZ_pE::initialized());
  CHECK_EQ(ZZ_pE::degree(), degree_before);
  CHECK_EQ(ZZ_pE::modulus().val(), modulus_poly_before);

  ZZ_pX H;
  SetCoeff(H, 1, 1);
  ZZ_pE expected;
  conv(expected, H);

  CHECK_EQ(pe, expected);
  CHECK(deg(rep(pe)) <= 1);
}

}  // namespace

int main() {
  try {
    RUN_TEST(TestUtilsBasics);
    RUN_TEST(TestConversionsAndPolyOps_Field);
    RUN_TEST(TestInverse_Field);
    RUN_TEST(TestFindPrimitivePoly);
    RUN_TEST(TestInterpolateForGR);
    RUN_TEST(TestHenselLift_Smoke);
    RUN_TEST(TestPrimitiveElement_Smoke);
  } catch (const exception& e) {
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
