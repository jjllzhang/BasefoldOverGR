#include "GaloisRing/HenselLift.hpp"
#include "GaloisRing/PrimitiveElement.hpp"
#include "GaloisRing/utils.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace NTL;
using namespace std;

namespace {

int g_failures = 0;

/*
    Minimal assertion helper used by this test runner.

    Usage:
      - Prefer CHECK(...), CHECK_EQ(...), CHECK_MSG(...) macros below.
      - On failure, increments a global failure counter and prints a message to stderr.
*/
void Check(bool condition, const string& message, const char* file, int line) {
    if (condition) return;
    cerr << file << ":" << line << " FAIL: " << message << "\n";
    g_failures++;
}

#define CHECK(cond) Check((cond), #cond, __FILE__, __LINE__)
#define CHECK_MSG(cond, msg) Check((cond), (msg), __FILE__, __LINE__)
#define CHECK_EQ(a, b) Check(((a) == (b)), string(#a) + " == " + string(#b), __FILE__, __LINE__)

/*
    Constructs a constant ZZ_pE element (value as degree-0 polynomial) in the current ZZ_pE context.
    Preconditions: ZZ_p::init(modulus) and ZZ_pE::init(F) are already set.
    Usage: ZZ_pE one = ConstZZpE(1);
*/
ZZ_pE ConstZZpE(long value) {
    ZZ_pX poly;
    SetCoeff(poly, 0, to_ZZ_p(value));
    ZZ_pE out;
    conv(out, poly);
    return out;
}

/*
    Serializes a ZZ_pX into vector<long> coefficients.
    Usage: auto coeffs = ToLongVec(poly);
*/
vector<long> ToLongVec(const ZZ_pX& poly) {
    vector<long> out;
    ZZpX2long(poly, out);
    return out;
}

/*
    Serializes a ZZ_pE element into vector<long> of length s (basis coefficients).
    Usage: auto coeffs = ToLongVec(a, s);
*/
vector<long> ToLongVec(const ZZ_pE& element, long s) {
    vector<long> out;
    ZzpE2Veclong(element, out, s);
    return out;
}

/*
    Flattens a ZZ_pEX into vector<long>, expanding each ZZ_pE coefficient into s longs.
    Usage: auto flat = ToLongVec(f, s);
*/
vector<long> ToLongVec(const ZZ_pEX& poly, long s) {
    vector<long> out;
    ZZpEX2long(poly, out, s);
    return out;
}

/*
    Tests "plain" helper functions that do not depend on any NTL modulus/extension context.
*/
void TestUtilsBasics() {
    // FindFactor
    CHECK_EQ(FindFactor(8), (vector<long>{1, 2, 4, 8}));

    // isPowerOfTwo / nextPowerOf2
    CHECK(isPowerOfTwo(1));
    CHECK(isPowerOfTwo(2));
    CHECK(!isPowerOfTwo(3));
    CHECK(!isPowerOfTwo(0));
    CHECK(!isPowerOfTwo(-4));
    CHECK_EQ(nextPowerOf2(1), 1);
    CHECK_EQ(nextPowerOf2(5), 8);
    CHECK_EQ(nextPowerOf2(8), 8);
    CHECK_EQ(nextPowerOf2(0), 1);

    // Veclong2String
    CHECK_EQ(Veclong2String({1, 0, 23}), string("1023"));

    // Pad / Trim
    CHECK_EQ(PadVectorToLength({1, 2, 3}, 5), (vector<long>{1, 2, 3, 0, 0}));
    CHECK_EQ(PadVectorToLength({1, 2, 3}, 2), (vector<long>{1, 2}));
    CHECK_EQ(TrimVector({1, 2, 0, 0}), (vector<long>{1, 2}));
    CHECK_EQ(TrimVector({0, 0, 0}), (vector<long>{}));
    CHECK_EQ(TrimVector({0}), (vector<long>{}));
    CHECK_EQ(TrimVector({}), (vector<long>{}));

    // SplitAndPadVector: only test the "well-formed" case where input.size() is divisible by numSegments.
    CHECK_EQ(SplitAndPadVector({1, 2, 3, 4, 5, 6, 7, 8}, /*segmentLength=*/6, /*numSegments=*/2),
             (vector<long>{1, 2, 3, 4, 0, 0, 5, 6, 7, 8, 0, 0}));

    // splitVector
    CHECK_EQ(splitVector({1, 2, 3, 4, 5, 6, 7}, 3),
             (vector<vector<long>>{{1, 2, 3}, {4, 5, 6}, {7, 0, 0}}));

    // nearestPerfectSquare (as implemented: returns the smallest perfect square strictly greater than num)
    CHECK_EQ(nearestPerfectSquare(0), 1);
    CHECK_EQ(nearestPerfectSquare(1), 4);
    CHECK_EQ(nearestPerfectSquare(5), 9);

    // print (smoke)
    vector<long> tmp = {1, 2, 3};
    print(tmp);
}

/*
    Tests conversion helpers between vector<long> and NTL types in a field setting:
      - ZZ_p with modulus p
      - ZZ_pE as an extension of degree s via an irreducible polynomial
*/
void TestConversionsAndPolyOps_Field() {
    // Work in GF(p) with extension degree s over an irreducible polynomial.
    const ZZ p = to_ZZ(7);
    const long s = 2;

    ZZ_pPush p_push(p);

    // Use F(x) = x^2 + 1, irreducible over GF(7).
    ZZ_pX F;
    SetCoeff(F, 2, 1);
    SetCoeff(F, 0, 1);
    ZZ_pEPush e_push(F);

    // long <-> ZZ_pX
    {
        const vector<long> coeffs = {1, 2, 6};
        const ZZ_pX poly = long2ZZpX(coeffs);
        CHECK_EQ(ToLongVec(poly), coeffs);
    }

    // long <-> ZZ_pE
    {
        const vector<long> coeffs = {3, 5};  // 3 + 5*x
        const ZZ_pE a = long2ZZpE(coeffs);
        CHECK_EQ(ToLongVec(a, s), coeffs);
    }

    // vec_ZZ_pE -> vector<long> / vector<string>
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

    // ZZ_pEX <-> vector<long>
    {
        const vector<long> packed = {1, 2, 3, 4, 5, 6};  // 3 coefficients, each uses s=2 longs
        ZZ_pEX poly;
        clear(poly);
        Long2ZZpEX(packed, poly, s);
        CHECK_EQ(ToLongVec(poly, s), packed);

        ZZ_pEX poly2;
        clear(poly2);
        Long2ZZpEX2(packed, poly2, s, /*n=*/3);
        CHECK_EQ(ToLongVec(poly2, s), packed);
    }

    // fillIrred / fillInterpolation
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

/*
    Tests Inv() and Inv2() in a field setting (k == 1), verifying a * inv(a) == 1.
*/
void TestInverse_Field() {
    const ZZ p = to_ZZ(7);
    const long s = 2;

    ZZ_pPush p_push(p);

    // F(x) = x^2 + 1 over GF(7)
    ZZ_pX F;
    SetCoeff(F, 2, 1);
    SetCoeff(F, 0, 1);
    ZZ_pEPush e_push(F);

    ZZ_pE x = ConstZZpE(0);
    {
        ZZ_pX H;
        SetCoeff(H, 1, 1);
        conv(x, H);
    }

    const ZZ_pE a = x + ConstZZpE(3);  // x + 3
    CHECK_MSG(a != 0, "a must be non-zero");

    const ZZ_pE a_inv_custom = Inv(a, s);
    ZZ_pE one;
    set(one);
    CHECK_EQ(a * a_inv_custom, one);

    const ZZ_pE a_inv2 = Inv2(a, F, p, s, /*k=*/1);
    CHECK_EQ(a * a_inv2, one);
    CHECK_EQ(a_inv2, inv(a));
}

/*
    Tests FindPrimitivePoly() at a small parameter size by checking the output is monic and irreducible.
*/
void TestFindPrimitivePoly() {
    const ZZ p = to_ZZ(2);
    const long n = 3;

    ZZ_pPush p_push(p);

    ZZ_pX g;
    clear(g);
    FindPrimitivePoly(g, p, n);

    CHECK_EQ(deg(g), n);
    CHECK_EQ(LeadCoeff(g), ZZ_p(1));
    CHECK_EQ(DetIrredTest(g), 1);
}

/*
    Tests interpolate_for_GR():
      - Field case (l == 1): must match NTL::interpolate exactly.
      - Ring case (l > 1): checks that the interpolated polynomial evaluates to the given values.
*/
void TestInterpolateForGR() {
    // Field case: l == 1 should match NTL's interpolate exactly.
    {
        const ZZ p = to_ZZ(7);
        const long s = 2;
        ZZ_pPush p_push(p);

        ZZ_pX F;
        SetCoeff(F, 2, 1);
        SetCoeff(F, 0, 1);
        ZZ_pEPush e_push(F);

        vec_ZZ_pE a, b;
        a.SetLength(3);
        b.SetLength(3);
        a[0] = ConstZZpE(0);
        a[1] = ConstZZpE(1);
        a[2] = ConstZZpE(2);
        b[0] = ConstZZpE(1);
        b[1] = ConstZZpE(4);
        b[2] = ConstZZpE(2);

        ZZ_pEX f_ntl;
        interpolate(f_ntl, a, b);

        ZZ_pEX f_custom;
        interpolate_for_GR(f_custom, a, b, p, /*l=*/1, s);

        CHECK_EQ(f_custom, f_ntl);
    }

    // Ring case: l > 1. Use s == 1 and choose points with distinct residues mod p
    // so that all required inverses exist (units).
    {
        const ZZ p = to_ZZ(5);
        const long l = 2;       // modulus p^l = 25
        const long s = 1;       // extension degree 1 (elements are constants)
        const ZZ mod = power(p, l);

        ZZ_pPush p_push(mod);

        // Degree-1 modulus polynomial for ZZ_pE: F(x) = x.
        ZZ_pX F;
        SetCoeff(F, 1, 1);
        ZZ_pEPush e_push(F);

        vec_ZZ_pE a, b;
        a.SetLength(3);
        b.SetLength(3);

        a[0] = ConstZZpE(0);
        a[1] = ConstZZpE(1);
        a[2] = ConstZZpE(2);

        // b = a^2 + 1 (in Z_25)
        for (long i = 0; i < 3; i++) {
            b[i] = a[i] * a[i] + ConstZZpE(1);
        }

        ZZ_pEX f_custom;
        interpolate_for_GR(f_custom, a, b, p, l, s);

        for (long i = 0; i < 3; i++) {
            CHECK_EQ(eval(f_custom, a[i]), b[i]);
        }
    }
}

/*
    A small Hensel lifting smoke test on f(x) = x^2 - 1 over p = 3.
    Verifies that the lifted factor divides f modulo p^(n+1) and reduces to the original factor modulo p.
*/
void TestHenselLift_Smoke() {
    // A small Hensel lifting smoke test on f(x) = x^2 - 1 over p=3.
    const ZZ p = to_ZZ(3);
    const long n = 1;  // lift to p^(n+1) = 9
    const ZZ mod = power(p, n + 1);

    ZZ_pPush p_push(p);

    ZZ_pX f_mod_p;
    SetCoeff(f_mod_p, 2, 1);
    SetCoeff(f_mod_p, 0, -1);

    ZZ_pX g_mod_p;
    SetCoeff(g_mod_p, 1, 1);
    SetCoeff(g_mod_p, 0, -1);  // x - 1

    // Store expected g mod p coefficients as ZZ for later comparison without changing moduli.
    const ZZ expected_c0_mod_p = rep(coeff(g_mod_p, 0));  // 2
    const ZZ expected_c1_mod_p = rep(coeff(g_mod_p, 1));  // 1

    ZZ_pX g_lift;
    clear(g_lift);
    HenselLift(g_lift, f_mod_p, g_mod_p, p, n);

    CHECK_EQ(ZZ_p::modulus(), mod);

    // Verify divisibility in Z/(p^(n+1)) [x] by rebuilding f in the current modulus.
    ZZ_pX f_lift;
    SetCoeff(f_lift, 2, 1);
    SetCoeff(f_lift, 0, -1);

    ZZ_pX q, r;
    DivRem(q, r, f_lift, g_lift);
    CHECK(IsZero(r));

    // Verify g_lift ≡ g (mod p) coefficient-wise.
    auto ReduceCoeffModP = [&](const ZZ_pX& poly, long i) -> ZZ {
        ZZ c = rep(coeff(poly, i));
        c %= p;
        if (c < 0) c += p;
        return c;
    };

    CHECK_EQ(ReduceCoeffModP(g_lift, 0), expected_c0_mod_p % p);
    CHECK_EQ(ReduceCoeffModP(g_lift, 1), expected_c1_mod_p % p);
}

/*
    A primitive element smoke test: for k == 1, FindPrimitiveElement returns x (mod F) for the hard-coded F.
*/
void TestPrimitiveElement_Smoke() {
    // With k==1, FindPrimitiveElement returns b = x mod F, avoiding composite modulus contexts.
    const ZZ p = to_ZZ(2);
    const long k = 1;
    const long s = 30;

    ZZ_pPush p_push(p);

    ZZ_pE pe = FindPrimitiveElement(p, k, s);

    ZZ_pX H;
    SetCoeff(H, 1, 1);
    ZZ_pE expected;
    conv(expected, H);

    CHECK_EQ(pe, expected);
    CHECK(deg(rep(pe)) <= 1);
}

}  // namespace

/*
    Entry point for the test runner.
    Returns 0 on success, 1 if any CHECK failed, and 2 on an unexpected exception.
*/
int main() {
    try {
        TestUtilsBasics();
        TestConversionsAndPolyOps_Field();
        TestInverse_Field();
        TestFindPrimitivePoly();
        TestInterpolateForGR();
        TestHenselLift_Smoke();
        TestPrimitiveElement_Smoke();
    } catch (const exception& e) {
        cerr << "Unhandled std::exception: " << e.what() << "\n";
        return 2;
    } catch (...) {
        cerr << "Unhandled non-std exception\n";
        return 2;
    }

    if (g_failures == 0) {
        cout << "All tests passed.\n";
        return 0;
    }

    cerr << g_failures << " test(s) failed.\n";
    return 1;
}
