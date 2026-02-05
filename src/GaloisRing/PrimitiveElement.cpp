#include "GaloisRing/PrimitiveElement.hpp"

using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pBak;
using NTL::ZZ_pE;
using NTL::ZZ_pEBak;
using NTL::ZZ_pX;
using NTL::conv;
using NTL::power;

/*
    Returns a candidate "primitive element" over the Galois ring GR(p^k, s)
   (implementation-specific).

    IMPORTANT:
      - The extension/modulus polynomial F is hard-coded below. You will usually
   need to replace its coefficients with a polynomial that matches your chosen
   (p, s).
      - This function re-initializes NTL contexts (ZZ_p::init and ZZ_pE::init)
   internally, but restores the incoming contexts on return.

    Usage: ZZ_pE alpha = FindPrimitiveElement(p, k, s);
*/
ZZ_pE FindPrimitiveElement(ZZ p, long k, long s) {
  ZZ_pBak modulus_bak;
  modulus_bak.save();
  ZZ_pEBak extension_bak;
  extension_bak.save();

  ZZ q1 = power(p, s);
  ZZ q3 = power(p, k);
  ZZ q6 = power(p, k - 1);

  long q4;
  conv(q4, q1 - ZZ(1));
  long q5;
  conv(q5, q6);

  ZZ_p::init(p);
  ZZ_pX F;

  /////  user need to modify here
  SetCoeff(F, 30, 1);
  SetCoeff(F, 16, 1);
  SetCoeff(F, 15, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  /////

  ZZ_p::init(q3);
  ZZ_pE::init(F);

  ZZ_pX H;
  SetCoeff(H, 1, 1);
  ZZ_pE b;
  conv(b, H);

  return power(b, q5);
}
