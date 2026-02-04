#ifndef GALOISRING_HENSELLIFT_HPP_
#define GALOISRING_HENSELLIFT_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/vector.h>

// Hensel-lifts a factor g of f from modulo p to modulo p^(n+1).
//
// Given polynomials f, g in Z_p[x] with g | f (mod p), this computes g_ in
// Z_{p^(n+1)}[x] such that g_ | f (mod p^(n+1)) and g_ ≡ g (mod p).
//
// Notes:
//   - This function calls ZZ_p::init(...) internally for p, p^2, ..., p^(n+1)
//   and thus changes
//     the global ZZ_p modulus context. Wrap the call in ZZ_pPush/ZZ_pBak if
//     needed.
//
// Call: HenselLift(g_lift, f_mod_p, g_mod_p, p, n);
void HenselLift(NTL::ZZ_pX &g_, const NTL::ZZ_pX &f, const NTL::ZZ_pX &g,
                const NTL::ZZ p, long n);

#endif  // GALOISRING_HENSELLIFT_HPP_
