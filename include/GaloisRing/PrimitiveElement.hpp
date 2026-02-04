#ifndef INCLUDE_GALOISRING_PRIMITIVEELEMENT_HPP_
#define INCLUDE_GALOISRING_PRIMITIVEELEMENT_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/vector.h>

// Returns a candidate "primitive element" for GR(p^k, s)
// (implementation-specific).
//
// Notes:
//   - The current implementation uses a hard-coded extension polynomial F in
//   PrimitiveElement.cpp.
//     You likely need to modify F to match your chosen (p, s).
//   - The function re-initializes ZZ_p/ZZ_pE contexts internally (changes
//   global NTL modulus contexts).
//
// Call: ZZ_pE alpha = FindPrimitiveElement(p, k, s);
NTL::ZZ_pE FindPrimitiveElement(NTL::ZZ p, long k, long s);

#endif  // INCLUDE_GALOISRING_PRIMITIVEELEMENT_HPP_
