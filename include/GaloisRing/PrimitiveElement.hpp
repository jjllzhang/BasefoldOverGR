#ifndef GALOISRING_PRIMITIVEELEMENT_HPP_
#define GALOISRING_PRIMITIVEELEMENT_HPP_

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
//   - The function re-initializes ZZ_p/ZZ_pE contexts internally, but preserves
//     the incoming contexts (it restores them on return). To use the returned
//     element, make sure the caller has set the matching GR(p^k, s) context.
//
// Call: ZZ_pE alpha = FindPrimitiveElement(p, k, s);
NTL::ZZ_pE FindPrimitiveElement(NTL::ZZ p, long k, long s);

#endif  // GALOISRING_PRIMITIVEELEMENT_HPP_
