#ifndef GALOISRING_BASIS_HPP_
#define GALOISRING_BASIS_HPP_

#include <NTL/ZZ_pE.h>

#include <vector>

// Shared basis-carrying data for GR(p^k, r) APIs.
//
// WP0 intentionally exposes only the storage type so public compiler APIs can
// freeze on the final setup/params shape before the generic shared basis
// algebra lands in WP1.
struct GaloisRingBasisData {
  std::vector<NTL::ZZ_pE> basis;
  std::vector<NTL::ZZ_pE> dual_basis;
};

#endif  // GALOISRING_BASIS_HPP_
