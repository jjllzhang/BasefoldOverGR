#ifndef GALOISRING_BASIS_HPP_
#define GALOISRING_BASIS_HPP_

#include <NTL/ZZ_p.h>
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

std::vector<NTL::ZZ_pE> BuildPolynomialBasis(long dimension);

void ValidateBasisShapeOrThrow(const std::vector<NTL::ZZ_pE> &basis,
                               const char *label,
                               const char *func_name);

void ValidateBasisDataOrThrow(const GaloisRingBasisData &basis,
                              const char *label,
                              const char *func_name);

bool IsBasisOverBaseRing(const std::vector<NTL::ZZ_pE> &basis);

std::vector<NTL::ZZ_pE> BuildDualBasisOrThrow(
    const std::vector<NTL::ZZ_pE> &basis);

bool IsBaseRingConstant(const NTL::ZZ_pE &value);

NTL::ZZ_p TraceToBaseRing(const NTL::ZZ_pE &element);

std::vector<NTL::ZZ_p> RecoverBasisCoordsOrThrow(
    const GaloisRingBasisData &basis, const NTL::ZZ_pE &element,
    const char *func_name);

NTL::ZZ_pE ComposeFromBasisCoordsOrThrow(
    const std::vector<NTL::ZZ_pE> &basis,
    const std::vector<NTL::ZZ_p> &coords, const char *func_name);

#endif  // GALOISRING_BASIS_HPP_
