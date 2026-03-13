#ifndef GALOISRING_FROBENIUSBASIS_HPP_
#define GALOISRING_FROBENIUSBASIS_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <vector>

// Frobenius/normal-basis helpers for GR(p^k, r).
//
// Context contract:
//   - The caller must have initialized ZZ_p with modulus p^k.
//   - The caller must have initialized ZZ_pE with an extension polynomial of
//     degree r matching extension_modulus.
//   - These helpers validate the active context and operate inside it; they do
//     not own context initialization.

struct FrobeniusBasisParams {
  NTL::ZZ p;
  long k = 0;
  long r = 0;
  long teichmuller_generator_max_trials = 1024;
  long affine_search_limit = 4;
};

struct NormalBasisData {
  std::vector<NTL::ZZ_pE> beta;
  std::vector<NTL::ZZ_pE> alpha;
};

struct FrobeniusBasisData {
  FrobeniusBasisParams params;
  NTL::ZZ modulus;
  NTL::ZZ_pX extension_modulus;
  NTL::ZZ_pE teichmuller_generator;
  NormalBasisData normal_basis;
};

void ValidateFrobeniusBasisContextsOrThrow(
    const FrobeniusBasisParams &params, const NTL::ZZ_pX &extension_modulus);

FrobeniusBasisData BuildFrobeniusBasisOrThrow(
    const FrobeniusBasisParams &params, const NTL::ZZ_pX &extension_modulus);

NormalBasisData FindNormalBasisOrThrow(
    const FrobeniusBasisParams &params,
    const NTL::ZZ_pE &teichmuller_generator);

std::vector<NTL::ZZ_pE> BuildDualBasisOrThrow(
    const std::vector<NTL::ZZ_pE> &basis);

void ValidateNormalBasisOrThrow(const NormalBasisData &normal_basis);

bool IsBaseRingConstant(const NTL::ZZ_pE &value);

NTL::ZZ_p TraceToBaseRing(const NTL::ZZ_pE &element);

std::vector<NTL::ZZ_p> RecoverNormalBasisCoords(
    const NormalBasisData &normal_basis, const NTL::ZZ_pE &element);

NTL::ZZ_pE ComposeFromNormalBasisCoords(
    const NormalBasisData &normal_basis,
    const std::vector<NTL::ZZ_p> &coords);

NTL::ZZ_pE ApplyFrobeniusTau(const NormalBasisData &normal_basis,
                             const NTL::ZZ_pE &element);

NTL::ZZ_pE ApplyFrobeniusSigma(const NormalBasisData &normal_basis,
                               const NTL::ZZ_pE &element);

#endif  // GALOISRING_FROBENIUSBASIS_HPP_
