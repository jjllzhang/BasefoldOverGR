#include "Compiler/Z2k/FrobeniusPCS.hpp"

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "GaloisRing/Basis.hpp"
#include "PCS/Common/Multilinear.hpp"
#include "PCS/Common/Profile.hpp"
#include "PCS/Common/Transcript.hpp"

using NTL::LogicError;
using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pX;
using std::size_t;
using std::string;
using std::vector;

namespace basefold {

NTL::vec_ZZ_pE PackZ2kTableToFrobeniusGREvalsWithTrustedParamsOrThrow(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table);

namespace {

struct PackedCommitInputs {
  NTL::vec_ZZ_pE t_packed_table;
  NTL::vec_ZZ_pE t_packed_monomial_coeffs;
};

struct OuterProveEvalResult {
  FrobeniusPCSOuterEvalProof proof;
  std::vector<FieldElement> rprime_suffix;
};

struct SuffixOrbitCache {
  std::vector<std::vector<ZZ_p>> suffix_coords_by_var;
  std::vector<std::vector<FieldElement>> sigma_points_by_i;
  NTL::vec_ZZ_pE base_eq_table;
  std::vector<ZZ_p> base_eq_coords_by_w_then_j;
};

struct OrbitBuildProfileTargets {
  std::uint64_t *recover_coords_ns = nullptr;
  std::uint64_t *recover_coords_calls = nullptr;
  std::uint64_t *sigma_points_ns = nullptr;
  std::uint64_t *sigma_points_calls = nullptr;
  std::uint64_t *base_eq_table_ns = nullptr;
  std::uint64_t *base_eq_table_calls = nullptr;
  std::uint64_t *base_eq_coords_ns = nullptr;
  std::uint64_t *base_eq_coords_calls = nullptr;
};

long RotateIndex(long index, long shift, long dimension) {
  if (dimension <= 0) {
    LogicError("RotateIndex: dimension must be positive");
  }
  const long normalized_shift =
      ((shift % dimension) + dimension) % dimension;
  const long rotated = index + normalized_shift;
  return (rotated >= dimension) ? (rotated - dimension) : rotated;
}

long Pow2LongOrThrow(long exponent, const char *what) {
  if (exponent < 0) {
    LogicError(what);
  }
  long out = 1;
  for (long i = 0; i < exponent; ++i) {
    if (out > std::numeric_limits<long>::max() / 2) {
      LogicError(what);
    }
    out *= 2;
  }
  return out;
}

long Log2ExactPowerOfTwoZZOrThrow(const ZZ &value, const char *what) {
  if (value <= 0) {
    LogicError(what);
  }
  ZZ tmp = value;
  long out = 0;
  while (tmp > 1) {
    if ((tmp % 2) != 0) {
      LogicError(what);
    }
    tmp /= 2;
    ++out;
  }
  return out;
}

ZZ_pE BaseRingConstant(const ZZ_p &value) {
  ZZ_pX poly;
  if (value != 0) {
    SetCoeff(poly, 0, value);
  }
  ZZ_pE out;
  NTL::conv(out, poly);
  return out;
}

ZZ_pE BaseRingConstant(long value) {
  return BaseRingConstant(NTL::to_ZZ_p(value));
}

FieldElement ComposeFromCoordsWithBasisRow(
    const std::vector<FieldElement> &basis_row,
    const std::vector<ZZ_p> &coords, const char *func_name) {
  if (basis_row.size() != coords.size()) {
    LogicError((string(func_name) +
                ": basis row length must equal coordinate length")
                   .c_str());
  }
  FieldElement out = FieldElement(0);
  for (long i = 0; i < static_cast<long>(coords.size()); ++i) {
    out += BaseRingConstant(coords[static_cast<std::size_t>(i)]) *
           basis_row[static_cast<std::size_t>(i)];
  }
  return out;
}

FrobeniusPCSPrecomputedTables BuildFrobeniusPCSPrecomputedTables(
    const NormalBasisData &normal_basis) {
  const long basis_dimension =
      static_cast<long>(normal_basis.beta.size());
  if (basis_dimension <= 0 ||
      static_cast<long>(normal_basis.alpha.size()) != basis_dimension) {
    LogicError("BuildFrobeniusPCSPrecomputedTables: invalid normal basis shape");
  }

  FrobeniusPCSPrecomputedTables out;
  out.tau_basis_rows.resize(static_cast<std::size_t>(basis_dimension));
  out.sigma_basis_rows.resize(static_cast<std::size_t>(basis_dimension));
  for (long power = 0; power < basis_dimension; ++power) {
    std::vector<FieldElement> &tau_row =
        out.tau_basis_rows[static_cast<std::size_t>(power)];
    std::vector<FieldElement> &sigma_row =
        out.sigma_basis_rows[static_cast<std::size_t>(power)];
    tau_row.resize(static_cast<std::size_t>(basis_dimension));
    sigma_row.resize(static_cast<std::size_t>(basis_dimension));
    for (long j = 0; j < basis_dimension; ++j) {
      tau_row[static_cast<std::size_t>(j)] =
          normal_basis.beta[static_cast<std::size_t>(
              RotateIndex(j, power, basis_dimension))];
      sigma_row[static_cast<std::size_t>(j)] =
          normal_basis.beta[static_cast<std::size_t>(
              RotateIndex(j, -power, basis_dimension))];
    }
  }

  out.tau_alpha_by_u_then_i.resize(static_cast<std::size_t>(basis_dimension));
  for (long u = 0; u < basis_dimension; ++u) {
    const std::vector<ZZ_p> alpha_coords =
        ::RecoverNormalBasisCoords(normal_basis,
                                   normal_basis.alpha[static_cast<std::size_t>(u)]);
    std::vector<FieldElement> &row =
        out.tau_alpha_by_u_then_i[static_cast<std::size_t>(u)];
    row.resize(static_cast<std::size_t>(basis_dimension));
    for (long power = 0; power < basis_dimension; ++power) {
      row[static_cast<std::size_t>(power)] = ComposeFromCoordsWithBasisRow(
          out.tau_basis_rows[static_cast<std::size_t>(power)], alpha_coords,
          "BuildFrobeniusPCSPrecomputedTables");
    }
  }
  return out;
}

void ValidatePrecomputedTablesOrThrow(const FrobeniusPCSParams &params) {
  const long basis_dimension =
      static_cast<long>(params.basis_data.normal_basis.beta.size());
  if (static_cast<long>(params.precomputed.tau_basis_rows.size()) !=
          basis_dimension ||
      static_cast<long>(params.precomputed.sigma_basis_rows.size()) !=
          basis_dimension ||
      static_cast<long>(params.precomputed.tau_alpha_by_u_then_i.size()) !=
          basis_dimension) {
    LogicError("ValidatePrecomputedTablesOrThrow: precomputed table row count mismatch");
  }
  for (long power = 0; power < basis_dimension; ++power) {
    if (static_cast<long>(
            params.precomputed.tau_basis_rows[static_cast<std::size_t>(power)]
                .size()) != basis_dimension ||
        static_cast<long>(
            params.precomputed.sigma_basis_rows[static_cast<std::size_t>(power)]
                .size()) != basis_dimension) {
      LogicError("ValidatePrecomputedTablesOrThrow: precomputed basis-row width mismatch");
    }
  }
  for (long u = 0; u < basis_dimension; ++u) {
    const std::vector<FieldElement> &row =
        params.precomputed.tau_alpha_by_u_then_i[static_cast<std::size_t>(u)];
    if (static_cast<long>(row.size()) != basis_dimension) {
      LogicError("ValidatePrecomputedTablesOrThrow: tau_alpha row width mismatch");
    }
    if (row[0] != params.basis_data.normal_basis.alpha[static_cast<std::size_t>(u)]) {
      LogicError("ValidatePrecomputedTablesOrThrow: tau_alpha power-0 entry mismatch");
    }
  }
  if (params.precomputed.tau_basis_rows[0] != params.basis_data.normal_basis.beta ||
      params.precomputed.sigma_basis_rows[0] !=
          params.basis_data.normal_basis.beta) {
    LogicError("ValidatePrecomputedTablesOrThrow: power-0 basis rows must equal beta");
  }
}

void ValidateBaseRingConstantOrThrow(const ZZ_pE &value, const char *label,
                                     long index, const char *func_name) {
  if (!::IsBaseRingConstant(value)) {
    LogicError((string(func_name) + ": " + label + "[" +
                std::to_string(index) + "] must be a base-ring constant")
                   .c_str());
  }
}

void ValidateBaseRingVectorOrThrow(const NTL::vec_ZZ_pE &values,
                                   const char *label,
                                   const char *func_name) {
  for (long i = 0; i < values.length(); ++i) {
    ValidateBaseRingConstantOrThrow(values[i], label, i, func_name);
  }
}

NTL::vec_ZZ_pE BooleanHypercubeTableToMonomialCoeffsInternal(
    const NTL::vec_ZZ_pE &table_values) {
  const long n = table_values.length();
  if (n <= 0 || (n & (n - 1)) != 0) {
    LogicError("BooleanHypercubeTableToMonomialCoeffsInternal: table length must be a power of two");
  }
  long dimension = 0;
  for (long tmp = n; tmp > 1; tmp >>= 1) {
    ++dimension;
  }

  NTL::vec_ZZ_pE coeffs = table_values;
  for (long bit = 0; bit < dimension; ++bit) {
    const long step = 1L << bit;
    for (long mask = 0; mask < n; ++mask) {
      if (mask & step) {
        coeffs[mask] -= coeffs[mask ^ step];
      }
    }
  }
  return coeffs;
}

PackedCommitInputs BuildPackedCommitInputs(const FrobeniusPCSParams &params,
                                           const NTL::vec_ZZ_pE &t_table) {
  PackedCommitInputs out;
  out.t_packed_table =
      PackZ2kTableToFrobeniusGREvalsWithTrustedParamsOrThrow(params, t_table);
  out.t_packed_monomial_coeffs =
      BooleanHypercubeTableToMonomialCoeffsInternal(out.t_packed_table);
  return out;
}

std::vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  std::vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] = BaseRingConstant((index >> i) & 1L);
  }
  return point;
}

HashTranscript MakeFrobeniusTranscript() {
  HashTranscriptConfig config;
  config.domain_separator = "FrobeniusPCS/v1";
  config.byte_order = TranscriptByteOrder::kLittleEndian;
  config.error_prefix = "FrobeniusHashTranscript";
  return HashTranscript(config);
}

void AbsorbPublicInput(HashTranscript &transcript,
                       const MerkleRoot &commitment,
                       const std::vector<FieldElement> &z,
                       const FieldElement &claimed_s) {
  transcript.AbsorbDigest(commitment);
  for (const FieldElement &zi : z) {
    transcript.AbsorbFieldElement(zi);
  }
  transcript.AbsorbFieldElement(claimed_s);
}

void ValidateEvalInputsOrThrow(const FrobeniusPCSParams &params,
                               const NTL::vec_ZZ_pE &t_table,
                               const std::vector<FieldElement> &z,
                               long num_queries, const char *func_name) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const long expected_t_length = Pow2LongOrThrow(
      params.ell, (std::string(func_name) + ": ell is too large for long").c_str());
  if (t_table.length() != expected_t_length) {
    LogicError((std::string(func_name) + ": t_table length must equal 2^ell")
                   .c_str());
  }
  ValidateBaseRingVectorOrThrow(t_table, "t_table", func_name);
  if (static_cast<long>(z.size()) != params.ell) {
    LogicError((std::string(func_name) + ": z dimension must equal ell").c_str());
  }
  if (num_queries < 0) {
    LogicError((std::string(func_name) + ": num_queries must be non-negative")
                   .c_str());
  }
}

void ValidateCommitArtifactsOrThrow(const FrobeniusPCSParams &params,
                                    const FrobeniusPCSCommitArtifacts &artifacts,
                                    const char *func_name) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const long expected_packed_length = Pow2LongOrThrow(
      params.ell_prime,
      (std::string(func_name) + ": ell_prime is too large for long").c_str());
  if (artifacts.t_packed_table.length() != expected_packed_length) {
    LogicError((std::string(func_name) +
                ": t_packed_table length must equal 2^(ell-kappa)")
                   .c_str());
  }
  if (artifacts.t_packed_monomial_coeffs.length() != expected_packed_length) {
    LogicError((std::string(func_name) +
                ": t_packed_monomial_coeffs length must equal 2^(ell-kappa)")
                   .c_str());
  }
  if (artifacts.commitment != artifacts.backend_commit_artifacts.commitment) {
    LogicError((std::string(func_name) +
                ": commitment must match backend_commit_artifacts.commitment")
                   .c_str());
  }
}

void ValidateOuterCommitArtifactsOrThrow(
    const FrobeniusPCSParams &params,
    const FrobeniusPCSOuterCommitArtifacts &artifacts,
    const char *func_name) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const long expected_packed_length = Pow2LongOrThrow(
      params.ell_prime,
      (std::string(func_name) + ": ell_prime is too large for long").c_str());
  if (artifacts.t_packed_table.length() != expected_packed_length) {
    LogicError((std::string(func_name) +
                ": t_packed_table length must equal 2^(ell-kappa)")
                   .c_str());
  }
  if (artifacts.t_packed_monomial_coeffs.length() != expected_packed_length) {
    LogicError((std::string(func_name) +
                ": t_packed_monomial_coeffs length must equal 2^(ell-kappa)")
                   .c_str());
  }
}

bool HasExpectedEvalProofShape(const FrobeniusPCSParams &params,
                               const FrobeniusPCSEvalProof &proof) {
  const long basis_dimension =
      static_cast<long>(params.basis_data.normal_basis.beta.size());
  return static_cast<long>(proof.s_by_i.size()) == basis_dimension &&
         static_cast<long>(proof.h_by_level.size()) == params.ell_prime;
}

bool HasExpectedOuterEvalProofShape(const FrobeniusPCSParams &params,
                                    const FrobeniusPCSOuterEvalProof &proof) {
  const long basis_dimension =
      static_cast<long>(params.basis_data.normal_basis.beta.size());
  return static_cast<long>(proof.s_by_i.size()) == basis_dimension &&
         static_cast<long>(proof.h_by_level.size()) == params.ell_prime;
}

bool HasCompatibleBackendEvalSubproof(const FrobeniusPCSParams &params,
                                      const Z2kPCSBackendEvalProof &proof) {
  return proof.vtable == params.backend.vtable && proof.payload &&
         proof.params_owner &&
         proof.params_owner.get() == params.backend.params.get();
}

std::vector<FieldElement> SlicePoint(const std::vector<FieldElement> &z,
                                     long begin, long count) {
  return std::vector<FieldElement>(z.begin() + begin, z.begin() + begin + count);
}

std::vector<std::vector<ZZ_p>> RecoverPointNormalBasisCoords(
    const FrobeniusPCSParams &params, const std::vector<FieldElement> &point) {
  std::vector<std::vector<ZZ_p>> coords_by_var(point.size());
  for (long i = 0; i < static_cast<long>(point.size()); ++i) {
    coords_by_var[static_cast<std::size_t>(i)] =
        ::RecoverNormalBasisCoords(params.basis_data.normal_basis,
                                   point[static_cast<std::size_t>(i)]);
  }
  return coords_by_var;
}

std::vector<FieldElement> ComposePointWithBasisRow(
    const std::vector<FieldElement> &basis_row,
    const std::vector<std::vector<ZZ_p>> &coords_by_var,
    const char *func_name) {
  std::vector<FieldElement> out(coords_by_var.size(), FieldElement(0));
  for (long i = 0; i < static_cast<long>(coords_by_var.size()); ++i) {
    out[static_cast<std::size_t>(i)] = ComposeFromCoordsWithBasisRow(
        basis_row, coords_by_var[static_cast<std::size_t>(i)], func_name);
  }
  return out;
}

const ZZ_p *BaseEqCoordsRowOrThrow(const SuffixOrbitCache &orbit_cache, long w,
                                   long basis_dimension,
                                   const char *func_name) {
  if (basis_dimension < 0) {
    LogicError((string(func_name) +
                ": basis dimension must be non-negative")
                   .c_str());
  }
  const size_t row_offset =
      static_cast<size_t>(w) * static_cast<size_t>(basis_dimension);
  if (row_offset + static_cast<size_t>(basis_dimension) >
      orbit_cache.base_eq_coords_by_w_then_j.size()) {
    LogicError((string(func_name) +
                ": base equality-coordinate row is out of bounds")
                   .c_str());
  }
  return orbit_cache.base_eq_coords_by_w_then_j.data() + row_offset;
}

SuffixOrbitCache BuildSuffixOrbitCache(const FrobeniusPCSParams &params,
                                       const std::vector<FieldElement> &r_suffix,
                                       bool build_eq_tables,
                                       const OrbitBuildProfileTargets *profile_targets =
                                           nullptr) {
  if (static_cast<long>(r_suffix.size()) != params.ell_prime) {
    LogicError("BuildSuffixOrbitCache: r_suffix dimension must equal ell_prime");
  }

  SuffixOrbitCache out;
  {
    ScopedTimer timer(profile_targets ? profile_targets->recover_coords_ns
                                      : nullptr,
                      profile_targets ? profile_targets->recover_coords_calls
                                      : nullptr);
    out.suffix_coords_by_var = RecoverPointNormalBasisCoords(params, r_suffix);
  }
  const long basis_dimension =
      static_cast<long>(params.precomputed.sigma_basis_rows.size());
  const long num_w = build_eq_tables
                         ? Pow2LongOrThrow(
                               params.ell_prime,
                               "BuildSuffixOrbitCache: ell_prime is too large for long")
                         : 0;
  out.sigma_points_by_i.resize(static_cast<std::size_t>(basis_dimension));
  if (build_eq_tables) {
    {
      ScopedTimer timer(profile_targets ? profile_targets->base_eq_table_ns
                                        : nullptr,
                        profile_targets ? profile_targets->base_eq_table_calls
                                        : nullptr);
      out.base_eq_table = EqualityTableFromPoint(r_suffix);
    }
    if (out.base_eq_table.length() != num_w) {
      LogicError("BuildSuffixOrbitCache: base equality-table width mismatch");
    }
    out.base_eq_coords_by_w_then_j.resize(
        static_cast<size_t>(num_w) * static_cast<size_t>(basis_dimension));
    {
      ScopedTimer timer(profile_targets ? profile_targets->base_eq_coords_ns
                                        : nullptr,
                        profile_targets ? profile_targets->base_eq_coords_calls
                                        : nullptr);
      for (long w = 0; w < num_w; ++w) {
        const std::vector<ZZ_p> coords =
            ::RecoverNormalBasisCoords(
                params.basis_data.normal_basis,
                out.base_eq_table[static_cast<std::size_t>(w)]);
        if (static_cast<long>(coords.size()) != basis_dimension) {
          LogicError(
              "BuildSuffixOrbitCache: base equality-coordinate width mismatch");
        }
        const size_t row_offset =
            static_cast<size_t>(w) * static_cast<size_t>(basis_dimension);
        for (long j = 0; j < basis_dimension; ++j) {
          out.base_eq_coords_by_w_then_j[row_offset +
                                         static_cast<size_t>(j)] =
              coords[static_cast<std::size_t>(j)];
        }
      }
    }
  }
  for (long i = 0; i < basis_dimension; ++i) {
    {
      ScopedTimer timer(profile_targets ? profile_targets->sigma_points_ns
                                        : nullptr,
                        profile_targets ? profile_targets->sigma_points_calls
                                        : nullptr);
      out.sigma_points_by_i[static_cast<std::size_t>(i)] =
          ComposePointWithBasisRow(
              params.precomputed.sigma_basis_rows[static_cast<std::size_t>(i)],
              out.suffix_coords_by_var, "BuildSuffixOrbitCache");
    }
  }
  return out;
}

std::vector<FieldElement> ComputeIndexedTauPowers(
    const FrobeniusPCSParams &params,
    const std::vector<FieldElement> &elements) {
  const long basis_dimension =
      static_cast<long>(params.precomputed.tau_basis_rows.size());
  if (static_cast<long>(elements.size()) != basis_dimension) {
    LogicError(
        "ComputeIndexedTauPowers: element count must equal basis dimension");
  }

  std::vector<FieldElement> out(static_cast<std::size_t>(basis_dimension),
                                FieldElement(0));
  for (long i = 0; i < basis_dimension; ++i) {
    const std::vector<ZZ_p> coords = ::RecoverNormalBasisCoords(
        params.basis_data.normal_basis, elements[static_cast<std::size_t>(i)]);
    out[static_cast<std::size_t>(i)] = ComposeFromCoordsWithBasisRow(
        params.precomputed.tau_basis_rows[static_cast<std::size_t>(i)], coords,
        "ComputeIndexedTauPowers");
  }
  return out;
}

std::vector<FieldElement> RecoverPartialEvaluationsFromSByI(
    const FrobeniusPCSParams &params,
    const std::vector<FieldElement> &s_by_i) {
  const long basis_dimension = static_cast<long>(
      params.precomputed.tau_alpha_by_u_then_i.size());
  if (static_cast<long>(s_by_i.size()) != basis_dimension) {
    LogicError(
        "RecoverPartialEvaluationsFromSByI: s_by_i size must equal basis dimension");
  }

  const std::vector<FieldElement> tau_s_by_i =
      ComputeIndexedTauPowers(params, s_by_i);
  std::vector<FieldElement> partials(static_cast<std::size_t>(basis_dimension),
                                     FieldElement(0));
  for (long u = 0; u < basis_dimension; ++u) {
    FieldElement acc = FieldElement(0);
    for (long i = 0; i < basis_dimension; ++i) {
      acc += params.precomputed.tau_alpha_by_u_then_i
                 [static_cast<std::size_t>(u)][static_cast<std::size_t>(i)] *
             tau_s_by_i[static_cast<std::size_t>(i)];
    }
    partials[static_cast<std::size_t>(u)] = acc;
  }
  return partials;
}

std::vector<FieldElement> ComputeDirectPartialEvaluations(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z_suffix) {
  const long basis_dimension =
      static_cast<long>(params.basis_data.normal_basis.beta.size());
  const long num_w = Pow2LongOrThrow(
      params.ell_prime,
      "ComputeDirectPartialEvaluations: ell_prime is too large for long");
  std::vector<FieldElement> partials(static_cast<std::size_t>(basis_dimension),
                                     FieldElement(0));
  for (long u = 0; u < basis_dimension; ++u) {
    NTL::vec_ZZ_pE slice;
    slice.SetLength(num_w);
    for (long w = 0; w < num_w; ++w) {
      slice[w] = t_table[u + w * basis_dimension];
    }
    const NTL::vec_ZZ_pE slice_monomial =
        BooleanHypercubeTableToMonomialCoeffsInternal(slice);
    partials[static_cast<std::size_t>(u)] =
        EvalMultilinearMonomialCoeffs(slice_monomial, z_suffix);
  }
  return partials;
}

FieldElement RecombineClaimFromPartials(const std::vector<FieldElement> &partials,
                                        const std::vector<FieldElement> &z_prefix) {
  FieldElement acc = FieldElement(0);
  for (long u = 0; u < static_cast<long>(partials.size()); ++u) {
    acc += partials[static_cast<std::size_t>(u)] *
           EqPolynomial(z_prefix, BooleanPointFromIndex(u,
                                                        static_cast<long>(z_prefix.size())));
  }
  return acc;
}

std::vector<FieldElement> ComputeSByIFromBaseEqCoords(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_packed_table,
    const SuffixOrbitCache &orbit_cache) {
  const long num_w = orbit_cache.base_eq_table.length();
  if (t_packed_table.length() != num_w) {
    LogicError(
        "ComputeSByIFromBaseEqCoords: packed table width mismatch");
  }
  const long basis_dimension =
      static_cast<long>(params.precomputed.sigma_basis_rows.size());
  if (static_cast<long>(orbit_cache.base_eq_coords_by_w_then_j.size()) !=
      num_w * basis_dimension) {
    LogicError(
        "ComputeSByIFromBaseEqCoords: base equality-coordinate shape mismatch");
  }
  std::vector<FieldElement> a_by_j(static_cast<std::size_t>(basis_dimension),
                                   FieldElement(0));
  for (long w = 0; w < num_w; ++w) {
    const ZZ_p *coords = BaseEqCoordsRowOrThrow(
        orbit_cache, w, basis_dimension, "ComputeSByIFromBaseEqCoords");
    const FieldElement &packed_eval = t_packed_table[w];
    for (long j = 0; j < basis_dimension; ++j) {
      if (coords[j] != 0) {
        a_by_j[static_cast<std::size_t>(j)] +=
            packed_eval * BaseRingConstant(coords[j]);
      }
    }
  }

  std::vector<FieldElement> s_by_i(static_cast<std::size_t>(basis_dimension),
                                   FieldElement(0));
  for (long i = 0; i < basis_dimension; ++i) {
    FieldElement acc = FieldElement(0);
    const std::vector<FieldElement> &sigma_basis_row =
        params.precomputed.sigma_basis_rows[static_cast<std::size_t>(i)];
    if (static_cast<long>(sigma_basis_row.size()) != basis_dimension) {
      LogicError(
          "ComputeSByIFromBaseEqCoords: sigma basis-row width mismatch");
    }
    for (long j = 0; j < basis_dimension; ++j) {
      acc += a_by_j[static_cast<std::size_t>(j)] *
             sigma_basis_row[static_cast<std::size_t>(j)];
    }
    s_by_i[static_cast<std::size_t>(i)] = acc;
  }
  return s_by_i;
}

std::vector<FieldElement> ComputeSByI(
    const FrobeniusPCSParams &params,
    const NTL::vec_ZZ_pE &t_packed_table,
    const NTL::vec_ZZ_pE &t_packed_monomial_coeffs,
    const SuffixOrbitCache &orbit_cache) {
  const long basis_dimension =
      static_cast<long>(orbit_cache.sigma_points_by_i.size());
  if (!orbit_cache.base_eq_coords_by_w_then_j.empty()) {
    if (static_cast<long>(params.precomputed.sigma_basis_rows.size()) !=
        basis_dimension) {
      LogicError("ComputeSByI: sigma basis-row count mismatch");
    }
    return ComputeSByIFromBaseEqCoords(params, t_packed_table, orbit_cache);
  }

  std::vector<FieldElement> s_by_i(static_cast<std::size_t>(basis_dimension),
                                   FieldElement(0));
  for (long i = 0; i < basis_dimension; ++i) {
    s_by_i[static_cast<std::size_t>(i)] =
        EvalMultilinearMonomialCoeffs(
            t_packed_monomial_coeffs,
            orbit_cache.sigma_points_by_i[static_cast<std::size_t>(i)]);
  }
  return s_by_i;
}

NTL::vec_ZZ_pE BuildBatchedGTable(const FrobeniusPCSParams &params,
                                  const NTL::vec_ZZ_pE &lambda_by_i,
                                  const SuffixOrbitCache &orbit_cache,
                                  long ell_prime) {
  const long num_w = orbit_cache.base_eq_table.length();
  const long basis_dimension =
      static_cast<long>(params.precomputed.sigma_basis_rows.size());
  if (lambda_by_i.length() != basis_dimension) {
    LogicError("BuildBatchedGTable: prefix equality table length mismatch");
  }
  const long expected_num_w = Pow2LongOrThrow(
      ell_prime, "BuildBatchedGTable: ell_prime is too large for long");
  if (num_w != expected_num_w) {
    LogicError("BuildBatchedGTable: base equality-table height mismatch");
  }
  if (static_cast<long>(orbit_cache.base_eq_coords_by_w_then_j.size()) !=
      num_w * basis_dimension) {
    LogicError("BuildBatchedGTable: base equality-coordinate shape mismatch");
  }

  std::vector<FieldElement> mu_by_j(static_cast<std::size_t>(basis_dimension),
                                    FieldElement(0));
  for (long i = 0; i < basis_dimension; ++i) {
    const std::vector<FieldElement> &sigma_basis_row =
        params.precomputed.sigma_basis_rows[static_cast<std::size_t>(i)];
    if (static_cast<long>(sigma_basis_row.size()) != basis_dimension) {
      LogicError("BuildBatchedGTable: sigma basis-row width mismatch");
    }
    const FieldElement &lambda_i = lambda_by_i[i];
    for (long j = 0; j < basis_dimension; ++j) {
      mu_by_j[static_cast<std::size_t>(j)] +=
          lambda_i * sigma_basis_row[static_cast<std::size_t>(j)];
    }
  }

  NTL::vec_ZZ_pE out;
  out.SetLength(num_w);
  for (long w = 0; w < num_w; ++w) {
    const ZZ_p *coords =
        BaseEqCoordsRowOrThrow(orbit_cache, w, basis_dimension,
                               "BuildBatchedGTable");
    NTL::clear(out[w]);
    for (long j = 0; j < basis_dimension; ++j) {
      if (coords[j] != 0) {
        out[w] += BaseRingConstant(coords[j]) *
                  mu_by_j[static_cast<std::size_t>(j)];
      }
    }
  }
  return out;
}

FieldElement ComputeInitialBatchedClaim(const std::vector<FieldElement> &s_by_i,
                                        const NTL::vec_ZZ_pE &lambda_by_i) {
  if (lambda_by_i.length() != static_cast<long>(s_by_i.size())) {
    LogicError(
        "ComputeInitialBatchedClaim: prefix equality table length mismatch");
  }
  FieldElement acc = FieldElement(0);
  for (long i = 0; i < static_cast<long>(s_by_i.size()); ++i) {
    acc += s_by_i[static_cast<std::size_t>(i)] * lambda_by_i[i];
  }
  return acc;
}

FieldElement ComputeFinalGStar(const NTL::vec_ZZ_pE &lambda_by_i,
                               const SuffixOrbitCache &orbit_cache,
                               const std::vector<FieldElement> &rprime_suffix) {
  const long basis_dimension =
      static_cast<long>(orbit_cache.sigma_points_by_i.size());
  if (lambda_by_i.length() != basis_dimension) {
    LogicError("ComputeFinalGStar: prefix equality table length mismatch");
  }

  FieldElement g_star = FieldElement(0);
  for (long i = 0; i < basis_dimension; ++i) {
    g_star += lambda_by_i[i] *
              EqPolynomial(
                  rprime_suffix,
                  orbit_cache.sigma_points_by_i[static_cast<std::size_t>(i)]);
  }
  return g_star;
}

FieldElement ComputeOriginalEvaluation(const NTL::vec_ZZ_pE &t_table,
                                       const std::vector<FieldElement> &z) {
  const NTL::vec_ZZ_pE t_monomial =
      BooleanHypercubeTableToMonomialCoeffsInternal(t_table);
  return EvalMultilinearMonomialCoeffs(t_monomial, z);
}

OuterProveEvalResult ProveOuterEvalFromCommitArtifactsInternal(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries,
    const FrobeniusPCSOuterCommitArtifacts &commit_artifacts,
    bool checked_path, const char *func_name) {
  Profile *prof = ActiveProfile();
  ScopedTimer total_timer(prof ? &prof->frobenius_outer_prove_total_ns
                               : nullptr,
                          prof ? &prof->frobenius_outer_prove_total_calls
                               : nullptr);
  if (checked_path) {
    ValidateEvalInputsOrThrow(params, t_table, z, num_queries, func_name);
    ValidateOuterCommitArtifactsOrThrow(params, commit_artifacts, func_name);
  }

  if (checked_path) {
    const FieldElement direct_eval = ComputeOriginalEvaluation(t_table, z);
    if (direct_eval != claimed_s) {
      LogicError(
          (std::string(func_name) + ": claimed_s must equal t(z)").c_str());
    }
  }

  const std::vector<FieldElement> z_prefix = SlicePoint(z, 0, params.kappa);
  const std::vector<FieldElement> r_suffix =
      SlicePoint(z, params.kappa, params.ell_prime);
  const SuffixOrbitCache suffix_orbit = [&] {
    OrbitBuildProfileTargets orbit_profile;
    if (prof != nullptr) {
      orbit_profile.recover_coords_ns =
          &prof->frobenius_outer_prove_orbit_recover_coords_ns;
      orbit_profile.recover_coords_calls =
          &prof->frobenius_outer_prove_orbit_recover_coords_calls;
      orbit_profile.sigma_points_ns =
          &prof->frobenius_outer_prove_orbit_sigma_points_ns;
      orbit_profile.sigma_points_calls =
          &prof->frobenius_outer_prove_orbit_sigma_points_calls;
      orbit_profile.base_eq_table_ns =
          &prof->frobenius_outer_prove_orbit_base_eq_table_ns;
      orbit_profile.base_eq_table_calls =
          &prof->frobenius_outer_prove_orbit_base_eq_table_calls;
      orbit_profile.base_eq_coords_ns =
          &prof->frobenius_outer_prove_orbit_base_eq_coords_ns;
      orbit_profile.base_eq_coords_calls =
          &prof->frobenius_outer_prove_orbit_base_eq_coords_calls;
    }
    ScopedTimer timer(prof ? &prof->frobenius_outer_prove_orbit_build_ns
                           : nullptr,
                      prof ? &prof->frobenius_outer_prove_orbit_build_calls
                           : nullptr);
    return BuildSuffixOrbitCache(params, r_suffix, params.ell_prime > 0,
                                 prof ? &orbit_profile : nullptr);
  }();

  OuterProveEvalResult out;
  {
    ScopedTimer timer(prof ? &prof->frobenius_outer_prove_compute_s_ns
                           : nullptr,
                      prof ? &prof->frobenius_outer_prove_compute_s_calls
                           : nullptr);
    out.proof.s_by_i = ComputeSByI(params, commit_artifacts.t_packed_table,
                                   commit_artifacts.t_packed_monomial_coeffs,
                                   suffix_orbit);
  }
  out.proof.h_by_level.resize(static_cast<std::size_t>(params.ell_prime));

  const std::vector<FieldElement> recovered_partials = [&] {
    ScopedTimer timer(
        prof ? &prof->frobenius_outer_prove_recover_partials_ns : nullptr,
        prof ? &prof->frobenius_outer_prove_recover_partials_calls : nullptr);
    return RecoverPartialEvaluationsFromSByI(params, out.proof.s_by_i);
  }();
  if (checked_path) {
    const std::vector<FieldElement> direct_partials =
        ComputeDirectPartialEvaluations(params, t_table, r_suffix);
    if (recovered_partials != direct_partials) {
      LogicError((std::string(func_name) +
                  ": recovered partial evaluations do not match Protocol 2 reconstruction")
                     .c_str());
    }
    if (RecombineClaimFromPartials(recovered_partials, z_prefix) != claimed_s) {
      LogicError((std::string(func_name) +
                  ": Equality Check 1 failed on honest witness")
                     .c_str());
    }
  }

  HashTranscript transcript = MakeFrobeniusTranscript();
  std::vector<FieldElement> rprime_prefix(
      static_cast<std::size_t>(params.kappa));
  NTL::vec_ZZ_pE lambda_by_i;
  FieldElement initial_claim = FieldElement(0);
  NTL::vec_ZZ_pE g_table;
  {
    ScopedTimer timer(prof ? &prof->frobenius_outer_prove_batch_prep_ns
                           : nullptr,
                      prof ? &prof->frobenius_outer_prove_batch_prep_calls
                           : nullptr);
    {
      ScopedTimer sub_timer(prof ? &prof->frobenius_outer_prove_batch_absorb_ns
                                 : nullptr,
                            prof ? &prof->frobenius_outer_prove_batch_absorb_calls
                                 : nullptr);
      AbsorbPublicInput(transcript, commitment, z, claimed_s);
      for (const FieldElement &s_i : out.proof.s_by_i) {
        transcript.AbsorbFieldElement(s_i);
      }
    }

    {
      ScopedTimer sub_timer(
          prof ? &prof->frobenius_outer_prove_batch_prefix_challenge_ns
               : nullptr,
          prof ? &prof->frobenius_outer_prove_batch_prefix_challenge_calls
               : nullptr);
      for (long i = 0; i < params.kappa; ++i) {
        rprime_prefix[static_cast<std::size_t>(i)] =
            transcript.ChallengeFieldElement("rprime/prefix/" +
                                             std::to_string(i));
      }
      lambda_by_i = EqualityTableFromPoint(rprime_prefix);
    }
    {
      ScopedTimer sub_timer(
          prof ? &prof->frobenius_outer_prove_batch_initial_claim_ns
               : nullptr,
          prof ? &prof->frobenius_outer_prove_batch_initial_claim_calls
               : nullptr);
      initial_claim = ComputeInitialBatchedClaim(out.proof.s_by_i, lambda_by_i);
    }
    if (params.ell_prime > 0) {
      ScopedTimer sub_timer(
          prof ? &prof->frobenius_outer_prove_batch_build_g_table_ns
               : nullptr,
          prof ? &prof->frobenius_outer_prove_batch_build_g_table_calls
               : nullptr);
      g_table = BuildBatchedGTable(params, lambda_by_i, suffix_orbit,
                                   params.ell_prime);
    }
  }

  out.rprime_suffix.resize(static_cast<std::size_t>(params.ell_prime));
  if (params.ell_prime > 0) {
    ScopedTimer timer(prof ? &prof->frobenius_outer_prove_sumcheck_ns
                           : nullptr,
                      prof ? &prof->frobenius_outer_prove_sumcheck_calls
                           : nullptr);
    ProductSumcheckProver sumcheck(commit_artifacts.t_packed_table, g_table);
    out.proof.h_by_level[static_cast<std::size_t>(params.ell_prime - 1)] =
        sumcheck.CurrentPolynomial();
    AbsorbQuadraticPoly(
        transcript,
        out.proof.h_by_level[static_cast<std::size_t>(params.ell_prime - 1)]);

    for (long i = params.ell_prime; i-- > 0;) {
      const FieldElement r_i = transcript.ChallengeFieldElement(
          "rprime/suffix/" + std::to_string(i));
      out.rprime_suffix[static_cast<std::size_t>(i)] = r_i;
      sumcheck.ReceiveChallenge(r_i);
      if (i > 0) {
        out.proof.h_by_level[static_cast<std::size_t>(i - 1)] =
            sumcheck.CurrentPolynomial();
        AbsorbQuadraticPoly(
            transcript,
            out.proof.h_by_level[static_cast<std::size_t>(i - 1)]);
      }
    }

    if (checked_path &&
        !CheckProductSumcheckChain(initial_claim, out.proof.h_by_level,
                                   out.rprime_suffix)) {
      LogicError((std::string(func_name) +
                  ": honest product sumcheck chain is inconsistent")
                     .c_str());
    }
  }

  {
    ScopedTimer timer(prof ? &prof->frobenius_outer_prove_final_check_ns
                           : nullptr,
                      prof ? &prof->frobenius_outer_prove_final_check_calls
                           : nullptr);
    out.proof.t_star = EvalMultilinearMonomialCoeffs(
        commit_artifacts.t_packed_monomial_coeffs, out.rprime_suffix);
    if (checked_path) {
      const FieldElement g_star =
          ComputeFinalGStar(lambda_by_i, suffix_orbit, out.rprime_suffix);
      const FieldElement final_sumcheck_claim =
          (params.ell_prime == 0)
              ? initial_claim
              : out.proof.h_by_level[0].Eval(out.rprime_suffix[0]);
      if (final_sumcheck_claim != out.proof.t_star * g_star) {
        LogicError((std::string(func_name) +
                    ": honest Equality Check 3 failed")
                       .c_str());
      }
    }
  }

  return out;
}

FrobeniusPCSEvalProof ProveEvalFromCommitArtifactsInternal(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const FrobeniusPCSCommitArtifacts &commit_artifacts,
    bool checked_path, const char *func_name) {
  if (checked_path) {
    ValidateCommitArtifactsOrThrow(params, commit_artifacts, func_name);
  }

  FrobeniusPCSOuterCommitArtifacts outer_commit_artifacts;
  outer_commit_artifacts.t_packed_table = commit_artifacts.t_packed_table;
  outer_commit_artifacts.t_packed_monomial_coeffs =
      commit_artifacts.t_packed_monomial_coeffs;
  OuterProveEvalResult outer = ProveOuterEvalFromCommitArtifactsInternal(
      params, t_table, commit_artifacts.commitment, z, claimed_s, num_queries,
      outer_commit_artifacts, checked_path, func_name);

  FrobeniusPCSEvalProof proof;
  proof.s_by_i = outer.proof.s_by_i;
  proof.h_by_level = outer.proof.h_by_level;
  proof.t_star = outer.proof.t_star;

  {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->z2k_backend_prove_ns : nullptr,
                      prof ? &prof->z2k_backend_prove_calls : nullptr);
    proof.backend_proof = Z2kPCSBackendProveEval(
        params.backend, commit_artifacts.t_packed_monomial_coeffs,
        outer.rprime_suffix, proof.t_star, num_queries,
        &commit_artifacts.backend_commit_artifacts);
  }
  return proof;
}

bool VerifyOuterEvalAndMaybeRecoverSuffix(
    const FrobeniusPCSParams &params, const MerkleRoot &commitment,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const FrobeniusPCSOuterEvalProof &proof,
    std::vector<FieldElement> *rprime_suffix_out) {
  Profile *prof = ActiveProfile();
  ScopedTimer total_timer(prof ? &prof->frobenius_outer_verify_total_ns
                               : nullptr,
                          prof ? &prof->frobenius_outer_verify_total_calls
                               : nullptr);
  if (static_cast<long>(z.size()) != params.ell) {
    return false;
  }
  if (num_queries < 0) {
    return false;
  }
  if (!HasExpectedOuterEvalProofShape(params, proof)) {
    return false;
  }

  const std::vector<FieldElement> z_prefix = SlicePoint(z, 0, params.kappa);
  const std::vector<FieldElement> r_suffix =
      SlicePoint(z, params.kappa, params.ell_prime);
  const SuffixOrbitCache suffix_orbit = [&] {
    ScopedTimer timer(prof ? &prof->frobenius_outer_verify_orbit_build_ns
                           : nullptr,
                      prof ? &prof->frobenius_outer_verify_orbit_build_calls
                           : nullptr);
    return BuildSuffixOrbitCache(params, r_suffix, /*build_eq_tables=*/false);
  }();

  const std::vector<FieldElement> recovered_partials = [&] {
    ScopedTimer timer(
        prof ? &prof->frobenius_outer_verify_recover_partials_ns : nullptr,
        prof ? &prof->frobenius_outer_verify_recover_partials_calls : nullptr);
    return RecoverPartialEvaluationsFromSByI(params, proof.s_by_i);
  }();
  if (RecombineClaimFromPartials(recovered_partials, z_prefix) != claimed_s) {
    return false;
  }

  HashTranscript transcript = MakeFrobeniusTranscript();
  std::vector<FieldElement> rprime_prefix(
      static_cast<std::size_t>(params.kappa));
  NTL::vec_ZZ_pE lambda_by_i;
  FieldElement initial_claim = FieldElement(0);
  {
    ScopedTimer timer(prof ? &prof->frobenius_outer_verify_prefix_replay_ns
                           : nullptr,
                      prof ? &prof->frobenius_outer_verify_prefix_replay_calls
                           : nullptr);
    AbsorbPublicInput(transcript, commitment, z, claimed_s);
    for (const FieldElement &s_i : proof.s_by_i) {
      transcript.AbsorbFieldElement(s_i);
    }

    for (long i = 0; i < params.kappa; ++i) {
      rprime_prefix[static_cast<std::size_t>(i)] = transcript.ChallengeFieldElement(
          "rprime/prefix/" + std::to_string(i));
    }
    lambda_by_i = EqualityTableFromPoint(rprime_prefix);
    initial_claim = ComputeInitialBatchedClaim(proof.s_by_i, lambda_by_i);
  }

  std::vector<FieldElement> rprime_suffix(
      static_cast<std::size_t>(params.ell_prime));
  if (params.ell_prime > 0) {
    ScopedTimer timer(prof ? &prof->frobenius_outer_verify_sumcheck_replay_ns
                           : nullptr,
                      prof ? &prof->frobenius_outer_verify_sumcheck_replay_calls
                           : nullptr);
    AbsorbQuadraticPoly(
        transcript,
        proof.h_by_level[static_cast<std::size_t>(params.ell_prime - 1)]);
    for (long i = params.ell_prime; i-- > 0;) {
      rprime_suffix[static_cast<std::size_t>(i)] = transcript.ChallengeFieldElement(
          "rprime/suffix/" + std::to_string(i));
      if (i > 0) {
        AbsorbQuadraticPoly(
            transcript,
            proof.h_by_level[static_cast<std::size_t>(i - 1)]);
      }
    }
  }

  if (!CheckProductSumcheckChain(initial_claim, proof.h_by_level,
                                 rprime_suffix)) {
    return false;
  }

  {
    ScopedTimer timer(prof ? &prof->frobenius_outer_verify_final_check_ns
                           : nullptr,
                      prof ? &prof->frobenius_outer_verify_final_check_calls
                           : nullptr);
    const FieldElement g_star =
        ComputeFinalGStar(lambda_by_i, suffix_orbit, rprime_suffix);
    const FieldElement final_sumcheck_claim =
        (params.ell_prime == 0)
            ? initial_claim
            : proof.h_by_level[0].Eval(rprime_suffix[0]);
    if (final_sumcheck_claim != proof.t_star * g_star) {
      return false;
    }
  }

  if (rprime_suffix_out != nullptr) {
    *rprime_suffix_out = std::move(rprime_suffix);
  }
  return true;
}

}  // namespace

void ValidateCurrentZ2kFrobeniusContextOrThrow(
    const ZZ &base_modulus, const ZZ_pX &extension_modulus, long kappa) {
  if (kappa < 1) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: kappa must be >= 1");
  }
  if (base_modulus <= 1) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: base_modulus must be > 1");
  }
  (void)Log2ExactPowerOfTwoZZOrThrow(
      base_modulus,
      "ValidateCurrentZ2kFrobeniusContextOrThrow: base_modulus must be a power of two");
  const long expected_degree = Pow2LongOrThrow(
      kappa,
      "ValidateCurrentZ2kFrobeniusContextOrThrow: kappa is too large for long");
  if (NTL::ZZ_p::modulus() != base_modulus) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: current ZZ_p modulus does not match base_modulus");
  }
  if (!NTL::ZZ_pE::initialized()) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: ZZ_pE context must be initialized");
  }
  if (NTL::ZZ_pE::degree() != expected_degree) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: current ZZ_pE degree must equal 2^kappa");
  }
  if (NTL::deg(extension_modulus) != expected_degree) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: extension_modulus degree must equal 2^kappa");
  }
  if (NTL::ZZ_pE::modulus() != extension_modulus) {
    LogicError(
        "ValidateCurrentZ2kFrobeniusContextOrThrow: current ZZ_pE modulus does not match extension_modulus");
  }
}

void ValidateFrobeniusPCSParamsOrThrow(const FrobeniusPCSParams &params) {
  if (params.kappa < 1) {
    LogicError("ValidateFrobeniusPCSParamsOrThrow: kappa must be >= 1");
  }
  if (params.ell < params.kappa) {
    LogicError("ValidateFrobeniusPCSParamsOrThrow: ell must be >= kappa");
  }

  ValidateCurrentZ2kFrobeniusContextOrThrow(params.base_modulus,
                                            params.extension_modulus,
                                            params.kappa);
  Z2kPCSBackendValidateParamsOrThrow(params.backend);

  const long base_log2 = Log2ExactPowerOfTwoZZOrThrow(
      params.base_modulus,
      "ValidateFrobeniusPCSParamsOrThrow: base_modulus must be a power of two");
  const long expected_degree = Pow2LongOrThrow(
      params.kappa,
      "ValidateFrobeniusPCSParamsOrThrow: kappa is too large for long");
  const long ell_prime = params.ell - params.kappa;
  if (params.ell_prime != ell_prime) {
    LogicError("ValidateFrobeniusPCSParamsOrThrow: ell_prime must equal ell-kappa");
  }

  if (params.basis_data.params.p != ZZ(2)) {
    LogicError("ValidateFrobeniusPCSParamsOrThrow: basis data must target p=2");
  }
  if (params.basis_data.params.k != base_log2) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: basis data k does not match base_modulus");
  }
  if (params.basis_data.params.r != expected_degree) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: basis data r does not match 2^kappa");
  }
  if (params.basis_data.modulus != params.base_modulus) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: basis data modulus does not match base_modulus");
  }
  if (params.basis_data.extension_modulus != params.extension_modulus) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: basis data extension modulus does not match params");
  }
  ::ValidateNormalBasisOrThrow(params.basis_data.normal_basis);

  const long basis_dimension =
      static_cast<long>(params.basis_data.normal_basis.beta.size());
  if (basis_dimension != expected_degree) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: beta basis dimension must equal 2^kappa");
  }
  if (static_cast<long>(params.basis_data.normal_basis.alpha.size()) !=
      expected_degree) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: alpha basis dimension must equal 2^kappa");
  }
  ValidatePrecomputedTablesOrThrow(params);

  const long expected_backend_message_length = Pow2LongOrThrow(
      ell_prime,
      "ValidateFrobeniusPCSParamsOrThrow: ell-kappa is too large for long");
  if (Z2kPCSBackendMessageLength(params.backend) !=
      expected_backend_message_length) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: backend message length must equal 2^(ell-kappa)");
  }
  if (Z2kPCSBackendPointDimension(params.backend) != ell_prime) {
    LogicError(
        "ValidateFrobeniusPCSParamsOrThrow: backend point dimension must equal ell-kappa");
  }
}

FrobeniusPCSParams FrobeniusPCSSetup(const FrobeniusPCSSetupInput &input) {
  if (input.kappa < 1) {
    LogicError("FrobeniusPCSSetup: kappa must be >= 1");
  }
  if (input.ell < input.kappa) {
    LogicError("FrobeniusPCSSetup: ell must be >= kappa");
  }

  const long base_log2 = Log2ExactPowerOfTwoZZOrThrow(
      input.base_modulus,
      "FrobeniusPCSSetup: base_modulus must be a power of two");
  const long expected_degree = Pow2LongOrThrow(
      input.kappa, "FrobeniusPCSSetup: kappa is too large for long");
  ValidateCurrentZ2kFrobeniusContextOrThrow(input.base_modulus,
                                            input.extension_modulus,
                                            input.kappa);
  Z2kPCSBackendValidateParamsOrThrow(input.backend);

  FrobeniusBasisParams basis_input;
  basis_input.p = ZZ(2);
  basis_input.k = base_log2;
  basis_input.r = expected_degree;
  basis_input.teichmuller_generator_max_trials =
      input.teichmuller_generator_max_trials;
  basis_input.affine_search_limit = input.affine_search_limit;

  FrobeniusPCSParams params;
  params.ell = input.ell;
  params.kappa = input.kappa;
  params.ell_prime = input.ell - input.kappa;
  params.base_modulus = input.base_modulus;
  params.extension_modulus = input.extension_modulus;
  if (input.use_provided_basis) {
    params.basis_data = ::BuildFrobeniusBasisFromProvidedNormalBasisOrThrow(
        basis_input, input.extension_modulus, input.provided_basis.normal_basis,
        input.provided_basis.has_teichmuller_generator,
        input.provided_basis.teichmuller_generator);
  } else {
    params.basis_data =
        ::BuildFrobeniusBasisOrThrow(basis_input, input.extension_modulus);
  }
  params.precomputed =
      BuildFrobeniusPCSPrecomputedTables(params.basis_data.normal_basis);
  params.backend = input.backend;
  ValidateFrobeniusPCSParamsOrThrow(params);
  return params;
}

NTL::vec_ZZ_pE PackZ2kTableToFrobeniusGREvals(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const long expected_length = Pow2LongOrThrow(
      params.ell, "PackZ2kTableToFrobeniusGREvals: ell is too large for long");
  if (t_table.length() != expected_length) {
    LogicError(
        "PackZ2kTableToFrobeniusGREvals: t_table length must equal 2^ell");
  }
  ValidateBaseRingVectorOrThrow(t_table, "t_table",
                                "PackZ2kTableToFrobeniusGREvals");

  return PackZ2kTableToFrobeniusGREvalsWithTrustedParamsOrThrow(params,
                                                                t_table);
}

NTL::vec_ZZ_pE PackZ2kTableToFrobeniusGREvalsWithTrustedParamsOrThrow(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table) {
  const long basis_dimension =
      static_cast<long>(params.basis_data.normal_basis.beta.size());
  const long packed_length = Pow2LongOrThrow(
      params.ell_prime,
      "PackZ2kTableToFrobeniusGREvalsWithTrustedParamsOrThrow: ell_prime is too large for long");

  const vector<ZZ_pE> &beta = params.basis_data.normal_basis.beta;

  NTL::vec_ZZ_pE packed;
  packed.SetLength(packed_length);
  for (long w = 0; w < packed_length; ++w) {
    ZZ_pE acc;
    NTL::clear(acc);
    for (long v = 0; v < basis_dimension; ++v) {
      acc += beta[static_cast<size_t>(v)] *
             t_table[v + w * basis_dimension];
    }
    packed[w] = acc;
  }
  return packed;
}

NTL::vec_ZZ_pE DecomposeGRElementToBaseCoeffsFrobeniusBasis(
    const FrobeniusPCSParams &params, const ZZ_pE &element) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const vector<ZZ_p> coords =
      ::RecoverNormalBasisCoords(params.basis_data.normal_basis, element);
  NTL::vec_ZZ_pE out;
  out.SetLength(static_cast<long>(coords.size()));
  for (long i = 0; i < static_cast<long>(coords.size()); ++i) {
    out[i] = BaseRingConstant(coords[static_cast<size_t>(i)]);
  }
  return out;
}

MerkleRoot FrobeniusPCSCommit(const FrobeniusPCSParams &params,
                              const NTL::vec_ZZ_pE &t_table) {
  return FrobeniusPCSBuildCommitArtifacts(params, t_table).commitment;
}

FrobeniusPCSOuterCommitArtifacts FrobeniusPCSBuildOuterCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const long expected_length = Pow2LongOrThrow(
      params.ell,
      "FrobeniusPCSBuildOuterCommitArtifacts: ell is too large for long");
  if (t_table.length() != expected_length) {
    LogicError(
        "FrobeniusPCSBuildOuterCommitArtifacts: t_table length must equal 2^ell");
  }
  ValidateBaseRingVectorOrThrow(
      t_table, "t_table", "FrobeniusPCSBuildOuterCommitArtifacts");
  return FrobeniusPCSBuildOuterCommitArtifactsUnchecked(params, t_table);
}

FrobeniusPCSOuterCommitArtifacts FrobeniusPCSBuildOuterCommitArtifactsUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table) {
  const PackedCommitInputs packed = BuildPackedCommitInputs(params, t_table);

  FrobeniusPCSOuterCommitArtifacts out;
  out.t_packed_table = packed.t_packed_table;
  out.t_packed_monomial_coeffs = packed.t_packed_monomial_coeffs;
  return out;
}

FrobeniusPCSCommitArtifacts FrobeniusPCSBuildCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  const long expected_length = Pow2LongOrThrow(
      params.ell,
      "FrobeniusPCSBuildCommitArtifacts: ell is too large for long");
  if (t_table.length() != expected_length) {
    LogicError("FrobeniusPCSBuildCommitArtifacts: t_table length must equal 2^ell");
  }
  ValidateBaseRingVectorOrThrow(
      t_table, "t_table", "FrobeniusPCSBuildCommitArtifacts");
  return FrobeniusPCSBuildCommitArtifactsUnchecked(params, t_table);
}

FrobeniusPCSCommitArtifacts FrobeniusPCSBuildCommitArtifactsUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table) {
  const FrobeniusPCSOuterCommitArtifacts outer =
      FrobeniusPCSBuildOuterCommitArtifactsUnchecked(params, t_table);

  FrobeniusPCSCommitArtifacts out;
  out.t_packed_table = outer.t_packed_table;
  out.t_packed_monomial_coeffs = outer.t_packed_monomial_coeffs;
  out.backend_commit_artifacts = Z2kPCSBackendBuildCommitArtifacts(
      params.backend, out.t_packed_monomial_coeffs);
  out.commitment = out.backend_commit_artifacts.commitment;
  return out;
}

FrobeniusPCSOuterEvalProof FrobeniusPCSProveOuterEval(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries) {
  const FrobeniusPCSOuterCommitArtifacts commit_artifacts =
      FrobeniusPCSBuildOuterCommitArtifacts(params, t_table);
  return FrobeniusPCSProveOuterEvalFromCommitArtifacts(
      params, t_table, commitment, z, claimed_s, num_queries, commit_artifacts);
}

FrobeniusPCSOuterEvalProof FrobeniusPCSProveOuterEvalUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries) {
  const FrobeniusPCSOuterCommitArtifacts commit_artifacts =
      FrobeniusPCSBuildOuterCommitArtifacts(params, t_table);
  return FrobeniusPCSProveOuterEvalFromCommitArtifactsUnchecked(
      params, t_table, commitment, z, claimed_s, num_queries, commit_artifacts);
}

FrobeniusPCSOuterEvalProof FrobeniusPCSProveOuterEvalFromCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries,
    const FrobeniusPCSOuterCommitArtifacts &commit_artifacts) {
  return ProveOuterEvalFromCommitArtifactsInternal(
             params, t_table, commitment, z, claimed_s, num_queries,
             commit_artifacts, /*checked_path=*/true,
             "FrobeniusPCSProveOuterEvalFromCommitArtifacts")
      .proof;
}

FrobeniusPCSOuterEvalProof FrobeniusPCSProveOuterEvalFromCommitArtifactsUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const MerkleRoot &commitment, const std::vector<FieldElement> &z,
    const FieldElement &claimed_s, long num_queries,
    const FrobeniusPCSOuterCommitArtifacts &commit_artifacts) {
  return ProveOuterEvalFromCommitArtifactsInternal(
             params, t_table, commitment, z, claimed_s, num_queries,
             commit_artifacts, /*checked_path=*/false,
             "FrobeniusPCSProveOuterEvalFromCommitArtifactsUnchecked")
      .proof;
}

FrobeniusPCSEvalProof FrobeniusPCSProveEval(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries) {
  const FrobeniusPCSCommitArtifacts commit_artifacts =
      FrobeniusPCSBuildCommitArtifacts(params, t_table);
  return FrobeniusPCSProveEvalFromCommitArtifacts(
      params, t_table, z, claimed_s, num_queries, commit_artifacts);
}

FrobeniusPCSEvalProof FrobeniusPCSProveEvalUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries) {
  const FrobeniusPCSCommitArtifacts commit_artifacts =
      FrobeniusPCSBuildCommitArtifacts(params, t_table);
  return FrobeniusPCSProveEvalFromCommitArtifactsUnchecked(
      params, t_table, z, claimed_s, num_queries, commit_artifacts);
}

FrobeniusPCSEvalProof FrobeniusPCSProveEvalFromCommitArtifacts(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const FrobeniusPCSCommitArtifacts &commit_artifacts) {
  return ProveEvalFromCommitArtifactsInternal(
      params, t_table, z, claimed_s, num_queries, commit_artifacts,
      /*checked_path=*/true, "FrobeniusPCSProveEvalFromCommitArtifacts");
}

FrobeniusPCSEvalProof FrobeniusPCSProveEvalFromCommitArtifactsUnchecked(
    const FrobeniusPCSParams &params, const NTL::vec_ZZ_pE &t_table,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const FrobeniusPCSCommitArtifacts &commit_artifacts) {
  return ProveEvalFromCommitArtifactsInternal(
      params, t_table, z, claimed_s, num_queries, commit_artifacts,
      /*checked_path=*/false,
      "FrobeniusPCSProveEvalFromCommitArtifactsUnchecked");
}

bool FrobeniusPCSVerifyEval(const FrobeniusPCSParams &params,
                            const MerkleRoot &commitment,
                            const std::vector<FieldElement> &z,
                            const FieldElement &claimed_s, long num_queries,
                            const FrobeniusPCSEvalProof &proof) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  if (!HasExpectedEvalProofShape(params, proof) ||
      !HasCompatibleBackendEvalSubproof(params, proof.backend_proof)) {
    return false;
  }

  FrobeniusPCSOuterEvalProof outer_proof;
  outer_proof.s_by_i = proof.s_by_i;
  outer_proof.h_by_level = proof.h_by_level;
  outer_proof.t_star = proof.t_star;

  std::vector<FieldElement> rprime_suffix;
  if (!VerifyOuterEvalAndMaybeRecoverSuffix(params, commitment, z, claimed_s,
                                            num_queries, outer_proof,
                                            &rprime_suffix)) {
    return false;
  }

  {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->z2k_backend_verify_ns : nullptr,
                      prof ? &prof->z2k_backend_verify_calls : nullptr);
    return Z2kPCSBackendVerifyEval(params.backend, commitment, rprime_suffix,
                                   proof.t_star, num_queries,
                                   proof.backend_proof);
  }
}

bool FrobeniusPCSVerifyEvalUnchecked(const FrobeniusPCSParams &params,
                                     const MerkleRoot &commitment,
                                     const std::vector<FieldElement> &z,
                                     const FieldElement &claimed_s,
                                     long num_queries,
                                     const FrobeniusPCSEvalProof &proof) {
  if (!HasExpectedEvalProofShape(params, proof) ||
      !HasCompatibleBackendEvalSubproof(params, proof.backend_proof)) {
    return false;
  }

  FrobeniusPCSOuterEvalProof outer_proof;
  outer_proof.s_by_i = proof.s_by_i;
  outer_proof.h_by_level = proof.h_by_level;
  outer_proof.t_star = proof.t_star;

  std::vector<FieldElement> rprime_suffix;
  if (!VerifyOuterEvalAndMaybeRecoverSuffix(params, commitment, z, claimed_s,
                                            num_queries, outer_proof,
                                            &rprime_suffix)) {
    return false;
  }

  {
    Profile *prof = ActiveProfile();
    ScopedTimer timer(prof ? &prof->z2k_backend_verify_ns : nullptr,
                      prof ? &prof->z2k_backend_verify_calls : nullptr);
    return Z2kPCSBackendVerifyEvalUnchecked(
        params.backend, commitment, rprime_suffix, proof.t_star, num_queries,
        proof.backend_proof);
  }
}

bool FrobeniusPCSVerifyOuterEval(const FrobeniusPCSParams &params,
                                 const MerkleRoot &commitment,
                                 const std::vector<FieldElement> &z,
                                 const FieldElement &claimed_s,
                                 long num_queries,
                                 const FrobeniusPCSOuterEvalProof &proof) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  return FrobeniusPCSVerifyOuterEvalUnchecked(params, commitment, z, claimed_s,
                                              num_queries, proof);
}

bool FrobeniusPCSVerifyOuterEvalUnchecked(
    const FrobeniusPCSParams &params, const MerkleRoot &commitment,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const FrobeniusPCSOuterEvalProof &proof) {
  return VerifyOuterEvalAndMaybeRecoverSuffix(params, commitment, z, claimed_s,
                                              num_queries, proof, nullptr);
}

bool FrobeniusPCSRecoverBackendEvaluationPoint(
    const FrobeniusPCSParams &params, const MerkleRoot &commitment,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const FrobeniusPCSOuterEvalProof &proof,
    std::vector<FieldElement> &backend_eval_point_out) {
  ValidateFrobeniusPCSParamsOrThrow(params);
  return FrobeniusPCSRecoverBackendEvaluationPointUnchecked(
      params, commitment, z, claimed_s, num_queries, proof,
      backend_eval_point_out);
}

bool FrobeniusPCSRecoverBackendEvaluationPointUnchecked(
    const FrobeniusPCSParams &params, const MerkleRoot &commitment,
    const std::vector<FieldElement> &z, const FieldElement &claimed_s,
    long num_queries, const FrobeniusPCSOuterEvalProof &proof,
    std::vector<FieldElement> &backend_eval_point_out) {
  return VerifyOuterEvalAndMaybeRecoverSuffix(params, commitment, z, claimed_s,
                                              num_queries, proof,
                                              &backend_eval_point_out);
}

}  // namespace basefold
