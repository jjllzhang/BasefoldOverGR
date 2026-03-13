#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/vec_ZZ_pE.h>

#include <cstdlib>
#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "GaloisRing/Basis.hpp"
#include "GaloisRing/FrobeniusBasis.hpp"
#include "GaloisRing/utils.hpp"
#include "PCS/Common/Multilinear.hpp"
#include "tests/test_common.hpp"

using NTL::SetCoeff;
using NTL::conv;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using std::cerr;
using std::cout;
using std::exception;
using std::string;
using std::vector;

int g_test_failure_count = 0;

namespace {

struct GRSpec {
  const char *name = nullptr;
  ZZ p;
  long k = 0;
  long r = 0;
  ZZ modulus;
  vector<long> extension_coeffs;
};

GRSpec MakeGR42Degree2Spec() {
  GRSpec spec;
  spec.name = "GR(4,2)";
  spec.p = to_ZZ(2);
  spec.k = 2;
  spec.r = 2;
  spec.modulus = to_ZZ(4);
  spec.extension_coeffs = {1, 1, 1};
  return spec;
}

GRSpec MakeGR42Degree4Spec() {
  GRSpec spec;
  spec.name = "GR(4,4)";
  spec.p = to_ZZ(2);
  spec.k = 2;
  spec.r = 4;
  spec.modulus = to_ZZ(4);
  spec.extension_coeffs = {1, 1, 0, 0, 1};
  return spec;
}

ZZ_pX BuildExtensionPolynomial(const GRSpec &spec) {
  ZZ_pX out;
  for (long i = 0; i < static_cast<long>(spec.extension_coeffs.size()); ++i) {
    SetCoeff(out, i, spec.extension_coeffs[static_cast<std::size_t>(i)]);
  }
  return out;
}

FrobeniusBasisParams MakeBasisParams(const GRSpec &spec) {
  FrobeniusBasisParams params;
  params.p = spec.p;
  params.k = spec.k;
  params.r = spec.r;
  params.teichmuller_generator_max_trials = 1024;
  params.affine_search_limit = 4;
  return params;
}

long PowInt(long base, long exp) {
  CHECK_MSG(base >= 0, "PowInt: base must be non-negative");
  CHECK_MSG(exp >= 0, "PowInt: exponent must be non-negative");
  long out = 1;
  for (long i = 0; i < exp; ++i) {
    out *= base;
  }
  return out;
}

long CurrentBaseModulusLong() {
  return NTL::conv<long>(ZZ_p::modulus());
}

ZZ_pE ConstFromZZp(const ZZ_p &value) {
  ZZ_pX poly;
  SetCoeff(poly, 0, value);
  ZZ_pE out;
  conv(out, poly);
  return out;
}

ZZ_pE ConstFromLong(long value) {
  return testutil::ConstZZpE(value);
}

vector<long> CoordsToLongs(const vector<ZZ_p> &coords) {
  vector<long> out(coords.size(), 0);
  for (long i = 0; i < static_cast<long>(coords.size()); ++i) {
    out[static_cast<std::size_t>(i)] =
        NTL::conv<long>(NTL::rep(coords[static_cast<std::size_t>(i)]));
  }
  return out;
}

vector<long> RotateRight(const vector<long> &v) {
  if (v.empty()) {
    return v;
  }
  vector<long> out = v;
  out[0] = v.back();
  for (long i = 1; i < static_cast<long>(v.size()); ++i) {
    out[static_cast<std::size_t>(i)] = v[static_cast<std::size_t>(i - 1)];
  }
  return out;
}

ZZ_pE ApplyTauPower(const NormalBasisData &normal_basis, const ZZ_pE &element,
                    long exp) {
  ZZ_pE out = element;
  for (long i = 0; i < exp; ++i) {
    out = ApplyFrobeniusTau(normal_basis, out);
  }
  return out;
}

vector<ZZ_pE> ApplySigmaPowerToPoint(const NormalBasisData &normal_basis,
                                     const vector<ZZ_pE> &point, long exp) {
  vector<ZZ_pE> out(point.size());
  for (long i = 0; i < static_cast<long>(point.size()); ++i) {
    ZZ_pE value = point[static_cast<std::size_t>(i)];
    for (long round = 0; round < exp; ++round) {
      value = ApplyFrobeniusSigma(normal_basis, value);
    }
    out[static_cast<std::size_t>(i)] = value;
  }
  return out;
}

vector<ZZ_pE> EnumerateAllExtensionElements(long degree) {
  const long modulus = CurrentBaseModulusLong();
  const long total = PowInt(modulus, degree);
  vector<ZZ_pE> elements;
  elements.reserve(static_cast<std::size_t>(total));

  for (long idx = 0; idx < total; ++idx) {
    long tmp = idx;
    vector<long> coeffs(static_cast<std::size_t>(degree), 0);
    for (long i = 0; i < degree; ++i) {
      coeffs[static_cast<std::size_t>(i)] = tmp % modulus;
      tmp /= modulus;
    }
    elements.push_back(LongVecToZZpE(coeffs));
  }
  return elements;
}

vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] = ConstFromLong((index >> i) & 1L);
  }
  return point;
}

vec_ZZ_pE BuildBaseRingTable(const vector<long> &values) {
  vec_ZZ_pE table;
  table.SetLength(static_cast<long>(values.size()));
  for (long i = 0; i < static_cast<long>(values.size()); ++i) {
    table[i] = ConstFromLong(values[static_cast<std::size_t>(i)]);
  }
  return table;
}

ZZ_pE EvalFromBooleanTable(const vec_ZZ_pE &table, long dimension,
                           const vector<ZZ_pE> &point) {
  CHECK_EQ(table.length(), PowInt(2, dimension));
  ZZ_pE acc = ZZ_pE(0);
  for (long idx = 0; idx < table.length(); ++idx) {
    acc += table[idx] *
           basefold::EqPolynomial(point, BooleanPointFromIndex(idx, dimension));
  }
  return acc;
}

vector<ZZ_pE> DirectPartialEvaluations(const vec_ZZ_pE &t_table, long kappa,
                                       long ell_prime,
                                       const vector<ZZ_pE> &r_suffix) {
  const long basis_dimension = PowInt(2, kappa);
  const long num_w = PowInt(2, ell_prime);
  vector<ZZ_pE> partials(static_cast<std::size_t>(basis_dimension), ZZ_pE(0));
  for (long u = 0; u < basis_dimension; ++u) {
    ZZ_pE acc = ZZ_pE(0);
    for (long w = 0; w < num_w; ++w) {
      acc += t_table[u + (w << kappa)] *
             basefold::EqPolynomial(r_suffix, BooleanPointFromIndex(w, ell_prime));
    }
    partials[static_cast<std::size_t>(u)] = acc;
  }
  return partials;
}

vec_ZZ_pE PackByBasis(const vec_ZZ_pE &t_table, const vector<ZZ_pE> &basis,
                      long kappa, long ell_prime) {
  const long basis_dimension = PowInt(2, kappa);
  const long num_w = PowInt(2, ell_prime);
  CHECK_EQ(static_cast<long>(basis.size()), basis_dimension);

  vec_ZZ_pE packed;
  packed.SetLength(num_w);
  for (long w = 0; w < num_w; ++w) {
    ZZ_pE acc = ZZ_pE(0);
    for (long v = 0; v < basis_dimension; ++v) {
      acc += basis[static_cast<std::size_t>(v)] * t_table[v + (w << kappa)];
    }
    packed[w] = acc;
  }
  return packed;
}

vector<ZZ_pE> RecoverPartialsFromFrobeniusOrbit(
    const NormalBasisData &normal_basis, const vec_ZZ_pE &t_prime_table,
    long ell_prime, const vector<ZZ_pE> &r_suffix) {
  const long basis_dimension = static_cast<long>(normal_basis.alpha.size());
  vector<ZZ_pE> s_by_i(static_cast<std::size_t>(basis_dimension), ZZ_pE(0));
  for (long i = 0; i < basis_dimension; ++i) {
    const vector<ZZ_pE> sigma_point =
        ApplySigmaPowerToPoint(normal_basis, r_suffix, i);
    s_by_i[static_cast<std::size_t>(i)] =
        EvalFromBooleanTable(t_prime_table, ell_prime, sigma_point);
  }

  vector<ZZ_pE> partials(static_cast<std::size_t>(basis_dimension), ZZ_pE(0));
  for (long u = 0; u < basis_dimension; ++u) {
    ZZ_pE acc = ZZ_pE(0);
    for (long i = 0; i < basis_dimension; ++i) {
      acc += ApplyTauPower(normal_basis,
                           normal_basis.alpha[static_cast<std::size_t>(u)] *
                               s_by_i[static_cast<std::size_t>(i)],
                           i);
    }
    partials[static_cast<std::size_t>(u)] = acc;
  }
  return partials;
}

ZZ_pE RecombineClaim(const vector<ZZ_pE> &partials,
                     const vector<ZZ_pE> &z_prefix) {
  ZZ_pE acc = ZZ_pE(0);
  for (long idx = 0; idx < static_cast<long>(partials.size()); ++idx) {
    acc += partials[static_cast<std::size_t>(idx)] *
           basefold::EqPolynomial(
               z_prefix,
               BooleanPointFromIndex(idx, static_cast<long>(z_prefix.size())));
  }
  return acc;
}

void ValidateFrobeniusSemantics(const GRSpec &spec) {
  testutil::PrintInfo(string("Frobenius semantics over ") + spec.name);
  ZZ_pPush mod_push(spec.modulus);
  const ZZ_pX extension_modulus = BuildExtensionPolynomial(spec);
  ZZ_pEPush ext_push(extension_modulus);

  const FrobeniusBasisData basis_data =
      BuildFrobeniusBasisOrThrow(MakeBasisParams(spec), extension_modulus);
  const NormalBasisData &normal_basis = basis_data.normal_basis;
  CHECK_MSG(!IsBaseRingConstant(basis_data.teichmuller_generator),
            "ValidateFrobeniusSemantics: Teichmuller generator must be non-trivial");
  ValidateNormalBasisOrThrow(normal_basis);

  const vector<ZZ_pE> all_elements = EnumerateAllExtensionElements(spec.r);
  for (const ZZ_pE &element : all_elements) {
    CHECK_EQ(ApplyTauPower(normal_basis, element, spec.r), element);

    ZZ_pE orbit_trace = ZZ_pE(0);
    for (long i = 0; i < spec.r; ++i) {
      orbit_trace += ApplyTauPower(normal_basis, element, i);
    }
    CHECK_EQ(orbit_trace, ConstFromZZp(TraceToBaseRing(element)));
    CHECK(IsBaseRingConstant(orbit_trace));

    const vector<ZZ_p> coords = RecoverNormalBasisCoords(normal_basis, element);
    CHECK_EQ(ComposeFromNormalBasisCoords(normal_basis, coords), element);

    const vector<ZZ_p> tau_coords = RecoverNormalBasisCoords(
        normal_basis, ApplyFrobeniusTau(normal_basis, element));
    CHECK_EQ(CoordsToLongs(tau_coords), RotateRight(CoordsToLongs(coords)));
  }

  for (long c = 0; c < CurrentBaseModulusLong(); ++c) {
    CHECK_EQ(ApplyFrobeniusTau(normal_basis, ConstFromLong(c)),
             ConstFromLong(c));
  }

  for (long i = 0; i < spec.r; ++i) {
    CHECK_EQ(normal_basis.beta[static_cast<std::size_t>((i + 1) % spec.r)],
             ApplyFrobeniusTau(
                 normal_basis,
                 normal_basis.beta[static_cast<std::size_t>(i)]));
    CHECK_EQ(normal_basis.alpha[static_cast<std::size_t>((i + 1) % spec.r)],
             ApplyFrobeniusTau(
                 normal_basis,
                 normal_basis.alpha[static_cast<std::size_t>(i)]));
  }
}

void ValidateProtocol2PartialRecovery(const GRSpec &spec, long kappa,
                                      const vector<long> &table_values,
                                      const vector<ZZ_pE> &z_prefix,
                                      const vector<ZZ_pE> &r_suffix) {
  testutil::PrintInfo(string("Protocol 2 partial recovery over ") + spec.name);
  ZZ_pPush mod_push(spec.modulus);
  const ZZ_pX extension_modulus = BuildExtensionPolynomial(spec);
  ZZ_pEPush ext_push(extension_modulus);

  const FrobeniusBasisData basis_data =
      BuildFrobeniusBasisOrThrow(MakeBasisParams(spec), extension_modulus);
  const NormalBasisData &normal_basis = basis_data.normal_basis;

  const long basis_dimension = PowInt(2, kappa);
  CHECK_EQ(basis_dimension, spec.r);
  const long ell_prime = static_cast<long>(r_suffix.size());
  const long ell = kappa + ell_prime;
  CHECK_EQ(static_cast<long>(table_values.size()), PowInt(2, ell));
  CHECK_EQ(static_cast<long>(z_prefix.size()), kappa);

  const vec_ZZ_pE t_table = BuildBaseRingTable(table_values);
  const vec_ZZ_pE t_prime_table =
      PackByBasis(t_table, normal_basis.beta, kappa, ell_prime);

  const vector<ZZ_pE> direct_partials =
      DirectPartialEvaluations(t_table, kappa, ell_prime, r_suffix);
  const vector<ZZ_pE> recovered_partials = RecoverPartialsFromFrobeniusOrbit(
      normal_basis, t_prime_table, ell_prime, r_suffix);
  CHECK_EQ(recovered_partials, direct_partials);

  vector<ZZ_pE> z = z_prefix;
  z.insert(z.end(), r_suffix.begin(), r_suffix.end());
  const ZZ_pE direct_eval = EvalFromBooleanTable(t_table, ell, z);
  const ZZ_pE recombined = RecombineClaim(recovered_partials, z_prefix);
  CHECK_EQ(recombined, direct_eval);
}

void TestFrobeniusSemantics_GR42Degree2() {
  ValidateFrobeniusSemantics(MakeGR42Degree2Spec());
}

void TestFrobeniusSemantics_GR42Degree4() {
  ValidateFrobeniusSemantics(MakeGR42Degree4Spec());
}

void TestProtocol2PartialRecovery_GR42Degree2() {
  const GRSpec spec = MakeGR42Degree2Spec();
  ZZ_pPush mod_push(spec.modulus);
  const ZZ_pX extension_modulus = BuildExtensionPolynomial(spec);
  ZZ_pEPush ext_push(extension_modulus);
  const FrobeniusBasisData basis_data =
      BuildFrobeniusBasisOrThrow(MakeBasisParams(spec), extension_modulus);

  const vector<ZZ_pE> z_prefix = {
      basis_data.teichmuller_generator + ConstFromLong(1)};
  const vector<ZZ_pE> r_suffix = {basis_data.teichmuller_generator,
                                  ConstFromLong(3)};
  const vector<long> table_values = {0, 1, 2, 3, 1, 0, 3, 2};
  ValidateProtocol2PartialRecovery(spec, /*kappa=*/1, table_values, z_prefix,
                                   r_suffix);
}

void TestProtocol2PartialRecovery_GR42Degree4() {
  const GRSpec spec = MakeGR42Degree4Spec();
  ZZ_pPush mod_push(spec.modulus);
  const ZZ_pX extension_modulus = BuildExtensionPolynomial(spec);
  ZZ_pEPush ext_push(extension_modulus);
  const FrobeniusBasisData basis_data =
      BuildFrobeniusBasisOrThrow(MakeBasisParams(spec), extension_modulus);

  const vector<ZZ_pE> z_prefix = {basis_data.teichmuller_generator +
                                      ConstFromLong(1),
                                  power(basis_data.teichmuller_generator, 2) +
                                      ConstFromLong(2)};
  const vector<ZZ_pE> r_suffix = {power(basis_data.teichmuller_generator, 3) +
                                      ConstFromLong(1),
                                  basis_data.teichmuller_generator +
                                      ConstFromLong(3)};
  const vector<long> table_values = {
      0, 1, 2, 3, 1, 2, 3, 0, 2, 0, 1, 3, 3, 1, 0, 2};
  ValidateProtocol2PartialRecovery(spec, /*kappa=*/2, table_values, z_prefix,
                                   r_suffix);
}

}  // namespace

int main() {
  try {
    RUN_TEST(TestFrobeniusSemantics_GR42Degree2);
    RUN_TEST(TestFrobeniusSemantics_GR42Degree4);
    RUN_TEST(TestProtocol2PartialRecovery_GR42Degree2);
    RUN_TEST(TestProtocol2PartialRecovery_GR42Degree4);
  } catch (const exception &e) {
    cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    cerr << "Unhandled non-std exception\n";
    return 2;
  }

  if (g_test_failure_count == 0) {
    cout << "\nAll tests passed.\n";
    return 0;
  }

  cerr << "\n" << g_test_failure_count << " test(s) failed.\n";
  return 1;
}
