#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <algorithm>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "PCS/BaseFold/BaseFoldPCS.hpp"
#include "PCS/Common/Multilinear.hpp"
#include "PCS/Common/Sumcheck.hpp"
#include "PCS/Common/Transcript.hpp"
#include "Compiler/Z2k/BaseFoldBackendAdapter.hpp"
#include "Compiler/Z2k/FrobeniusPCS.hpp"
#include "Compiler/Z2k/FrobeniusProofSerialize.hpp"
#include "tests/test_common.hpp"

using NTL::conv;
using NTL::mat_ZZ_pE;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using std::cerr;
using std::cout;
using std::exception;
using std::function;
using std::make_shared;
using std::string;
using std::vector;

int g_test_failure_count = 0;

namespace {

basefold::FoldableCodeParams BuildParamsGR42(const ZZ &p, const ZZ_pE &alpha) {
  const long c = 2;
  const long k0 = 1;
  const long d = 2;
  const long n0 = c * k0;

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = p;
  params.zeta = alpha;

  mat_ZZ_pE G0;
  G0.SetDims(k0, n0);
  const ZZ_pE one = testutil::ConstZZpE(1);
  G0[0][0] = one;
  G0[0][1] = alpha;
  params.G0 = G0;

  params.diag_T.resize(static_cast<std::size_t>(d));
  params.diag_T[0].SetLength(2);
  params.diag_T[0][0] = one;
  params.diag_T[0][1] = alpha;

  params.diag_T[1].SetLength(4);
  params.diag_T[1][0] = alpha + one;
  params.diag_T[1][1] = one;
  params.diag_T[1][2] = alpha;
  params.diag_T[1][3] = alpha + one;
  return params;
}

basefold::FoldableCodeParams BuildParamsGR42D0(const ZZ &p,
                                               const ZZ_pE &alpha) {
  const long c = 2;
  const long k0 = 1;
  const long d = 0;
  const long n0 = c * k0;

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = p;
  params.zeta = alpha;

  mat_ZZ_pE G0;
  G0.SetDims(k0, n0);
  const ZZ_pE one = testutil::ConstZZpE(1);
  G0[0][0] = one;
  G0[0][1] = alpha;
  params.G0 = G0;
  return params;
}

vec_ZZ_pE BuildBaseRingCoeffVector(const vector<long> &coeffs) {
  vec_ZZ_pE out;
  out.SetLength(static_cast<long>(coeffs.size()));
  for (long i = 0; i < static_cast<long>(coeffs.size()); ++i) {
    out[i] = testutil::ConstZZpE(coeffs[static_cast<std::size_t>(i)]);
  }
  return out;
}

long Pow2ForTest(long exponent) {
  CHECK_MSG(exponent >= 0, "Pow2ForTest: exponent must be non-negative");
  if (exponent < 0) {
    return 0;
  }
  return 1L << exponent;
}

vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] =
        testutil::ConstZZpE((index >> i) & 1L);
  }
  return point;
}

ZZ_pE EvalFromBooleanTable(const vec_ZZ_pE &table, long dimension,
                           const vector<ZZ_pE> &point) {
  CHECK_EQ(table.length(), Pow2ForTest(dimension));
  ZZ_pE acc = ZZ_pE(0);
  for (long idx = 0; idx < table.length(); ++idx) {
    acc += table[idx] *
           basefold::EqPolynomial(point, BooleanPointFromIndex(idx, dimension));
  }
  return acc;
}

std::uint64_t FixedFieldElementBytesForCurrentContext() {
  const ZZ modulus_minus_one = NTL::ZZ_p::modulus() - 1;
  const long bits = NTL::NumBits(modulus_minus_one);
  CHECK_MSG(bits > 0,
            "FixedFieldElementBytesForCurrentContext: invalid modulus bit width");
  const std::uint64_t coeff_bytes = static_cast<std::uint64_t>((bits + 7) / 8);
  return coeff_bytes * static_cast<std::uint64_t>(ZZ_pE::degree());
}

ZZ_pE ApplyTauPower(const basefold::FrobeniusPCSParams &params,
                    const ZZ_pE &element, long exp) {
  ZZ_pE out = element;
  for (long i = 0; i < exp; ++i) {
    out = ApplyFrobeniusTau(params.basis_data.normal_basis, out);
  }
  return out;
}

vector<ZZ_pE> ApplySigmaPowerToPoint(const basefold::FrobeniusPCSParams &params,
                                     const vector<ZZ_pE> &point, long exp) {
  vector<ZZ_pE> out = point;
  for (ZZ_pE &value : out) {
    for (long i = 0; i < exp; ++i) {
      value = ApplyFrobeniusSigma(params.basis_data.normal_basis, value);
    }
  }
  return out;
}

vector<ZZ_pE> RecoverPartialsFromProof(
    const basefold::FrobeniusPCSParams &params,
    const vector<ZZ_pE> &s_by_i) {
  const vector<ZZ_pE> &alpha = params.basis_data.normal_basis.alpha;
  CHECK_EQ(static_cast<long>(s_by_i.size()), static_cast<long>(alpha.size()));

  vector<ZZ_pE> partials(alpha.size(), ZZ_pE(0));
  for (long u = 0; u < static_cast<long>(alpha.size()); ++u) {
    ZZ_pE acc = ZZ_pE(0);
    for (long i = 0; i < static_cast<long>(alpha.size()); ++i) {
      acc += ApplyTauPower(params, alpha[static_cast<std::size_t>(u)] *
                                       s_by_i[static_cast<std::size_t>(i)],
                           i);
    }
    partials[static_cast<std::size_t>(u)] = acc;
  }
  return partials;
}

vector<ZZ_pE> DirectPartialEvaluations(const vec_ZZ_pE &t_table, long kappa,
                                       long ell_prime,
                                       const vector<ZZ_pE> &r_suffix) {
  const long basis_dimension = Pow2ForTest(kappa);
  const long num_w = Pow2ForTest(ell_prime);
  vector<ZZ_pE> partials(static_cast<std::size_t>(basis_dimension), ZZ_pE(0));
  for (long u = 0; u < basis_dimension; ++u) {
    vec_ZZ_pE slice;
    slice.SetLength(num_w);
    for (long w = 0; w < num_w; ++w) {
      slice[w] = t_table[u + (w << kappa)];
    }
    partials[static_cast<std::size_t>(u)] =
        EvalFromBooleanTable(slice, ell_prime, r_suffix);
  }
  return partials;
}

basefold::HashTranscript MakeFrobeniusTranscriptForTest() {
  basefold::HashTranscriptConfig config;
  config.domain_separator = "FrobeniusPCS/v1";
  config.byte_order = basefold::TranscriptByteOrder::kLittleEndian;
  config.error_prefix = "FrobeniusHashTranscript";
  return basefold::HashTranscript(config);
}

struct FrobeniusChallengeTrace {
  vector<ZZ_pE> rprime_prefix;
  vector<ZZ_pE> rprime_suffix;
};

FrobeniusChallengeTrace ReplayFrobeniusChallenges(
    const basefold::MerkleRoot &commitment, const vector<ZZ_pE> &z,
    const ZZ_pE &claimed_s, const basefold::FrobeniusPCSOuterEvalProof &proof,
    long kappa, long ell_prime) {
  basefold::HashTranscript transcript = MakeFrobeniusTranscriptForTest();
  transcript.AbsorbDigest(commitment);
  for (const ZZ_pE &zi : z) {
    transcript.AbsorbFieldElement(zi);
  }
  transcript.AbsorbFieldElement(claimed_s);
  for (const ZZ_pE &s_i : proof.s_by_i) {
    transcript.AbsorbFieldElement(s_i);
  }

  FrobeniusChallengeTrace trace;
  trace.rprime_prefix.resize(static_cast<std::size_t>(kappa));
  for (long i = 0; i < kappa; ++i) {
    trace.rprime_prefix[static_cast<std::size_t>(i)] =
        transcript.ChallengeFieldElement("rprime/prefix/" + std::to_string(i));
  }

  trace.rprime_suffix.resize(static_cast<std::size_t>(ell_prime));
  if (ell_prime == 0) {
    return trace;
  }

  AbsorbQuadraticPoly(
      transcript, proof.h_by_level[static_cast<std::size_t>(ell_prime - 1)]);
  for (long i = ell_prime; i-- > 0;) {
    trace.rprime_suffix[static_cast<std::size_t>(i)] =
        transcript.ChallengeFieldElement("rprime/suffix/" + std::to_string(i));
    if (i > 0) {
      AbsorbQuadraticPoly(
          transcript, proof.h_by_level[static_cast<std::size_t>(i - 1)]);
    }
  }
  return trace;
}

FrobeniusChallengeTrace ReplayFrobeniusChallenges(
    const basefold::MerkleRoot &commitment, const vector<ZZ_pE> &z,
    const ZZ_pE &claimed_s, const basefold::FrobeniusPCSEvalProof &proof,
    long kappa, long ell_prime) {
  basefold::FrobeniusPCSOuterEvalProof outer;
  outer.s_by_i = proof.s_by_i;
  outer.h_by_level = proof.h_by_level;
  outer.t_star = proof.t_star;
  return ReplayFrobeniusChallenges(commitment, z, claimed_s, outer, kappa,
                                   ell_prime);
}

basefold::Z2kPCSBackendEvalProof MutateBaseFoldBackendSubproof(
    const basefold::Z2kPCSBackendEvalProof &opaque_proof) {
  const auto basefold_ptr =
      std::static_pointer_cast<const basefold::BaseFoldPCSEvalProof>(
          opaque_proof.payload);
  CHECK(static_cast<bool>(basefold_ptr));

  basefold::BaseFoldPCSEvalProof mutated = *basefold_ptr;
  if (!mutated.h_by_level.empty()) {
    mutated.h_by_level[0].a0 += testutil::ConstZZpE(1);
  } else if (mutated.pi0_codeword.length() > 0) {
    mutated.pi0_codeword[0] += testutil::ConstZZpE(1);
  } else {
    CHECK_MSG(false, "MutateBaseFoldBackendSubproof: no mutable payload found");
  }

  basefold::Z2kPCSBackendEvalProof out = opaque_proof;
  out.payload =
      std::static_pointer_cast<const void>(make_shared<basefold::BaseFoldPCSEvalProof>(
          std::move(mutated)));
  return out;
}

void ExpectChildFailureContains(const function<void()> &fn, const string &needle,
                                const string &label) {
#if defined(__unix__) || defined(__APPLE__)
  int pipe_fds[2];
  CHECK_MSG(pipe(pipe_fds) == 0, label + ": pipe() failed");
  if (g_test_failure_count != 0) {
    return;
  }

  cout.flush();
  cerr.flush();
  const pid_t pid = fork();
  CHECK_MSG(pid >= 0, label + ": fork() failed");
  if (g_test_failure_count != 0) {
    close(pipe_fds[0]);
    close(pipe_fds[1]);
    return;
  }

  if (pid == 0) {
    close(pipe_fds[0]);
    dup2(pipe_fds[1], STDOUT_FILENO);
    dup2(pipe_fds[1], STDERR_FILENO);
    close(pipe_fds[1]);
    fn();
    _exit(0);
  }

  close(pipe_fds[1]);
  string child_output;
  char buffer[256];
  while (true) {
    const ssize_t bytes_read = read(pipe_fds[0], buffer, sizeof(buffer));
    if (bytes_read <= 0) {
      break;
    }
    child_output.append(buffer, static_cast<std::size_t>(bytes_read));
  }
  close(pipe_fds[0]);

  int status = 0;
  CHECK_MSG(waitpid(pid, &status, 0) == pid, label + ": waitpid() failed");
  if (g_test_failure_count != 0) {
    return;
  }

  const bool child_failed =
      (WIFEXITED(status) && WEXITSTATUS(status) != 0) || WIFSIGNALED(status);
  CHECK_MSG(child_failed, label + ": child unexpectedly succeeded");
  CHECK_MSG(child_output.find(needle) != string::npos,
            label + ": child output did not contain expected text: " +
                child_output);
#else
  (void)fn;
  (void)needle;
  testutil::PrintInfo(label + ": skipped negative child-process assertion on this platform");
#endif
}

basefold::FrobeniusPCSParams BuildFrobeniusParams(const ZZ &p,
                                                  const ZZ &base_modulus,
                                                  const ZZ_pX &F,
                                                  const ZZ_pE &alpha,
                                                  long ell,
                                                  long kappa) {
  basefold::FrobeniusPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = base_modulus;
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
  return basefold::FrobeniusPCSSetup(input);
}

basefold::FrobeniusPCSParams BuildFrobeniusParamsD0(const ZZ &p,
                                                    const ZZ &base_modulus,
                                                    const ZZ_pX &F,
                                                    const ZZ_pE &alpha,
                                                    long ell,
                                                    long kappa) {
  basefold::FrobeniusPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = base_modulus;
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42D0(p, alpha));
  return basefold::FrobeniusPCSSetup(input);
}

void TestFrobeniusPCSSetup_BuildsExpectedParams_GR42() {
  testutil::PrintInfo("Frobenius Phase 2: setup computes ell_prime and basis data");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);

  CHECK_EQ(params.ell_prime, 2L);
  CHECK_EQ(static_cast<long>(params.basis_data.normal_basis.beta.size()), 2L);
  CHECK_EQ(static_cast<long>(params.basis_data.normal_basis.alpha.size()), 2L);
  CHECK_EQ(basefold::Z2kPCSBackendMessageLength(params.backend), 4L);
  CHECK_EQ(basefold::Z2kPCSBackendPointDimension(params.backend), 2L);
}

void TestFrobeniusPCSSetup_RejectsBackendLengthMismatch() {
  testutil::PrintInfo("Frobenius Phase 2: setup rejects backend message-length mismatch");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F;
        SetCoeff(F, 2, 1);
        SetCoeff(F, 1, 1);
        SetCoeff(F, 0, 1);
        ZZ_pEPush ext_push(F);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        basefold::FrobeniusPCSSetupInput input;
        input.ell = 3;
        input.kappa = 1;
        input.base_modulus = modulus;
        input.extension_modulus = F;
        input.backend =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42D0(p, alpha));
        (void)basefold::FrobeniusPCSSetup(input);
      },
      "backend message length must equal 2^(ell-kappa)",
      "TestFrobeniusPCSSetup_RejectsBackendLengthMismatch");
}

void TestFrobeniusPCSSetup_RejectsDegreeMismatch() {
  testutil::PrintInfo("Frobenius Phase 2: setup rejects kappa/degree mismatch");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F;
        SetCoeff(F, 2, 1);
        SetCoeff(F, 1, 1);
        SetCoeff(F, 0, 1);
        ZZ_pEPush ext_push(F);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        basefold::FrobeniusPCSSetupInput input;
        input.ell = 4;
        input.kappa = 2;
        input.base_modulus = modulus;
        input.extension_modulus = F;
        input.backend =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
        (void)basefold::FrobeniusPCSSetup(input);
      },
      "current ZZ_pE degree must equal 2^kappa",
      "TestFrobeniusPCSSetup_RejectsDegreeMismatch");
}

void TestPackZ2kTableToFrobeniusGREvals_RoundTripsDegree2() {
  testutil::PrintInfo("Frobenius Phase 2: packing round-trips via dual basis in GR(4,2)");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});

  const vec_ZZ_pE packed =
      basefold::PackZ2kTableToFrobeniusGREvals(params, t_table);
  CHECK_EQ(packed.length(), 4L);
  for (long w = 0; w < packed.length(); ++w) {
    const vec_ZZ_pE unpacked =
        basefold::DecomposeGRElementToBaseCoeffsFrobeniusBasis(params,
                                                               packed[w]);
    CHECK_EQ(unpacked.length(), 2L);
    for (long u = 0; u < unpacked.length(); ++u) {
      CHECK_EQ(unpacked[u], t_table[u + (w << params.kappa)]);
    }
  }
}

void TestPackZ2kTableToFrobeniusGREvals_RoundTripsDegree4() {
  testutil::PrintInfo("Frobenius Phase 2: packing round-trips via dual basis in GR(4,4)");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 4, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/4, /*kappa=*/2);
  const vec_ZZ_pE t_table = BuildBaseRingCoeffVector(
      {0, 1, 2, 3, 1, 2, 3, 0, 2, 0, 1, 3, 3, 1, 0, 2});

  const vec_ZZ_pE packed =
      basefold::PackZ2kTableToFrobeniusGREvals(params, t_table);
  CHECK_EQ(packed.length(), 4L);
  for (long w = 0; w < packed.length(); ++w) {
    const vec_ZZ_pE unpacked =
        basefold::DecomposeGRElementToBaseCoeffsFrobeniusBasis(params,
                                                               packed[w]);
    CHECK_EQ(unpacked.length(), 4L);
    for (long u = 0; u < unpacked.length(); ++u) {
      CHECK_EQ(unpacked[u], t_table[u + (w << params.kappa)]);
    }
  }
}

void TestPackZ2kTableToFrobeniusGREvals_RejectsNonBaseRingInputs() {
  testutil::PrintInfo("Frobenius Phase 2: packing rejects non-constant Z2k inputs");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F;
        SetCoeff(F, 2, 1);
        SetCoeff(F, 1, 1);
        SetCoeff(F, 0, 1);
        ZZ_pEPush ext_push(F);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        const basefold::FrobeniusPCSParams params =
            BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
        vec_ZZ_pE t_table =
            BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});
        t_table[3] = alpha;
        (void)basefold::PackZ2kTableToFrobeniusGREvals(params, t_table);
      },
      "must be a base-ring constant",
      "TestPackZ2kTableToFrobeniusGREvals_RejectsNonBaseRingInputs");
}

void TestFrobeniusPCSCommit_MatchesDirectBackendCommit() {
  testutil::PrintInfo("Frobenius Phase 2: compiler commitment matches direct backend commitment");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});

  const basefold::MerkleRoot compiler_commitment =
      basefold::FrobeniusPCSCommit(params, t_table);
  const basefold::FrobeniusPCSOuterCommitArtifacts outer =
      basefold::FrobeniusPCSBuildOuterCommitArtifacts(params, t_table);
  const basefold::MerkleRoot direct_backend_commitment =
      basefold::Z2kPCSBackendCommit(params.backend, outer.t_packed_monomial_coeffs);

  CHECK_EQ(compiler_commitment, direct_backend_commitment);
}

void TestFrobeniusPCSBuildCommitArtifacts_CachesPackedRepresentations() {
  testutil::PrintInfo("Frobenius Phase 2: commit artifacts cache packed representations and backend artifacts");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 2, 1, 0, 1, 2, 3, 0});

  const basefold::FrobeniusPCSCommitArtifacts artifacts =
      basefold::FrobeniusPCSBuildCommitArtifacts(params, t_table);
  const basefold::FrobeniusPCSOuterCommitArtifacts outer =
      basefold::FrobeniusPCSBuildOuterCommitArtifacts(params, t_table);
  const basefold::MerkleRoot direct_backend_commitment =
      basefold::Z2kPCSBackendCommit(params.backend, outer.t_packed_monomial_coeffs);

  CHECK_EQ(artifacts.t_packed_table, outer.t_packed_table);
  CHECK_EQ(artifacts.t_packed_monomial_coeffs, outer.t_packed_monomial_coeffs);
  CHECK_EQ(artifacts.commitment, direct_backend_commitment);
  CHECK_EQ(artifacts.backend_commit_artifacts.commitment,
           direct_backend_commitment);
}

void TestFrobeniusPCSProveEval_ProducesConsistentOuterProof() {
  testutil::PrintInfo("Frobenius Phase 3: prove path matches Protocol 2 outer proof checks");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});
  const ZZ_pE omega = params.basis_data.teichmuller_generator;
  const vector<ZZ_pE> z = {omega + testutil::ConstZZpE(1), omega,
                           omega + testutil::ConstZZpE(1)};
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, params.ell, z);

  const basefold::FrobeniusPCSCommitArtifacts artifacts =
      basefold::FrobeniusPCSBuildCommitArtifacts(params, t_table);
  const basefold::FrobeniusPCSEvalProof proof =
      basefold::FrobeniusPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK_EQ(static_cast<long>(proof.s_by_i.size()), Pow2ForTest(params.kappa));
  CHECK_EQ(static_cast<long>(proof.h_by_level.size()), params.ell_prime);

  const vector<ZZ_pE> r_suffix(z.begin() + params.kappa, z.end());
  const vector<ZZ_pE> recovered_partials =
      RecoverPartialsFromProof(params, proof.s_by_i);
  const vector<ZZ_pE> direct_partials =
      DirectPartialEvaluations(t_table, params.kappa, params.ell_prime, r_suffix);
  CHECK_EQ(recovered_partials, direct_partials);

  const FrobeniusChallengeTrace trace = ReplayFrobeniusChallenges(
      artifacts.commitment, z, claimed_s, proof, params.kappa, params.ell_prime);
  ZZ_pE initial_claim = ZZ_pE(0);
  for (long i = 0; i < static_cast<long>(proof.s_by_i.size()); ++i) {
    initial_claim += proof.s_by_i[static_cast<std::size_t>(i)] *
                     basefold::EqPolynomial(
                         trace.rprime_prefix,
                         BooleanPointFromIndex(i, params.kappa));
  }
  CHECK(basefold::CheckProductSumcheckChain(initial_claim, proof.h_by_level,
                                            trace.rprime_suffix));

  CHECK_EQ(proof.t_star, EvalFromBooleanTable(artifacts.t_packed_table,
                                              params.ell_prime,
                                              trace.rprime_suffix));

  ZZ_pE g_star = ZZ_pE(0);
  for (long i = 0; i < Pow2ForTest(params.kappa); ++i) {
    const vector<ZZ_pE> sigma_point =
        ApplySigmaPowerToPoint(params, r_suffix, i);
    g_star += basefold::EqPolynomial(
                  trace.rprime_prefix, BooleanPointFromIndex(i, params.kappa)) *
              basefold::EqPolynomial(trace.rprime_suffix, sigma_point);
  }
  const ZZ_pE final_claim = proof.h_by_level.empty()
                                ? initial_claim
                                : proof.h_by_level[0].Eval(trace.rprime_suffix[0]);
  CHECK_EQ(final_claim, proof.t_star * g_star);

  CHECK(basefold::Z2kPCSBackendVerifyEval(
      params.backend, artifacts.commitment, trace.rprime_suffix, proof.t_star,
      /*num_queries=*/2, proof.backend_proof));
}

void TestFrobeniusPCSProveEval_DirectAndArtifactPathsAgreeOnOuterMessages() {
  testutil::PrintInfo("Frobenius Phase 3: direct and artifact-backed prove paths agree on outer messages");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const ZZ_pE omega = params.basis_data.teichmuller_generator;
  const vector<ZZ_pE> z = {omega, omega + testutil::ConstZZpE(1),
                           testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, params.ell, z);

  const basefold::FrobeniusPCSCommitArtifacts artifacts =
      basefold::FrobeniusPCSBuildCommitArtifacts(params, t_table);
  const basefold::FrobeniusPCSEvalProof direct =
      basefold::FrobeniusPCSProveEval(params, t_table, z, claimed_s,
                                      /*num_queries=*/2);
  const basefold::FrobeniusPCSEvalProof cached =
      basefold::FrobeniusPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK_EQ(direct.s_by_i, cached.s_by_i);
  CHECK_EQ(static_cast<long>(direct.h_by_level.size()),
           static_cast<long>(cached.h_by_level.size()));
  for (long i = 0; i < static_cast<long>(direct.h_by_level.size()); ++i) {
    CHECK_EQ(direct.h_by_level[static_cast<std::size_t>(i)].a0,
             cached.h_by_level[static_cast<std::size_t>(i)].a0);
    CHECK_EQ(direct.h_by_level[static_cast<std::size_t>(i)].a1,
             cached.h_by_level[static_cast<std::size_t>(i)].a1);
    CHECK_EQ(direct.h_by_level[static_cast<std::size_t>(i)].a2,
             cached.h_by_level[static_cast<std::size_t>(i)].a2);
  }
  CHECK_EQ(direct.t_star, cached.t_star);
  CHECK_EQ(basefold::Z2kPCSBackendEvalProofSizeBytes(params.backend,
                                                     direct.backend_proof),
           basefold::Z2kPCSBackendEvalProofSizeBytes(params.backend,
                                                     cached.backend_proof));
}

void TestFrobeniusPCSVerifyEval_AcceptsHonestProof() {
  testutil::PrintInfo("Frobenius Phase 3: verifier accepts an honest proof");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});
  const ZZ_pE omega = params.basis_data.teichmuller_generator;
  const vector<ZZ_pE> z = {omega + testutil::ConstZZpE(1), omega,
                           omega + testutil::ConstZZpE(1)};
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, params.ell, z);

  const basefold::FrobeniusPCSCommitArtifacts artifacts =
      basefold::FrobeniusPCSBuildCommitArtifacts(params, t_table);
  const basefold::FrobeniusPCSEvalProof proof =
      basefold::FrobeniusPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK(basefold::FrobeniusPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s, /*num_queries=*/2, proof));
}

void TestFrobeniusPCSVerifyEval_AcceptsHonestProofFromDirectProvePath() {
  testutil::PrintInfo("Frobenius Phase 3: direct prove path verifies end-to-end");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const ZZ_pE omega = params.basis_data.teichmuller_generator;
  const vector<ZZ_pE> z = {omega, omega + testutil::ConstZZpE(1),
                           testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, params.ell, z);

  const basefold::MerkleRoot commitment =
      basefold::FrobeniusPCSCommit(params, t_table);
  const basefold::FrobeniusPCSEvalProof proof =
      basefold::FrobeniusPCSProveEval(params, t_table, z, claimed_s,
                                      /*num_queries=*/2);

  CHECK(basefold::FrobeniusPCSVerifyEval(
      params, commitment, z, claimed_s, /*num_queries=*/2, proof));
}

void TestFrobeniusPCSPaperAPI_AcceptsHonestProof() {
  testutil::PrintInfo("Frobenius Phase 3: staged setup/commit/prove/verify API matches the direct flow");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  basefold::FrobeniusPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
  const basefold::FrobeniusPCSSetupOutput api =
      basefold::FrobeniusPCSSetupProtocol(input);

  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const ZZ_pE omega = api.params.basis_data.teichmuller_generator;
  const vector<ZZ_pE> z = {omega, omega + testutil::ConstZZpE(1),
                           testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, api.params.ell, z);

  const basefold::FrobeniusPCSCommittedWitness committed =
      api.prover.Commit(t_table);
  CHECK_EQ(committed.commitment, basefold::FrobeniusPCSCommit(api.params, t_table));
  CHECK_EQ(committed.commitment, committed.commit_artifacts.commitment);

  const basefold::FrobeniusPCSEvalProof proof =
      api.prover.Prove(committed, z, claimed_s, /*num_queries=*/2);
  CHECK(api.verifier.Verify(committed.commitment, z, claimed_s,
                            /*num_queries=*/2, proof));
}

void TestFrobeniusPCSVerifyEval_RejectsTampering() {
  testutil::PrintInfo("Frobenius Phase 3: verifier rejects wrong claim, outer tampering, and backend tampering");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const ZZ_pE omega = params.basis_data.teichmuller_generator;
  const vector<ZZ_pE> z = {omega, omega + testutil::ConstZZpE(1),
                           testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, params.ell, z);

  const basefold::FrobeniusPCSCommitArtifacts artifacts =
      basefold::FrobeniusPCSBuildCommitArtifacts(params, t_table);
  const basefold::FrobeniusPCSEvalProof proof =
      basefold::FrobeniusPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK(basefold::FrobeniusPCSVerifyEval(params, artifacts.commitment, z,
                                         claimed_s, /*num_queries=*/2, proof));

  CHECK(!basefold::FrobeniusPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s + testutil::ConstZZpE(1),
      /*num_queries=*/2, proof));

  basefold::FrobeniusPCSEvalProof tampered_s_by_i = proof;
  tampered_s_by_i.s_by_i[0] += testutil::ConstZZpE(1);
  CHECK(!basefold::FrobeniusPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      tampered_s_by_i));

  basefold::FrobeniusPCSEvalProof tampered_h = proof;
  tampered_h.h_by_level[0].a0 += testutil::ConstZZpE(1);
  CHECK(!basefold::FrobeniusPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      tampered_h));

  basefold::FrobeniusPCSEvalProof tampered_t_star = proof;
  tampered_t_star.t_star += testutil::ConstZZpE(1);
  CHECK(!basefold::FrobeniusPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      tampered_t_star));

  basefold::FrobeniusPCSEvalProof tampered_backend = proof;
  tampered_backend.backend_proof =
      MutateBaseFoldBackendSubproof(proof.backend_proof);
  CHECK(!basefold::FrobeniusPCSVerifyEval(
      params, artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      tampered_backend));
}

void TestFrobeniusPCSVerifyEval_DimensionZeroUsesNoSumcheckRounds() {
  testutil::PrintInfo("Frobenius Phase 3: verifier handles ell_prime=0 as a zero-round sumcheck");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParamsD0(p, modulus, F, alpha, /*ell=*/1, /*kappa=*/1);
  const vec_ZZ_pE t_table = BuildBaseRingCoeffVector({1, 3});
  const ZZ_pE omega = params.basis_data.teichmuller_generator;
  const vector<ZZ_pE> z = {omega + testutil::ConstZZpE(1)};
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, params.ell, z);

  const basefold::FrobeniusPCSCommitArtifacts artifacts =
      basefold::FrobeniusPCSBuildCommitArtifacts(params, t_table);
  const basefold::FrobeniusPCSEvalProof proof =
      basefold::FrobeniusPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);

  CHECK(proof.h_by_level.empty());
  CHECK(basefold::FrobeniusPCSVerifyEval(params, artifacts.commitment, z,
                                         claimed_s, /*num_queries=*/2, proof));
}

void TestFrobeniusPCSOuterProveVerify_AcceptsHonestProof() {
  testutil::PrintInfo("Frobenius Phase 3: prover/verifier accepts an honest outer-only proof");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const ZZ_pE omega = params.basis_data.teichmuller_generator;
  const vector<ZZ_pE> z = {omega, omega + testutil::ConstZZpE(1),
                           testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, params.ell, z);

  const basefold::FrobeniusPCSCommitArtifacts composed_artifacts =
      basefold::FrobeniusPCSBuildCommitArtifacts(params, t_table);
  const basefold::FrobeniusPCSOuterCommitArtifacts outer_artifacts =
      basefold::FrobeniusPCSBuildOuterCommitArtifacts(params, t_table);
  CHECK_EQ(outer_artifacts.t_packed_table, composed_artifacts.t_packed_table);
  CHECK_EQ(outer_artifacts.t_packed_monomial_coeffs,
           composed_artifacts.t_packed_monomial_coeffs);

  const basefold::FrobeniusPCSOuterEvalProof outer_proof =
      basefold::FrobeniusPCSProveOuterEvalFromCommitArtifacts(
          params, t_table, composed_artifacts.commitment, z, claimed_s,
          /*num_queries=*/2, outer_artifacts);
  CHECK(basefold::FrobeniusPCSVerifyOuterEval(
      params, composed_artifacts.commitment, z, claimed_s, /*num_queries=*/2,
      outer_proof));

  const basefold::FrobeniusPCSEvalProof composed_proof =
      basefold::FrobeniusPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, composed_artifacts);
  CHECK_EQ(outer_proof.s_by_i, composed_proof.s_by_i);
  CHECK_EQ(outer_proof.h_by_level.size(), composed_proof.h_by_level.size());
  for (long i = 0; i < static_cast<long>(outer_proof.h_by_level.size()); ++i) {
    CHECK_EQ(outer_proof.h_by_level[static_cast<std::size_t>(i)].a0,
             composed_proof.h_by_level[static_cast<std::size_t>(i)].a0);
    CHECK_EQ(outer_proof.h_by_level[static_cast<std::size_t>(i)].a1,
             composed_proof.h_by_level[static_cast<std::size_t>(i)].a1);
    CHECK_EQ(outer_proof.h_by_level[static_cast<std::size_t>(i)].a2,
             composed_proof.h_by_level[static_cast<std::size_t>(i)].a2);
  }
  CHECK_EQ(outer_proof.t_star, composed_proof.t_star);
}

void TestFrobeniusProofSerialize_ComposedSizeMatchesBytes() {
  testutil::PrintInfo("Frobenius Phase 4: serializer bytes match outer and composed proof-size accounting");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FrobeniusPCSParams params =
      BuildFrobeniusParams(p, modulus, F, alpha, /*ell=*/3, /*kappa=*/1);
  const vec_ZZ_pE t_table =
      BuildBaseRingCoeffVector({3, 1, 0, 2, 1, 2, 3, 0});
  const vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), alpha,
                           testutil::ConstZZpE(3)};
  const ZZ_pE claimed_s = EvalFromBooleanTable(t_table, params.ell, z);

  const basefold::FrobeniusPCSCommitArtifacts artifacts =
      basefold::FrobeniusPCSBuildCommitArtifacts(params, t_table);
  const basefold::FrobeniusPCSEvalProof proof =
      basefold::FrobeniusPCSProveEvalFromCommitArtifacts(
          params, t_table, z, claimed_s, /*num_queries=*/2, artifacts);
  const basefold::FrobeniusPCSOuterEvalProof outer_proof =
      basefold::FrobeniusPCSProveOuterEvalFromCommitArtifacts(
          params, t_table, artifacts.commitment, z, claimed_s,
          /*num_queries=*/2,
          basefold::FrobeniusPCSBuildOuterCommitArtifacts(params, t_table));

  const std::uint64_t field_elem_bytes = FixedFieldElementBytesForCurrentContext();
  const std::uint64_t expected_outer_bytes =
      1U + 8U +
      static_cast<std::uint64_t>(outer_proof.s_by_i.size()) * field_elem_bytes +
      8U + static_cast<std::uint64_t>(outer_proof.h_by_level.size()) * 3U *
                field_elem_bytes +
      field_elem_bytes;

  const basefold::Bytes outer_bytes =
      basefold::SerializeFrobeniusPCSOuterProofFixedBytes(params, outer_proof);
  const std::uint64_t outer_size =
      basefold::FrobeniusPCSOuterProofSizeBytes(params, outer_proof);
  CHECK_EQ(outer_bytes.size(), static_cast<std::size_t>(outer_size));
  CHECK_EQ(outer_size, expected_outer_bytes);

  const basefold::Bytes backend_bytes =
      basefold::Z2kPCSBackendSerializeEvalProof(params.backend,
                                                proof.backend_proof);
  const std::uint64_t backend_size =
      basefold::Z2kPCSBackendEvalProofSizeBytes(params.backend,
                                                proof.backend_proof);
  CHECK_EQ(backend_bytes.size(), static_cast<std::size_t>(backend_size));

  const basefold::Bytes composed_bytes =
      basefold::SerializeFrobeniusPCSEvalProofFixedBytes(params, proof);
  const std::uint64_t composed_size =
      basefold::FrobeniusPCSEvalProofSizeBytes(params, proof);
  CHECK_EQ(composed_bytes.size(), static_cast<std::size_t>(composed_size));
  CHECK_EQ(composed_size, outer_size + 8U + backend_size);
  CHECK(std::equal(outer_bytes.begin(), outer_bytes.end(),
                   composed_bytes.begin()));
}

}  // namespace

int main() {
  try {
    RUN_TEST(TestFrobeniusPCSSetup_BuildsExpectedParams_GR42);
    RUN_TEST(TestFrobeniusPCSSetup_RejectsBackendLengthMismatch);
    RUN_TEST(TestFrobeniusPCSSetup_RejectsDegreeMismatch);
    RUN_TEST(TestPackZ2kTableToFrobeniusGREvals_RoundTripsDegree2);
    RUN_TEST(TestPackZ2kTableToFrobeniusGREvals_RoundTripsDegree4);
    RUN_TEST(TestPackZ2kTableToFrobeniusGREvals_RejectsNonBaseRingInputs);
    RUN_TEST(TestFrobeniusPCSCommit_MatchesDirectBackendCommit);
    RUN_TEST(TestFrobeniusPCSBuildCommitArtifacts_CachesPackedRepresentations);
    RUN_TEST(TestFrobeniusPCSProveEval_ProducesConsistentOuterProof);
    RUN_TEST(TestFrobeniusPCSProveEval_DirectAndArtifactPathsAgreeOnOuterMessages);
    RUN_TEST(TestFrobeniusPCSVerifyEval_AcceptsHonestProof);
    RUN_TEST(TestFrobeniusPCSVerifyEval_AcceptsHonestProofFromDirectProvePath);
    RUN_TEST(TestFrobeniusPCSPaperAPI_AcceptsHonestProof);
    RUN_TEST(TestFrobeniusPCSVerifyEval_RejectsTampering);
    RUN_TEST(TestFrobeniusPCSVerifyEval_DimensionZeroUsesNoSumcheckRounds);
    RUN_TEST(TestFrobeniusPCSOuterProveVerify_AcceptsHonestProof);
    RUN_TEST(TestFrobeniusProofSerialize_ComposedSizeMatchesBytes);
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
