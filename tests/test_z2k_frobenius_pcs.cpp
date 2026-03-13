#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <exception>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "Compiler/Z2k/BaseFoldBackendAdapter.hpp"
#include "Compiler/Z2k/FrobeniusPCS.hpp"
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
