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

#include "BaseFold/Multilinear.hpp"
#include "BaseFold/Z2kPCSBackend.hpp"
#include "BaseFold/Z2kRingSwitchPCS.hpp"
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

basefold::FoldableCodeParams BuildParamsGR42Variant(const ZZ &p,
                                                    const ZZ_pE &alpha) {
  basefold::FoldableCodeParams params = BuildParamsGR42(p, alpha);
  const ZZ_pE one = testutil::ConstZZpE(1);
  params.zeta = alpha + one;
  params.G0[0][1] = one;
  params.diag_T[0][1] = one;
  params.diag_T[1][0] = one;
  params.diag_T[1][2] = alpha + one;
  return params;
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

void TestBaseFoldBackendAdapter_Smoke() {
  testutil::PrintInfo("Z2k backend: BaseFold adapter commit/prove/verify works");

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

  const basefold::FoldableCodeParams params = BuildParamsGR42(p, alpha);
  const basefold::Z2kPCSBackendHandle backend =
      basefold::MakeBaseFoldZ2kPCSBackend(params);

  CHECK_EQ(string(basefold::Z2kPCSBackendName(backend)), string("basefold"));
  CHECK_EQ(basefold::Z2kPCSBackendMessageLength(backend), 4L);
  CHECK_EQ(basefold::Z2kPCSBackendPointDimension(backend), 2L);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(4);
  f_coeffs[0] = testutil::ConstZZpE(0);
  f_coeffs[1] = testutil::ConstZZpE(1);
  f_coeffs[2] = alpha;
  f_coeffs[3] = testutil::ConstZZpE(2);

  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                testutil::ConstZZpE(3)};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
  const basefold::MerkleRoot commitment =
      basefold::Z2kPCSBackendCommit(backend, f_coeffs);
  const basefold::Z2kPCSBackendCommitArtifacts commit_artifacts =
      basefold::Z2kPCSBackendBuildCommitArtifacts(backend, f_coeffs);
  const basefold::Z2kPCSBackendEvalProof proof = basefold::Z2kPCSBackendProveEval(
      backend, f_coeffs, z, y, /*num_queries=*/3, &commit_artifacts);

  CHECK(basefold::Z2kPCSBackendVerifyEval(backend, commitment, z, y,
                                          /*num_queries=*/3, proof));
  CHECK_GT(basefold::Z2kPCSBackendEvalProofSizeBytes(backend, proof), 0U);
}

void TestBaseFoldBackendAdapter_RejectsSameDegreeContextSwitch() {
  testutil::PrintInfo("Z2k backend: handle rejects same-degree context switches");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F1;
        SetCoeff(F1, 2, 1);
        SetCoeff(F1, 1, 1);
        SetCoeff(F1, 0, 1);
        ZZ_pEPush ext_push1(F1);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        const basefold::Z2kPCSBackendHandle backend =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

        ZZ_pX F2;
        SetCoeff(F2, 2, 1);
        SetCoeff(F2, 1, 3);
        SetCoeff(F2, 0, 1);
        ZZ_pEPush ext_push2(F2);

        basefold::RingSwitchPCSSetupInput input;
        input.ell = 3;
        input.kappa = 1;
        input.base_modulus = modulus;
        input.extension_modulus = F2;
        input.alpha_basis = basefold::ActivePolynomialBasisDescriptor();
        input.beta_basis = basefold::ActivePolynomialBasisDescriptor();
        input.backend = backend;
        (void)basefold::RingSwitchPCSSetup(input);
      },
      "current ZZ_pE modulus does not match backend context",
      "TestBaseFoldBackendAdapter_RejectsSameDegreeContextSwitch");
}

void TestBaseFoldBackendAdapter_RejectsArtifactsFromDifferentHandle() {
  testutil::PrintInfo("Z2k backend: commit artifacts are bound to one backend instance");

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

        const basefold::Z2kPCSBackendHandle backend_a =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
        const basefold::Z2kPCSBackendHandle backend_b =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42Variant(p, alpha));

        vec_ZZ_pE f_coeffs;
        f_coeffs.SetLength(4);
        f_coeffs[0] = testutil::ConstZZpE(0);
        f_coeffs[1] = testutil::ConstZZpE(1);
        f_coeffs[2] = alpha;
        f_coeffs[3] = testutil::ConstZZpE(2);
        const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                      testutil::ConstZZpE(3)};
        const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

        const basefold::Z2kPCSBackendCommitArtifacts artifacts_a =
            basefold::Z2kPCSBackendBuildCommitArtifacts(backend_a, f_coeffs);
        (void)basefold::Z2kPCSBackendProveEval(backend_b, f_coeffs, z, y,
                                               /*num_queries=*/3, &artifacts_a);
      },
      "commit artifacts must belong to the same backend params instance",
      "TestBaseFoldBackendAdapter_RejectsArtifactsFromDifferentHandle");
}

void TestBaseFoldBackendAdapter_RejectsProofReuseAcrossHandles() {
  testutil::PrintInfo("Z2k backend: proofs are bound to one backend instance");

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

        const basefold::Z2kPCSBackendHandle backend_a =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
        const basefold::Z2kPCSBackendHandle backend_b =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42Variant(p, alpha));

        vec_ZZ_pE f_coeffs;
        f_coeffs.SetLength(4);
        f_coeffs[0] = testutil::ConstZZpE(0);
        f_coeffs[1] = testutil::ConstZZpE(1);
        f_coeffs[2] = alpha;
        f_coeffs[3] = testutil::ConstZZpE(2);
        const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                      testutil::ConstZZpE(3)};
        const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
        const basefold::MerkleRoot commitment =
            basefold::Z2kPCSBackendCommit(backend_a, f_coeffs);
        const basefold::Z2kPCSBackendEvalProof proof_a =
            basefold::Z2kPCSBackendProveEval(backend_a, f_coeffs, z, y,
                                             /*num_queries=*/3);
        (void)basefold::Z2kPCSBackendVerifyEval(backend_b, commitment, z, y,
                                                /*num_queries=*/3, proof_a);
      },
      "proof must belong to the same backend params instance",
      "TestBaseFoldBackendAdapter_RejectsProofReuseAcrossHandles");
}

void TestBaseFoldBackendAdapter_RejectsProofSizeAfterContextSwitch() {
  testutil::PrintInfo("Z2k backend: proof-size counting rejects context drift");

  ExpectChildFailureContains(
      [&]() {
        const ZZ p = to_ZZ(2);
        const ZZ modulus = to_ZZ(4);
        ZZ_pPush mod_push(modulus);

        ZZ_pX F1;
        SetCoeff(F1, 2, 1);
        SetCoeff(F1, 1, 1);
        SetCoeff(F1, 0, 1);
        ZZ_pEPush ext_push1(F1);

        ZZ_pX xpoly;
        SetCoeff(xpoly, 1, 1);
        ZZ_pE alpha;
        conv(alpha, xpoly);

        const basefold::Z2kPCSBackendHandle backend =
            basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
        vec_ZZ_pE f_coeffs;
        f_coeffs.SetLength(4);
        f_coeffs[0] = testutil::ConstZZpE(0);
        f_coeffs[1] = testutil::ConstZZpE(1);
        f_coeffs[2] = alpha;
        f_coeffs[3] = testutil::ConstZZpE(2);
        const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                      testutil::ConstZZpE(3)};
        const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
        const basefold::Z2kPCSBackendEvalProof proof =
            basefold::Z2kPCSBackendProveEval(backend, f_coeffs, z, y,
                                             /*num_queries=*/3);

        ZZ_pX F2;
        SetCoeff(F2, 2, 1);
        SetCoeff(F2, 1, 3);
        SetCoeff(F2, 0, 1);
        ZZ_pEPush ext_push2(F2);

        (void)basefold::Z2kPCSBackendEvalProofSizeBytes(backend, proof);
      },
      "current ZZ_pE modulus does not match backend context",
      "TestBaseFoldBackendAdapter_RejectsProofSizeAfterContextSwitch");
}

void TestActivePolynomialBasisDescriptor_UsesCurrentDegree() {
  testutil::PrintInfo("Ring-switch setup: active polynomial basis reflects ZZ_pE degree");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  const basefold::RingSwitchBasisDescriptor basis =
      basefold::ActivePolynomialBasisDescriptor();
  CHECK(basis.kind == basefold::RingSwitchBasisKind::kPolynomial);
  CHECK_EQ(basis.dimension, 2L);
}

void TestValidateCurrentZ2kRingContext_SucceedsForGR42() {
  testutil::PrintInfo("Ring-switch setup: current context validation accepts GR(4,2)");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  basefold::ValidateCurrentZ2kRingContextOrThrow(modulus, F, /*kappa=*/1);
}

void TestValidateCurrentZ2kRingContext_RejectsReducibleMod2Polynomial() {
  testutil::PrintInfo("Ring-switch setup: current context validation rejects non-basic extension");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F_bad;
  SetCoeff(F_bad, 2, 1);
  SetCoeff(F_bad, 0, 1);
  ZZ_pEPush ext_push(F_bad);

  ExpectChildFailureContains(
      [&]() {
        basefold::ValidateCurrentZ2kRingContextOrThrow(modulus, F_bad,
                                                       /*kappa=*/1);
      },
      "basic irreducible modulo 2",
      "TestValidateCurrentZ2kRingContext_RejectsReducibleMod2Polynomial");
}

void TestRingSwitchSetup_Succeeds() {
  testutil::PrintInfo("Ring-switch setup: lightweight Setup validates and returns params");

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

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.alpha_basis = basefold::ActivePolynomialBasisDescriptor();
  input.beta_basis = basefold::ActivePolynomialBasisDescriptor();
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  const basefold::RingSwitchPCSParams params = basefold::RingSwitchPCSSetup(input);
  CHECK_EQ(params.ell, 3L);
  CHECK_EQ(params.kappa, 1L);
  CHECK_EQ(params.ell_prime, 2L);
  CHECK(params.alpha_basis.kind == basefold::RingSwitchBasisKind::kPolynomial);
  CHECK(params.beta_basis.kind == basefold::RingSwitchBasisKind::kPolynomial);
  CHECK_EQ(params.alpha_basis.dimension, 2L);
  CHECK_EQ(params.beta_basis.dimension, 2L);
  CHECK_EQ(basefold::Z2kPCSBackendPointDimension(params.backend), 2L);
}

void TestRingSwitchSetup_RejectsMismatchedBaseModulus() {
  testutil::PrintInfo("Ring-switch setup: Setup rejects mismatched base modulus");

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

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = to_ZZ(8);
  input.extension_modulus = F;
  input.alpha_basis = basefold::ActivePolynomialBasisDescriptor();
  input.beta_basis = basefold::ActivePolynomialBasisDescriptor();
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "current ZZ_p modulus does not match base_modulus",
      "TestRingSwitchSetup_RejectsMismatchedBaseModulus");
}

void TestRingSwitchSetup_RejectsMismatchedExtensionModulus() {
  testutil::PrintInfo("Ring-switch setup: Setup rejects mismatched extension modulus");

  const ZZ p = to_ZZ(2);
  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  ZZ_pX F_bad;
  SetCoeff(F_bad, 2, 1);
  SetCoeff(F_bad, 0, 1);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F_bad;
  input.alpha_basis = basefold::ActivePolynomialBasisDescriptor();
  input.beta_basis = basefold::ActivePolynomialBasisDescriptor();
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "current ZZ_pE modulus does not match extension_modulus",
      "TestRingSwitchSetup_RejectsMismatchedExtensionModulus");
}

void TestRingSwitchSetup_RejectsNonPolynomialBasis() {
  testutil::PrintInfo("Ring-switch setup: Setup rejects unsupported basis kinds");

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

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 3;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.alpha_basis.kind = static_cast<basefold::RingSwitchBasisKind>(1);
  input.alpha_basis.dimension = 2;
  input.beta_basis = basefold::ActivePolynomialBasisDescriptor();
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "alpha_basis must be the active polynomial basis",
      "TestRingSwitchSetup_RejectsNonPolynomialBasis");
}

void TestRingSwitchSetup_RejectsBackendDimensionMismatch() {
  testutil::PrintInfo("Ring-switch setup: Setup rejects backend dimension mismatch");

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

  basefold::RingSwitchPCSSetupInput input;
  input.ell = 4;
  input.kappa = 1;
  input.base_modulus = modulus;
  input.extension_modulus = F;
  input.alpha_basis = basefold::ActivePolynomialBasisDescriptor();
  input.beta_basis = basefold::ActivePolynomialBasisDescriptor();
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));

  ExpectChildFailureContains(
      [&]() { (void)basefold::RingSwitchPCSSetup(input); },
      "backend message length must equal 2^(ell-kappa)",
      "TestRingSwitchSetup_RejectsBackendDimensionMismatch");
}

}  // namespace

int main() {
  try {
    RUN_TEST(TestBaseFoldBackendAdapter_Smoke);
    RUN_TEST(TestBaseFoldBackendAdapter_RejectsSameDegreeContextSwitch);
    RUN_TEST(TestBaseFoldBackendAdapter_RejectsArtifactsFromDifferentHandle);
    RUN_TEST(TestBaseFoldBackendAdapter_RejectsProofReuseAcrossHandles);
    RUN_TEST(TestBaseFoldBackendAdapter_RejectsProofSizeAfterContextSwitch);
    RUN_TEST(TestActivePolynomialBasisDescriptor_UsesCurrentDegree);
    RUN_TEST(TestValidateCurrentZ2kRingContext_SucceedsForGR42);
    RUN_TEST(TestValidateCurrentZ2kRingContext_RejectsReducibleMod2Polynomial);
    RUN_TEST(TestRingSwitchSetup_Succeeds);
    RUN_TEST(TestRingSwitchSetup_RejectsMismatchedBaseModulus);
    RUN_TEST(TestRingSwitchSetup_RejectsMismatchedExtensionModulus);
    RUN_TEST(TestRingSwitchSetup_RejectsNonPolynomialBasis);
    RUN_TEST(TestRingSwitchSetup_RejectsBackendDimensionMismatch);
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
