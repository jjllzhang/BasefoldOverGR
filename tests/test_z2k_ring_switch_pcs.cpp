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
#include "BaseFold/Sumcheck.hpp"
#include "BaseFold/Z2kPCSBackend.hpp"
#include "BaseFold/Z2kRingSwitchPCS.hpp"
#include "GaloisRing/utils.hpp"
#include "tests/test_common.hpp"

using NTL::conv;
using NTL::clear;
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

long Pow2ForTest(long exponent) {
  CHECK_MSG(exponent >= 0, "Pow2ForTest: exponent must be non-negative");
  if (exponent < 0) {
    return 0;
  }
  return 1L << exponent;
}

vec_ZZ_pE BuildBaseRingCoeffVector(const std::vector<long> &coeffs) {
  vec_ZZ_pE out;
  out.SetLength(static_cast<long>(coeffs.size()));
  for (long i = 0; i < static_cast<long>(coeffs.size()); ++i) {
    out[i] = testutil::ConstZZpE(coeffs[static_cast<std::size_t>(i)]);
  }
  return out;
}

std::vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  std::vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] =
        testutil::ConstZZpE((index >> i) & 1L);
  }
  return point;
}

ZZ_pE PolynomialBasisElement(long degree) {
  ZZ_pX poly;
  SetCoeff(poly, degree, 1);
  ZZ_pE out;
  conv(out, poly);
  return out;
}

ZZ_pE SumOfPointwiseProducts(const vec_ZZ_pE &f_table,
                             const vec_ZZ_pE &g_table) {
  CHECK_EQ(f_table.length(), g_table.length());
  ZZ_pE acc;
  clear(acc);
  for (long i = 0; i < f_table.length(); ++i) {
    acc += f_table[i] * g_table[i];
  }
  return acc;
}

basefold::RingSwitchPCSParams BuildRingSwitchParamsGR42(long ell, long kappa,
                                                        const ZZ &base_modulus,
                                                        const ZZ_pX &F,
                                                        const ZZ &p,
                                                        const ZZ_pE &alpha) {
  basefold::RingSwitchPCSSetupInput input;
  input.ell = ell;
  input.kappa = kappa;
  input.base_modulus = base_modulus;
  input.extension_modulus = F;
  input.alpha_basis = basefold::ActivePolynomialBasisDescriptor();
  input.beta_basis = basefold::ActivePolynomialBasisDescriptor();
  input.backend = basefold::MakeBaseFoldZ2kPCSBackend(BuildParamsGR42(p, alpha));
  return basefold::RingSwitchPCSSetup(input);
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

void TestPackZ2kCoeffsToGREvals_RoundTripsSmallExample() {
  testutil::PrintInfo("Ring-switch WP1: packing follows the polynomial-basis block layout");

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

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_coeffs =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});

  const vec_ZZ_pE packed = basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
  CHECK_EQ(packed.length(), 4L);
  CHECK_EQ(packed[0], LongVecToZZpE({0, 1}));
  CHECK_EQ(packed[1], LongVecToZZpE({2, 3}));
  CHECK_EQ(packed[2], LongVecToZZpE({1, 0}));
  CHECK_EQ(packed[3], LongVecToZZpE({3, 2}));

  for (long w = 0; w < packed.length(); ++w) {
    const vec_ZZ_pE decomposed =
        basefold::DecomposeGRElementToBaseCoeffsPolynomialBasis(params,
                                                                packed[w]);
    CHECK_EQ(decomposed.length(), Pow2ForTest(params.kappa));
    for (long v = 0; v < decomposed.length(); ++v) {
      CHECK_EQ(decomposed[v], t_coeffs[v + (w << params.kappa)]);
    }
  }
}

void TestBooleanHypercubeTableToMonomialCoeffs_EvaluatesTableCorrectly() {
  testutil::PrintInfo("Ring-switch WP1: Boolean-table to monomial conversion matches table semantics");

  const ZZ modulus = to_ZZ(4);
  ZZ_pPush mod_push(modulus);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush ext_push(F);

  const vec_ZZ_pE table = BuildBaseRingCoeffVector({1, 2, 3, 0});
  const vec_ZZ_pE monomial =
      basefold::BooleanHypercubeTableToMonomialCoeffs(table);

  CHECK_EQ(basefold::EvalMultilinearMonomialCoeffs(
               monomial, BooleanPointFromIndex(/*index=*/0, /*dimension=*/2)),
           table[0]);
  CHECK_EQ(basefold::EvalMultilinearMonomialCoeffs(
               monomial, BooleanPointFromIndex(/*index=*/1, /*dimension=*/2)),
           table[1]);
  CHECK_EQ(basefold::EvalMultilinearMonomialCoeffs(
               monomial, BooleanPointFromIndex(/*index=*/2, /*dimension=*/2)),
           table[2]);
  CHECK_EQ(basefold::EvalMultilinearMonomialCoeffs(
               monomial, BooleanPointFromIndex(/*index=*/3, /*dimension=*/2)),
           table[3]);
}

void TestPackZ2kCoeffsToGREvals_RejectsNonBaseRingInputs() {
  testutil::PrintInfo("Ring-switch WP1: packing rejects non-constant Z2k inputs");

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

        const basefold::RingSwitchPCSParams params = BuildRingSwitchParamsGR42(
            /*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
        vec_ZZ_pE t_coeffs = BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});
        t_coeffs[3] = alpha;
        (void)basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
      },
      "must be a base-ring constant",
      "TestPackZ2kCoeffsToGREvals_RejectsNonBaseRingInputs");
}

void TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValues() {
  testutil::PrintInfo("Ring-switch WP1: A_{u||w} reconstructs suffix equalities");

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

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const std::vector<ZZ_pE> r_suffix = {alpha + testutil::ConstZZpE(1),
                                       testutil::ConstZZpE(3)};
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);

  CHECK_EQ(tensor.basis_dimension, Pow2ForTest(params.kappa));
  CHECK_EQ(tensor.ell_prime, params.ell_prime);
  CHECK_EQ(tensor.a_by_u_then_w.length(), Pow2ForTest(params.ell));
  CHECK_EQ(tensor.r_table.length(), Pow2ForTest(params.ell));
  CHECK_EQ(tensor.r_monomial_coeffs.length(), Pow2ForTest(params.ell));
  CHECK_EQ(tensor.r_monomial_coeffs,
           basefold::BooleanHypercubeTableToMonomialCoeffs(tensor.r_table));

  const ZZ_pE x = PolynomialBasisElement(1);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long w = 0; w < num_w; ++w) {
    ZZ_pE reconstructed = tensor.a_by_u_then_w[w];
    reconstructed += tensor.a_by_u_then_w[num_w + w] * x;
    const ZZ_pE expected =
        basefold::EqPolynomial(r_suffix, BooleanPointFromIndex(w, params.ell_prime));
    CHECK_EQ(reconstructed, expected);
  }
}

void TestBuildRingSwitchComponentTensor_RecoversPartialEvaluations() {
  testutil::PrintInfo("Ring-switch WP1: s_u recovers the Appendix C.1 partial evaluations");

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

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const vec_ZZ_pE t_coeffs =
      BuildBaseRingCoeffVector({0, 1, 2, 3, 1, 0, 3, 2});
  const vec_ZZ_pE packed = basefold::PackZ2kCoeffsToGREvals(params, t_coeffs);
  const std::vector<ZZ_pE> r_suffix = {alpha + testutil::ConstZZpE(1),
                                       testutil::ConstZZpE(3)};
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);

  vec_ZZ_pE s_by_u;
  s_by_u.SetLength(tensor.basis_dimension);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long u = 0; u < tensor.basis_dimension; ++u) {
    ZZ_pE acc;
    clear(acc);
    for (long w = 0; w < num_w; ++w) {
      acc += tensor.a_by_u_then_w[u * num_w + w] * packed[w];
    }
    s_by_u[u] = acc;
  }

  const ZZ_pE x = PolynomialBasisElement(1);
  for (long v = 0; v < tensor.basis_dimension; ++v) {
    ZZ_pE recovered_partial;
    clear(recovered_partial);
    for (long u = 0; u < tensor.basis_dimension; ++u) {
      const vec_ZZ_pE s_u_coeffs =
          basefold::DecomposeGRElementToBaseCoeffsPolynomialBasis(params,
                                                                  s_by_u[u]);
      const ZZ_pE alpha_u = (u == 0) ? testutil::ConstZZpE(1) : x;
      recovered_partial += s_u_coeffs[v] * alpha_u;
    }

    vec_ZZ_pE slice;
    slice.SetLength(num_w);
    for (long w = 0; w < num_w; ++w) {
      slice[w] = t_coeffs[v + (w << params.kappa)];
    }
    const vec_ZZ_pE slice_monomial =
        basefold::BooleanHypercubeTableToMonomialCoeffs(slice);
    const ZZ_pE direct_partial =
        basefold::EvalMultilinearMonomialCoeffs(slice_monomial, r_suffix);
    CHECK_EQ(recovered_partial, direct_partial);
  }
}

void TestBuildRingSwitchComponentTensor_RCoeffsEvaluateAsExpected() {
  testutil::PrintInfo("Ring-switch WP1: verifier polynomial r evaluates with v||w flattening");

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

  const basefold::RingSwitchPCSParams params =
      BuildRingSwitchParamsGR42(/*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
  const std::vector<ZZ_pE> r_suffix = {alpha + testutil::ConstZZpE(1),
                                       testutil::ConstZZpE(3)};
  const basefold::RingSwitchComponentTensor tensor =
      basefold::BuildRingSwitchComponentTensor(params, r_suffix);

  const std::vector<ZZ_pE> r_prime = {alpha, testutil::ConstZZpE(1),
                                      alpha + testutil::ConstZZpE(1)};
  const ZZ_pE got =
      basefold::EvalMultilinearMonomialCoeffs(tensor.r_monomial_coeffs, r_prime);

  const std::vector<ZZ_pE> r_prefix = {r_prime[0]};
  const std::vector<ZZ_pE> r_suffix_prime = {r_prime[1], r_prime[2]};
  ZZ_pE expected;
  clear(expected);
  const long num_w = Pow2ForTest(params.ell_prime);
  for (long w = 0; w < num_w; ++w) {
    const ZZ_pE eq_w =
        basefold::EqPolynomial(r_suffix_prime, BooleanPointFromIndex(w, params.ell_prime));
    ZZ_pE inner;
    clear(inner);
    for (long u = 0; u < tensor.basis_dimension; ++u) {
      const ZZ_pE eq_u =
          basefold::EqPolynomial(r_prefix, BooleanPointFromIndex(u, params.kappa));
      inner += tensor.a_by_u_then_w[u * num_w + w] * eq_u;
    }
    expected += inner * eq_w;
  }
  CHECK_EQ(got, expected);
}

void TestBuildRingSwitchComponentTensor_RejectsWrongSuffixDimension() {
  testutil::PrintInfo("Ring-switch WP1: component tensor rejects wrong suffix dimension");

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

        const basefold::RingSwitchPCSParams params = BuildRingSwitchParamsGR42(
            /*ell=*/3, /*kappa=*/1, modulus, F, p, alpha);
        const std::vector<ZZ_pE> wrong_suffix = {alpha};
        (void)basefold::BuildRingSwitchComponentTensor(params, wrong_suffix);
      },
      "r_suffix dimension must equal ell_prime",
      "TestBuildRingSwitchComponentTensor_RejectsWrongSuffixDimension");
}

void TestProductSumcheckProver_BooleanTablesPasses() {
  testutil::PrintInfo("Ring-switch WP2: product sumcheck passes on honest Boolean tables");

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

  vec_ZZ_pE f_table;
  f_table.SetLength(4);
  f_table[0] = testutil::ConstZZpE(1);
  f_table[1] = alpha;
  f_table[2] = testutil::ConstZZpE(3);
  f_table[3] = alpha + testutil::ConstZZpE(1);

  vec_ZZ_pE g_table;
  g_table.SetLength(4);
  g_table[0] = alpha + testutil::ConstZZpE(2);
  g_table[1] = testutil::ConstZZpE(1);
  g_table[2] = alpha;
  g_table[3] = testutil::ConstZZpE(0);

  const ZZ_pE initial_claim = SumOfPointwiseProducts(f_table, g_table);
  basefold::ProductSumcheckProver prover(f_table, g_table);

  std::vector<basefold::QuadraticPoly> h_by_level(2);
  h_by_level[1] = prover.CurrentPolynomial();

  std::vector<ZZ_pE> r_by_level = {alpha + testutil::ConstZZpE(1),
                                   testutil::ConstZZpE(3)};
  prover.ReceiveChallenge(r_by_level[1]);
  h_by_level[0] = prover.CurrentPolynomial();
  prover.ReceiveChallenge(r_by_level[0]);

  CHECK(basefold::CheckProductSumcheckChain(initial_claim, h_by_level,
                                            r_by_level));

  const vec_ZZ_pE f_monomial =
      basefold::BooleanHypercubeTableToMonomialCoeffs(f_table);
  const vec_ZZ_pE g_monomial =
      basefold::BooleanHypercubeTableToMonomialCoeffs(g_table);
  const ZZ_pE expected_final =
      basefold::EvalMultilinearMonomialCoeffs(f_monomial, r_by_level) *
      basefold::EvalMultilinearMonomialCoeffs(g_monomial, r_by_level);
  CHECK_EQ(h_by_level[0].Eval(r_by_level[0]), expected_final);
}

void TestProductSumcheckChain_RejectsTamperedPolynomial() {
  testutil::PrintInfo("Ring-switch WP2: product sumcheck chain rejects tampering");

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

  vec_ZZ_pE f_table;
  f_table.SetLength(4);
  f_table[0] = testutil::ConstZZpE(0);
  f_table[1] = testutil::ConstZZpE(1);
  f_table[2] = alpha;
  f_table[3] = testutil::ConstZZpE(2);

  vec_ZZ_pE g_table;
  g_table.SetLength(4);
  g_table[0] = alpha + testutil::ConstZZpE(1);
  g_table[1] = testutil::ConstZZpE(3);
  g_table[2] = testutil::ConstZZpE(1);
  g_table[3] = alpha;

  const ZZ_pE initial_claim = SumOfPointwiseProducts(f_table, g_table);
  basefold::ProductSumcheckProver prover(f_table, g_table);

  std::vector<basefold::QuadraticPoly> h_by_level(2);
  h_by_level[1] = prover.CurrentPolynomial();

  const std::vector<ZZ_pE> r_by_level = {alpha, testutil::ConstZZpE(1)};
  prover.ReceiveChallenge(r_by_level[1]);
  h_by_level[0] = prover.CurrentPolynomial();

  CHECK(basefold::CheckProductSumcheckChain(initial_claim, h_by_level,
                                            r_by_level));

  std::vector<basefold::QuadraticPoly> tampered = h_by_level;
  tampered[0].a0 += testutil::ConstZZpE(1);
  CHECK(!basefold::CheckProductSumcheckChain(initial_claim, tampered,
                                             r_by_level));
}

void TestProductSumcheckProver_FromMonomialCoeffsMatchesTables() {
  testutil::PrintInfo("Ring-switch WP2: monomial helper matches table-native product sumcheck");

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

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(4);
  f_coeffs[0] = testutil::ConstZZpE(1);
  f_coeffs[1] = alpha;
  f_coeffs[2] = testutil::ConstZZpE(2);
  f_coeffs[3] = alpha + testutil::ConstZZpE(1);

  vec_ZZ_pE g_coeffs;
  g_coeffs.SetLength(4);
  g_coeffs[0] = alpha + testutil::ConstZZpE(1);
  g_coeffs[1] = testutil::ConstZZpE(3);
  g_coeffs[2] = testutil::ConstZZpE(0);
  g_coeffs[3] = testutil::ConstZZpE(1);

  vec_ZZ_pE f_table;
  vec_ZZ_pE g_table;
  f_table.SetLength(4);
  g_table.SetLength(4);
  for (long i = 0; i < 4; ++i) {
    const std::vector<ZZ_pE> point = BooleanPointFromIndex(i, 2);
    f_table[i] = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, point);
    g_table[i] = basefold::EvalMultilinearMonomialCoeffs(g_coeffs, point);
  }

  basefold::ProductSumcheckProver from_tables(f_table, g_table);
  basefold::ProductSumcheckProver from_coeffs =
      basefold::ProductSumcheckProver::FromMonomialCoeffs(f_coeffs, g_coeffs);

  CHECK_EQ(from_tables.CurrentPolynomial().a0, from_coeffs.CurrentPolynomial().a0);
  CHECK_EQ(from_tables.CurrentPolynomial().a1, from_coeffs.CurrentPolynomial().a1);
  CHECK_EQ(from_tables.CurrentPolynomial().a2, from_coeffs.CurrentPolynomial().a2);

  const std::vector<ZZ_pE> r_by_level = {alpha + testutil::ConstZZpE(1),
                                         testutil::ConstZZpE(2)};
  from_tables.ReceiveChallenge(r_by_level[1]);
  from_coeffs.ReceiveChallenge(r_by_level[1]);

  CHECK_EQ(from_tables.CurrentPolynomial().a0, from_coeffs.CurrentPolynomial().a0);
  CHECK_EQ(from_tables.CurrentPolynomial().a1, from_coeffs.CurrentPolynomial().a1);
  CHECK_EQ(from_tables.CurrentPolynomial().a2, from_coeffs.CurrentPolynomial().a2);

  from_tables.ReceiveChallenge(r_by_level[0]);
  from_coeffs.ReceiveChallenge(r_by_level[0]);
}

void TestProductSumcheckProver_DimensionZeroUsesNoRounds() {
  testutil::PrintInfo("Ring-switch WP2: product sumcheck has an explicit d=0 no-round contract");

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

  vec_ZZ_pE f_table;
  f_table.SetLength(1);
  f_table[0] = alpha + testutil::ConstZZpE(1);

  vec_ZZ_pE g_table;
  g_table.SetLength(1);
  g_table[0] = alpha;

  basefold::ProductSumcheckProver prover(f_table, g_table);
  CHECK_EQ(prover.Dimension(), 0L);
  CHECK_EQ(prover.RemainingVars(), 0L);

  std::vector<basefold::QuadraticPoly> h_by_level;
  std::vector<ZZ_pE> r_by_level;
  CHECK(basefold::CheckProductSumcheckChain(
      f_table[0] * g_table[0], h_by_level, r_by_level));
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
    RUN_TEST(TestPackZ2kCoeffsToGREvals_RoundTripsSmallExample);
    RUN_TEST(TestBooleanHypercubeTableToMonomialCoeffs_EvaluatesTableCorrectly);
    RUN_TEST(TestPackZ2kCoeffsToGREvals_RejectsNonBaseRingInputs);
    RUN_TEST(TestBuildRingSwitchComponentTensor_ReconstructsSuffixEqualityValues);
    RUN_TEST(TestBuildRingSwitchComponentTensor_RecoversPartialEvaluations);
    RUN_TEST(TestBuildRingSwitchComponentTensor_RCoeffsEvaluateAsExpected);
    RUN_TEST(TestBuildRingSwitchComponentTensor_RejectsWrongSuffixDimension);
    RUN_TEST(TestProductSumcheckProver_BooleanTablesPasses);
    RUN_TEST(TestProductSumcheckChain_RejectsTamperedPolynomial);
    RUN_TEST(TestProductSumcheckProver_FromMonomialCoeffsMatchesTables);
    RUN_TEST(TestProductSumcheckProver_DimensionZeroUsesNoRounds);
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
