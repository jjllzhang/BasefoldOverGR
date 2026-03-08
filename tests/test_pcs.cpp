#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "BaseFold/BaseFoldPCS.hpp"
#include "BaseFold/Multilinear.hpp"
#include "BaseFold/ProofSize.hpp"
#include "tests/test_common.hpp"

using NTL::conv;
using NTL::mat_ZZ_pE;
using NTL::SetCoeff;
using NTL::to_ZZ;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZ_pE;
using NTL::ZZ_pEX;
using NTL::ZZ_pEPush;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using std::cerr;
using std::cout;
using std::exception;
using std::ostringstream;
using std::string;

int g_failures = 0;

namespace {

basefold::FoldableCodeParams BuildParamsGF4_k0_1(const ZZ &p, const ZZ_pE &alpha) {
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

  const ZZ_pE one = testutil::ConstZZpE(1);

  mat_ZZ_pE G0;
  G0.SetDims(k0, n0);
  G0[0][0] = one;
  G0[0][1] = alpha;
  params.G0 = G0;

  params.diag_T.resize(static_cast<std::size_t>(d));
  params.diag_T[0].SetLength(/*n0=*/2);
  params.diag_T[0][0] = one;
  params.diag_T[0][1] = one;

  params.diag_T[1].SetLength(/*n1=*/4);
  for (long i = 0; i < 4; ++i) {
    params.diag_T[1][i] = one;
  }

  return params;
}

basefold::FoldableCodeParams BuildParamsGF4_k0_2(const ZZ &p, const ZZ_pE &alpha) {
  const long c = 2;
  const long k0 = 2;
  const long d = 2;
  const long n0 = c * k0;

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = p;
  params.zeta = alpha;

  const ZZ_pE one = testutil::ConstZZpE(1);

  mat_ZZ_pE G0;
  G0.SetDims(k0, n0);
  G0[0][0] = one;
  G0[0][1] = ZZ_pE(0);
  G0[0][2] = one;
  G0[0][3] = ZZ_pE(0);

  G0[1][0] = ZZ_pE(0);
  G0[1][1] = one;
  G0[1][2] = ZZ_pE(0);
  G0[1][3] = one;
  params.G0 = G0;

  params.diag_T.resize(static_cast<std::size_t>(d));
  for (long level = 0; level < d; ++level) {
    const long ni = c * k0 * (1L << level);
    params.diag_T[static_cast<std::size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      params.diag_T[static_cast<std::size_t>(level)][i] = one;
    }
  }

  return params;
}

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

  params.diag_T[0].SetLength(/*n0=*/2);
  params.diag_T[0][0] = one;
  params.diag_T[0][1] = alpha;

  params.diag_T[1].SetLength(/*n1=*/4);
  params.diag_T[1][0] = alpha + one;
  params.diag_T[1][1] = one;
  params.diag_T[1][2] = alpha;
  params.diag_T[1][3] = alpha + one;

  return params;
}

basefold::FoldableCodeParams BuildParamsGR42_k0_2(const ZZ &p, const ZZ_pE &alpha) {
  const long c = 2;
  const long k0 = 2;
  const long d = 2;
  const long n0 = c * k0;

  basefold::FoldableCodeParams params;
  params.c = c;
  params.k0 = k0;
  params.d = d;
  params.p = p;
  params.zeta = alpha;

  const ZZ_pE one = testutil::ConstZZpE(1);

  mat_ZZ_pE G0;
  G0.SetDims(k0, n0);
  G0[0][0] = one;
  G0[0][1] = ZZ_pE(0);
  G0[0][2] = one;
  G0[0][3] = ZZ_pE(0);

  G0[1][0] = ZZ_pE(0);
  G0[1][1] = one;
  G0[1][2] = ZZ_pE(0);
  G0[1][3] = one;
  params.G0 = G0;

  params.diag_T.resize(static_cast<std::size_t>(d));
  for (long level = 0; level < d; ++level) {
    const long ni = c * k0 * (1L << level);
    params.diag_T[static_cast<std::size_t>(level)].SetLength(ni);
    for (long i = 0; i < ni; ++i) {
      params.diag_T[static_cast<std::size_t>(level)][i] = one;
    }
  }

  return params;
}

void TestPCS_EvalProof_GF4() {
  testutil::PrintInfo("PCS: eval proof verifies over GF(2^2)");

  const ZZ p = to_ZZ(2);
  ZZ_pPush p_push(p);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FoldableCodeParams params = BuildParamsGF4_k0_1(p, alpha);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  f_coeffs[0] = testutil::ConstZZpE(0);
  f_coeffs[1] = testutil::ConstZZpE(1);
  f_coeffs[2] = alpha;
  f_coeffs[3] = alpha + testutil::ConstZZpE(1);

  const std::vector<ZZ_pE> z = {alpha, alpha + testutil::ConstZZpE(1)};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  const basefold::MerkleRoot C = basefold::BaseFoldPCSCommit(f_coeffs, params);
  const basefold::BaseFoldPCSCommitArtifacts commit_artifacts =
      basefold::BaseFoldPCSBuildCommitArtifacts(f_coeffs, params);
  CHECK(commit_artifacts.root_d == C);
  const long num_queries = 3;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);
  const basefold::BaseFoldPCSEvalProof proof_from_committed =
      basefold::BaseFoldPCSProveEvalFromCommittedTopOracle(
          f_coeffs, z, y, num_queries, params, commit_artifacts);

  CHECK(static_cast<long>(proof.query_multiproofs.size()) == params.d + 1);
  CHECK(basefold::BaseFoldPCSVerifyEval(C, z, y, num_queries, proof, params));
  CHECK(basefold::BaseFoldPCSVerifyEval(C, z, y, num_queries,
                                        proof_from_committed, params));

  const ZZ_pE y_bad = y + testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(C, z, y_bad, num_queries, proof, params));

  basefold::BaseFoldPCSEvalProof proof_tampered = proof;
  CHECK(proof_tampered.query_multiproofs[0].values.length() > 0);
  proof_tampered.query_multiproofs[0].values[0] += testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(C, z, y, num_queries, proof_tampered,
                                         params));
}

void TestPCS_EvalProof_GF4_k0_2() {
  testutil::PrintInfo("PCS: eval proof verifies over GF(2^2) (k0=2)");

  const ZZ p = to_ZZ(2);
  ZZ_pPush p_push(p);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FoldableCodeParams params = BuildParamsGF4_k0_2(p, alpha);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  for (long i = 0; i < f_coeffs.length(); ++i) {
    f_coeffs[i] = testutil::ConstZZpE(i % 2);
  }

  const std::vector<ZZ_pE> z = {alpha, alpha + testutil::ConstZZpE(1), alpha};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  const basefold::MerkleRoot C = basefold::BaseFoldPCSCommit(f_coeffs, params);
  const long num_queries = 3;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);

  CHECK(basefold::BaseFoldPCSVerifyEval(C, z, y, num_queries, proof, params));

  const ZZ_pE y_bad = y + testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(C, z, y_bad, num_queries, proof, params));
}

void TestPCS_EvalProof_GR42() {
  testutil::PrintInfo("PCS: eval proof verifies over GR(4,2)");

  const ZZ p = to_ZZ(2);
  const ZZ mod = ZZ(4);
  ZZ_pPush mod_push(mod);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FoldableCodeParams params = BuildParamsGR42(p, alpha);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  f_coeffs[0] = testutil::ConstZZpE(0);
  f_coeffs[1] = testutil::ConstZZpE(1);
  f_coeffs[2] = alpha;
  f_coeffs[3] = testutil::ConstZZpE(2);

  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), testutil::ConstZZpE(3)};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  const basefold::MerkleRoot C = basefold::BaseFoldPCSCommit(f_coeffs, params);
  const long num_queries = 4;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);

  CHECK(static_cast<long>(proof.query_multiproofs.size()) == params.d + 1);
  CHECK(basefold::BaseFoldPCSVerifyEval(C, z, y, num_queries, proof, params));

  const ZZ_pE y_bad = y + testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(C, z, y_bad, num_queries, proof, params));
}

void TestPCS_EvalProof_GR42_k0_2() {
  testutil::PrintInfo("PCS: eval proof verifies over GR(4,2) (k0=2)");

  const ZZ p = to_ZZ(2);
  const ZZ mod = ZZ(4);
  ZZ_pPush mod_push(mod);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FoldableCodeParams params = BuildParamsGR42_k0_2(p, alpha);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  for (long i = 0; i < f_coeffs.length(); ++i) {
    f_coeffs[i] = testutil::ConstZZpE(i % 4);
  }

  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), testutil::ConstZZpE(3), alpha};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  const basefold::MerkleRoot C = basefold::BaseFoldPCSCommit(f_coeffs, params);
  const long num_queries = 4;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);

  CHECK(basefold::BaseFoldPCSVerifyEval(C, z, y, num_queries, proof, params));

  const ZZ_pE y_bad = y + testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(C, z, y_bad, num_queries, proof, params));
}

void TestPCS_EvalProof_ExtChallengeConfig_GF4() {
  testutil::PrintInfo("PCS: challenge-config path verifies over GF(2^2)");

  const ZZ p = to_ZZ(2);
  ZZ_pPush p_push(p);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FoldableCodeParams params = BuildParamsGF4_k0_1(p, alpha);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  f_coeffs[0] = testutil::ConstZZpE(0);
  f_coeffs[1] = testutil::ConstZZpE(1);
  f_coeffs[2] = alpha;
  f_coeffs[3] = alpha + testutil::ConstZZpE(1);

  const std::vector<ZZ_pE> z = {alpha, alpha + testutil::ConstZZpE(1)};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
  const basefold::MerkleRoot C = basefold::BaseFoldPCSCommit(f_coeffs, params);
  const basefold::BaseFoldPCSCommitArtifacts commit_artifacts =
      basefold::BaseFoldPCSBuildCommitArtifacts(f_coeffs, params);

  basefold::BaseFoldPCSChallengeConfig cfg;
  cfg.use_extension_challenges = true;
  ZZ_pEX E;
  SetCoeff(E, 0, alpha);
  SetCoeff(E, 1, testutil::ConstZZpE(1));
  SetCoeff(E, 2, testutil::ConstZZpE(1));
  cfg.challenge_extension_modulus = E;

  const long num_queries = 3;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEvalWithChallengeConfig(f_coeffs, z, y,
                                                        num_queries, params, cfg);
  const basefold::BaseFoldPCSEvalProof proof_from_committed =
      basefold::BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle(
          f_coeffs, z, y, num_queries, params, commit_artifacts, cfg);

  CHECK(proof.extension.enabled);
  CHECK(proof.extension.r_by_level.empty());
  CHECK(static_cast<long>(proof.extension.roots_by_level.size()) == params.d);
  CHECK(static_cast<long>(proof.extension.query_multiproofs.size()) == params.d);
  CHECK(!proof.extension.base_top_query_multiproof.queried_indices.empty());
  CHECK(proof.extension.base_top_query_multiproof.values.length() > 0);

  CHECK(basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof, params, cfg));
  CHECK(basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_from_committed, params, cfg));

  CHECK(!basefold::BaseFoldPCSVerifyEval(C, z, y, num_queries, proof, params));

  basefold::BaseFoldPCSEvalProof proof_ext_tampered = proof;
  CHECK(!proof_ext_tampered.extension.query_multiproofs.empty());
  CHECK(!proof_ext_tampered.extension.query_multiproofs[0].values.empty());
  ZZ_pE coeff1 = NTL::coeff(
      proof_ext_tampered.extension.query_multiproofs[0].values[0], 1);
  coeff1 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_ext_tampered.extension.query_multiproofs[0].values[0], 1,
                coeff1);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_ext_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_top_tampered = proof;
  CHECK(proof_top_tampered.extension.base_top_query_multiproof.values.length() > 0);
  proof_top_tampered.extension.base_top_query_multiproof.values[0] +=
      testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_top_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_root_tampered = proof;
  proof_root_tampered.extension.roots_by_level[0][0] ^= static_cast<basefold::Byte>(0x01);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_root_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_h_tampered = proof;
  ZZ_pE h0_a0_coeff0 =
      NTL::coeff(proof_h_tampered.extension.h_by_level[0].a0, 0);
  h0_a0_coeff0 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_h_tampered.extension.h_by_level[0].a0, 0, h0_a0_coeff0);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_h_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_msg0_tampered = proof;
  ZZ_pE msg0_coeff0 = NTL::coeff(proof_msg0_tampered.extension.msg0_coeffs[0], 0);
  msg0_coeff0 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_msg0_tampered.extension.msg0_coeffs[0], 0, msg0_coeff0);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_msg0_tampered, params, cfg));

  basefold::BaseFoldPCSChallengeConfig cfg_bad = cfg;
  ZZ_pEX E_bad;
  SetCoeff(E_bad, 0, alpha + testutil::ConstZZpE(1));
  SetCoeff(E_bad, 1, testutil::ConstZZpE(1));
  SetCoeff(E_bad, 2, testutil::ConstZZpE(1));
  cfg_bad.challenge_extension_modulus = E_bad;
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof, params, cfg_bad));

}

void TestPCS_EvalProof_ExtChallengeConfig_GR42() {
  testutil::PrintInfo("PCS: challenge-config path verifies over GR(4,2)");

  const ZZ p = to_ZZ(2);
  const ZZ mod = ZZ(4);
  ZZ_pPush mod_push(mod);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  ZZ_pX xpoly;
  SetCoeff(xpoly, 1, 1);
  ZZ_pE alpha;
  conv(alpha, xpoly);

  const basefold::FoldableCodeParams params = BuildParamsGR42(p, alpha);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  f_coeffs[0] = testutil::ConstZZpE(0);
  f_coeffs[1] = testutil::ConstZZpE(1);
  f_coeffs[2] = alpha;
  f_coeffs[3] = testutil::ConstZZpE(2);

  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                testutil::ConstZZpE(3)};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
  const basefold::MerkleRoot C = basefold::BaseFoldPCSCommit(f_coeffs, params);

  basefold::BaseFoldPCSChallengeConfig cfg;
  cfg.use_extension_challenges = true;
  ZZ_pEX E;
  SetCoeff(E, 0, alpha);
  SetCoeff(E, 1, testutil::ConstZZpE(1));
  SetCoeff(E, 2, testutil::ConstZZpE(1));
  cfg.challenge_extension_modulus = E;

  const long num_queries = 4;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEvalWithChallengeConfig(
          f_coeffs, z, y, num_queries, params, cfg);

  CHECK(proof.extension.enabled);
  CHECK(proof.extension.r_by_level.empty());
  CHECK(static_cast<long>(proof.extension.roots_by_level.size()) == params.d);
  CHECK(static_cast<long>(proof.extension.query_multiproofs.size()) == params.d);
  CHECK(!proof.extension.base_top_query_multiproof.queried_indices.empty());
  CHECK(proof.extension.base_top_query_multiproof.values.length() > 0);
  CHECK(basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_open_tampered = proof;
  CHECK(!proof_open_tampered.extension.query_multiproofs.empty());
  CHECK(!proof_open_tampered.extension.query_multiproofs[0].values.empty());
  ZZ_pE coeff1 = NTL::coeff(
      proof_open_tampered.extension.query_multiproofs[0].values[0], 1);
  coeff1 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_open_tampered.extension.query_multiproofs[0].values[0], 1,
                coeff1);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_open_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_top_tampered = proof;
  CHECK(proof_top_tampered.extension.base_top_query_multiproof.values.length() > 0);
  proof_top_tampered.extension.base_top_query_multiproof.values[0] +=
      testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_top_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_root_tampered = proof;
  proof_root_tampered.extension.roots_by_level[0][0] ^=
      static_cast<basefold::Byte>(0x01);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_root_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_h_tampered = proof;
  ZZ_pE h0_a0_coeff0 =
      NTL::coeff(proof_h_tampered.extension.h_by_level[0].a0, 0);
  h0_a0_coeff0 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_h_tampered.extension.h_by_level[0].a0, 0, h0_a0_coeff0);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_h_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_msg0_tampered = proof;
  ZZ_pE msg0_coeff0 = NTL::coeff(proof_msg0_tampered.extension.msg0_coeffs[0], 0);
  msg0_coeff0 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_msg0_tampered.extension.msg0_coeffs[0], 0, msg0_coeff0);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      C, z, y, num_queries, proof_msg0_tampered, params, cfg));

}

void TestPCS_ProofSizeFixedWidth_GF4_HandCheck() {
  testutil::PrintInfo("PCS: fixed-width proof size matches hand calculation");

  const ZZ p = to_ZZ(2);
  ZZ_pPush p_push(p);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  basefold::BaseFoldPCSEvalProof proof;
  proof.commitments.roots_by_level.resize(1);
  proof.h_by_level.resize(1);
  proof.pi0_full.SetLength(2);
  proof.query_multiproofs.resize(1);
  proof.query_multiproofs[0].queried_indices = {0, 1};
  proof.query_multiproofs[0].values.SetLength(2);
  proof.query_multiproofs[0].sibling_hashes.resize(3);
  proof.extension.enabled = false;

  const std::uint64_t expected_bytes = 192;
  const std::uint64_t actual_bytes =
      basefold::BaseFoldPCSEvalProofSizeBytes(proof);
  CHECK_EQ(actual_bytes, expected_bytes);
  CHECK_EQ(basefold::BaseFoldPCSEvalProofSizeKB(proof),
           static_cast<double>(expected_bytes) / 1024.0);

  basefold::BaseFoldPCSEvalProof proof_without_indices = proof;
  proof_without_indices.query_multiproofs[0].queried_indices.clear();
  CHECK_EQ(basefold::BaseFoldPCSEvalProofSizeBytes(proof_without_indices),
           expected_bytes);
}

void TestPCS_ProofSizeFixedWidth_ExtensionWidthDerivation() {
  testutil::PrintInfo("PCS: fixed-width proof size derives extension width");

  const ZZ p = to_ZZ(2);
  ZZ_pPush p_push(p);

  ZZ_pX F;
  SetCoeff(F, 2, 1);
  SetCoeff(F, 1, 1);
  SetCoeff(F, 0, 1);
  ZZ_pEPush e_push(F);

  basefold::BaseFoldPCSEvalProof proof;
  proof.extension.enabled = true;
  proof.extension.msg0_coeffs.resize(1);
  proof.extension.r_by_level.resize(2);

  basefold::BaseFoldProofSizeOptions degree2;
  degree2.challenge_ext_degree = 2;
  CHECK_EQ(basefold::BaseFoldPCSEvalProofSizeBytes(proof, degree2),
           static_cast<std::uint64_t>(94));
  CHECK_EQ(basefold::BaseFoldPCSEvalProofSizeKB(proof, degree2),
           94.0 / 1024.0);

  basefold::BaseFoldProofSizeOptions degree3;
  degree3.challenge_ext_degree = 3;
  CHECK_EQ(basefold::BaseFoldPCSEvalProofSizeBytes(proof, degree3),
           static_cast<std::uint64_t>(96));

  degree3.include_version_byte = false;
  CHECK_EQ(basefold::BaseFoldPCSEvalProofSizeBytes(proof, degree3),
           static_cast<std::uint64_t>(95));
}

}  // namespace

int main() {
  try {
    RUN_TEST(TestPCS_EvalProof_GF4);
    RUN_TEST(TestPCS_EvalProof_GF4_k0_2);
    RUN_TEST(TestPCS_EvalProof_GR42);
    RUN_TEST(TestPCS_EvalProof_GR42_k0_2);
    RUN_TEST(TestPCS_EvalProof_ExtChallengeConfig_GF4);
    RUN_TEST(TestPCS_EvalProof_ExtChallengeConfig_GR42);
    RUN_TEST(TestPCS_ProofSizeFixedWidth_GF4_HandCheck);
    RUN_TEST(TestPCS_ProofSizeFixedWidth_ExtensionWidthDerivation);
  } catch (const exception &e) {
    cerr << "Unhandled std::exception: " << e.what() << "\n";
    return 2;
  } catch (...) {
    cerr << "Unhandled non-std exception\n";
    return 2;
  }

  if (g_failures == 0) {
    cout << "\nAll tests passed.\n";
    return 0;
  }

  cerr << "\n" << g_failures << " test(s) failed.\n";
  return 1;
}
