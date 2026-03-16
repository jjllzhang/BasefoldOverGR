#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <exception>
#include <functional>
#include <iostream>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

#include "bench/bench_basefold_pcs_artifact_common.hpp"
#include "PCS/BaseFold/BaseFoldPCS.hpp"
#include "PCS/BaseFold/ProofDeserialize.hpp"
#include "PCS/BaseFold/ProofSerialize.hpp"
#include "PCS/Common/MerkleMultiproofReplay.hpp"
#include "PCS/Common/Multilinear.hpp"
#include "PCS/BaseFold/ProofSize.hpp"
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
namespace fs = std::filesystem;

int g_test_failure_count = 0;

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

basefold::FoldableCodeParams BuildParamsGR_4_2(const ZZ &p,
                                               const ZZ_pE &alpha) {
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

basefold::FoldableCodeParams BuildParamsGR_4_2_k0_2(const ZZ &p,
                                                    const ZZ_pE &alpha) {
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

std::vector<ZZ_pE> BooleanPointFromIndex(long index, long dimension) {
  std::vector<ZZ_pE> point(static_cast<std::size_t>(dimension));
  for (long i = 0; i < dimension; ++i) {
    point[static_cast<std::size_t>(i)] =
        testutil::ConstZZpE((index >> i) & 1L);
  }
  return point;
}

void TestPCS_EqualityTableFromPoint_MatchesEqPolynomialDefinition() {
  testutil::PrintInfo(
      "PCS common: EqualityTableFromPoint matches the EqPolynomial definition");

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

  const std::vector<ZZ_pE> point = {
      alpha, alpha + testutil::ConstZZpE(1), testutil::ConstZZpE(1)};
  const vec_ZZ_pE table = basefold::EqualityTableFromPoint(point);
  CHECK_EQ(table.length(), 8);
  for (long idx = 0; idx < table.length(); ++idx) {
    CHECK_EQ(table[idx],
             basefold::EqPolynomial(point, BooleanPointFromIndex(idx, 3)));
  }

  const vec_ZZ_pE empty_table =
      basefold::EqualityTableFromPoint(std::vector<ZZ_pE>{});
  CHECK_EQ(empty_table.length(), 1);
  CHECK_EQ(empty_table[0], testutil::ConstZZpE(1));
}

void TestPCS_MultiproofValueCache_ReplaysExpectedPositions() {
  testutil::PrintInfo(
      "PCS common: multiproof value cache returns the same positions as queried indices");

  const std::vector<long> queried_indices = {2, 5, 9};
  std::vector<long> values = {20, 50, 90};
  basefold::multiproof_replay::ValuePositionCache cache;
  CHECK(basefold::multiproof_replay::BuildValuePositionCache(
      queried_indices, static_cast<long>(values.size()), &cache));

  const long *found = basefold::multiproof_replay::FindMultiproofValue(
      cache, 5, [&](std::size_t pos) { return &values[pos]; });
  CHECK(found != nullptr);
  CHECK_EQ(*found, 50);

  const long *missing = basefold::multiproof_replay::FindMultiproofValue(
      cache, 6, [&](std::size_t pos) { return &values[pos]; });
  CHECK(missing == nullptr);

  CHECK(!basefold::multiproof_replay::BuildValuePositionCache(
      queried_indices, /*value_count=*/2, &cache));
  CHECK(!basefold::multiproof_replay::BuildValuePositionCache(
      std::vector<long>{2, 2, 5}, /*value_count=*/3, &cache));
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

  const basefold::MerkleRoot commitment_root =
      basefold::BaseFoldPCSCommit(f_coeffs, params);
  const basefold::BaseFoldPCSCommitArtifacts commit_artifacts =
      basefold::BaseFoldPCSBuildCommitArtifacts(f_coeffs, params);
  CHECK(commit_artifacts.root_d == commitment_root);
  CHECK(commit_artifacts.base_sumcheck_precomputation.valid);
  CHECK_EQ(commit_artifacts.base_sumcheck_precomputation.f_eval_table.length(),
           f_coeffs.length());
  const long num_queries = 3;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);
  const basefold::BaseFoldPCSEvalProof proof_from_committed =
      basefold::BaseFoldPCSProveEvalFromCommittedTopOracle(
          f_coeffs, z, y, num_queries, params, commit_artifacts);
  basefold::BaseFoldPCSCommitArtifacts commit_artifacts_without_sumcheck =
      commit_artifacts;
  commit_artifacts_without_sumcheck.base_sumcheck_precomputation = {};
  const basefold::BaseFoldPCSEvalProof proof_from_committed_without_sumcheck =
      basefold::BaseFoldPCSProveEvalFromCommittedTopOracle(
          f_coeffs, z, y, num_queries, params,
          commit_artifacts_without_sumcheck);

  CHECK(static_cast<long>(proof.query_multiproofs.size()) == params.d + 1);
  CHECK(basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                        proof, params));
  CHECK(basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                        proof_from_committed, params));
  CHECK(basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                        proof_from_committed_without_sumcheck,
                                        params));
  basefold::FixedProofEncodingOptions encoding_options;
  encoding_options.include_version_byte = true;
  CHECK_EQ(basefold::SerializeBaseFoldPCSEvalProofFixedBytes(
               proof_from_committed, encoding_options),
           basefold::SerializeBaseFoldPCSEvalProofFixedBytes(
               proof_from_committed_without_sumcheck, encoding_options));
  basefold::BaseFoldPCSEvalProof proof_without_indices = proof;
  for (basefold::MerkleMultiproof &multiproof :
       proof_without_indices.query_multiproofs) {
    multiproof.queried_indices.clear();
  }
  CHECK(basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                        proof_without_indices, params));

  const ZZ_pE y_bad = y + testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(commitment_root, z, y_bad, num_queries,
                                         proof, params));

  basefold::BaseFoldPCSEvalProof proof_tampered = proof;
  CHECK(proof_tampered.query_multiproofs[0].values.length() > 0);
  proof_tampered.query_multiproofs[0].values[0] += testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                         proof_tampered, params));
}

void TestPCS_PaperAPI_GF4() {
  testutil::PrintInfo("PCS: staged setup/commit/prove/verify API matches the legacy flow");

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
  const basefold::BaseFoldPCSSetupOutput api = basefold::BaseFoldPCSSetup(params);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  f_coeffs[0] = testutil::ConstZZpE(0);
  f_coeffs[1] = testutil::ConstZZpE(1);
  f_coeffs[2] = alpha;
  f_coeffs[3] = alpha + testutil::ConstZZpE(1);

  const std::vector<ZZ_pE> z = {alpha, alpha + testutil::ConstZZpE(1)};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
  const long num_queries = 3;

  const basefold::BaseFoldPCSCommittedWitness committed =
      api.prover.Commit(f_coeffs);
  CHECK_EQ(committed.commitment, basefold::BaseFoldPCSCommit(f_coeffs, params));
  CHECK_EQ(committed.commitment, committed.commit_artifacts.root_d);

  const basefold::BaseFoldPCSEvalProof proof =
      api.prover.Prove(committed, z, y, num_queries);
  CHECK(api.verifier.Verify(committed.commitment, z, y, num_queries, proof));
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

  const basefold::MerkleRoot commitment_root =
      basefold::BaseFoldPCSCommit(f_coeffs, params);
  const long num_queries = 3;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);

  CHECK(basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                        proof, params));

  const ZZ_pE y_bad = y + testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(commitment_root, z, y_bad, num_queries,
                                         proof, params));
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

  const basefold::FoldableCodeParams params = BuildParamsGR_4_2(p, alpha);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  f_coeffs[0] = testutil::ConstZZpE(0);
  f_coeffs[1] = testutil::ConstZZpE(1);
  f_coeffs[2] = alpha;
  f_coeffs[3] = testutil::ConstZZpE(2);

  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), testutil::ConstZZpE(3)};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  const basefold::MerkleRoot commitment_root =
      basefold::BaseFoldPCSCommit(f_coeffs, params);
  const long num_queries = 4;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);

  CHECK(static_cast<long>(proof.query_multiproofs.size()) == params.d + 1);
  CHECK(basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                        proof, params));

  const ZZ_pE y_bad = y + testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(commitment_root, z, y_bad, num_queries,
                                         proof, params));
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

  const basefold::FoldableCodeParams params = BuildParamsGR_4_2_k0_2(p, alpha);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  for (long i = 0; i < f_coeffs.length(); ++i) {
    f_coeffs[i] = testutil::ConstZZpE(i % 4);
  }

  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1), testutil::ConstZZpE(3), alpha};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  const basefold::MerkleRoot commitment_root =
      basefold::BaseFoldPCSCommit(f_coeffs, params);
  const long num_queries = 4;
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);

  CHECK(basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                        proof, params));

  const ZZ_pE y_bad = y + testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEval(commitment_root, z, y_bad, num_queries,
                                         proof, params));
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
  const basefold::MerkleRoot commitment_root =
      basefold::BaseFoldPCSCommit(f_coeffs, params);
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

  CHECK(proof.extension.has_extension_payload);
  CHECK(proof.extension.r_by_level.empty());
  CHECK(static_cast<long>(proof.extension.roots_by_level.size()) == params.d);
  CHECK(static_cast<long>(proof.extension.query_multiproofs.size()) == params.d);
  CHECK(!proof.extension.base_top_query_multiproof.queried_indices.empty());
  CHECK(proof.extension.base_top_query_multiproof.values.length() > 0);

  CHECK(basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof, params, cfg));
  CHECK(basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_from_committed, params, cfg));
  basefold::BaseFoldPCSEvalProof proof_without_metadata = proof;
  proof_without_metadata.extension.r_by_level.resize(
      static_cast<std::size_t>(params.d));
  proof_without_metadata.extension.base_top_query_multiproof.queried_indices
      .clear();
  for (basefold::ExtensionMerkleMultiproof &multiproof :
       proof_without_metadata.extension.query_multiproofs) {
    multiproof.queried_indices.clear();
  }
  CHECK(basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_without_metadata, params, cfg));

  CHECK(!basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                         proof, params));

  basefold::BaseFoldPCSEvalProof proof_ext_tampered = proof;
  CHECK(!proof_ext_tampered.extension.query_multiproofs.empty());
  CHECK(!proof_ext_tampered.extension.query_multiproofs[0].values.empty());
  ZZ_pE coeff1 = NTL::coeff(
      proof_ext_tampered.extension.query_multiproofs[0].values[0], 1);
  coeff1 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_ext_tampered.extension.query_multiproofs[0].values[0], 1,
                coeff1);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_ext_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_top_tampered = proof;
  CHECK(proof_top_tampered.extension.base_top_query_multiproof.values.length() > 0);
  proof_top_tampered.extension.base_top_query_multiproof.values[0] +=
      testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_top_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_root_tampered = proof;
  proof_root_tampered.extension.roots_by_level[0][0] ^= static_cast<basefold::Byte>(0x01);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_root_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_h_tampered = proof;
  ZZ_pE h0_a0_coeff0 =
      NTL::coeff(proof_h_tampered.extension.h_by_level[0].a0, 0);
  h0_a0_coeff0 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_h_tampered.extension.h_by_level[0].a0, 0, h0_a0_coeff0);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_h_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_msg0_tampered = proof;
  ZZ_pE msg0_coeff0 = NTL::coeff(proof_msg0_tampered.extension.msg0_coeffs[0], 0);
  msg0_coeff0 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_msg0_tampered.extension.msg0_coeffs[0], 0, msg0_coeff0);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_msg0_tampered, params, cfg));

  basefold::BaseFoldPCSChallengeConfig cfg_bad = cfg;
  ZZ_pEX E_bad;
  SetCoeff(E_bad, 0, alpha + testutil::ConstZZpE(1));
  SetCoeff(E_bad, 1, testutil::ConstZZpE(1));
  SetCoeff(E_bad, 2, testutil::ConstZZpE(1));
  cfg_bad.challenge_extension_modulus = E_bad;
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof, params, cfg_bad));

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

  const basefold::FoldableCodeParams params = BuildParamsGR_4_2(p, alpha);

  vec_ZZ_pE f_coeffs;
  f_coeffs.SetLength(basefold::MessageLength(params));
  f_coeffs[0] = testutil::ConstZZpE(0);
  f_coeffs[1] = testutil::ConstZZpE(1);
  f_coeffs[2] = alpha;
  f_coeffs[3] = testutil::ConstZZpE(2);

  const std::vector<ZZ_pE> z = {alpha + testutil::ConstZZpE(1),
                                testutil::ConstZZpE(3)};
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);
  const basefold::MerkleRoot commitment_root =
      basefold::BaseFoldPCSCommit(f_coeffs, params);

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

  CHECK(proof.extension.has_extension_payload);
  CHECK(proof.extension.r_by_level.empty());
  CHECK(static_cast<long>(proof.extension.roots_by_level.size()) == params.d);
  CHECK(static_cast<long>(proof.extension.query_multiproofs.size()) == params.d);
  CHECK(!proof.extension.base_top_query_multiproof.queried_indices.empty());
  CHECK(proof.extension.base_top_query_multiproof.values.length() > 0);
  CHECK(basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_open_tampered = proof;
  CHECK(!proof_open_tampered.extension.query_multiproofs.empty());
  CHECK(!proof_open_tampered.extension.query_multiproofs[0].values.empty());
  ZZ_pE coeff1 = NTL::coeff(
      proof_open_tampered.extension.query_multiproofs[0].values[0], 1);
  coeff1 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_open_tampered.extension.query_multiproofs[0].values[0], 1,
                coeff1);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_open_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_top_tampered = proof;
  CHECK(proof_top_tampered.extension.base_top_query_multiproof.values.length() > 0);
  proof_top_tampered.extension.base_top_query_multiproof.values[0] +=
      testutil::ConstZZpE(1);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_top_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_root_tampered = proof;
  proof_root_tampered.extension.roots_by_level[0][0] ^=
      static_cast<basefold::Byte>(0x01);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_root_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_h_tampered = proof;
  ZZ_pE h0_a0_coeff0 =
      NTL::coeff(proof_h_tampered.extension.h_by_level[0].a0, 0);
  h0_a0_coeff0 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_h_tampered.extension.h_by_level[0].a0, 0, h0_a0_coeff0);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_h_tampered, params, cfg));

  basefold::BaseFoldPCSEvalProof proof_msg0_tampered = proof;
  ZZ_pE msg0_coeff0 = NTL::coeff(proof_msg0_tampered.extension.msg0_coeffs[0], 0);
  msg0_coeff0 += testutil::ConstZZpE(1);
  NTL::SetCoeff(proof_msg0_tampered.extension.msg0_coeffs[0], 0, msg0_coeff0);
  CHECK(!basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, proof_msg0_tampered, params, cfg));

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
  proof.pi0_codeword.SetLength(2);
  proof.query_multiproofs.resize(1);
  proof.query_multiproofs[0].queried_indices = {0, 1};
  proof.query_multiproofs[0].values.SetLength(2);
  proof.query_multiproofs[0].sibling_hashes.resize(3);
  proof.extension.has_extension_payload = false;

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
  proof.extension.has_extension_payload = true;
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

void TestPCS_ProofFixedWidthRoundTrip_GF4() {
  testutil::PrintInfo(
      "PCS: fixed-width proof deserialize round-trip verifies over GF(2^2)");

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
  const basefold::MerkleRoot commitment_root =
      basefold::BaseFoldPCSCommit(f_coeffs, params);
  const long num_queries = 3;

  const basefold::BaseFoldPCSEvalProof proof =
      basefold::BaseFoldPCSProveEval(f_coeffs, z, y, num_queries, params);

  basefold::FixedProofEncodingOptions options;
  const basefold::Bytes proof_bytes =
      basefold::SerializeBaseFoldPCSEvalProofFixedBytes(proof, options);
  const basefold::BaseFoldPCSEvalProof decoded =
      basefold::DeserializeBaseFoldPCSEvalProofFixedBytes(proof_bytes, options);

  CHECK(decoded.query_multiproofs.size() == proof.query_multiproofs.size());
  for (const basefold::MerkleMultiproof &multiproof : decoded.query_multiproofs) {
    CHECK(multiproof.queried_indices.empty());
  }
  CHECK(!decoded.extension.has_extension_payload);
  CHECK(basefold::BaseFoldPCSVerifyEval(commitment_root, z, y, num_queries,
                                        decoded, params));
}

void TestPCS_ProofFixedWidthRoundTrip_ExtChallenge_GF4() {
  testutil::PrintInfo(
      "PCS: fixed-width proof deserialize round-trip verifies over GF(2^2) with extension challenges");

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
  const basefold::MerkleRoot commitment_root =
      basefold::BaseFoldPCSCommit(f_coeffs, params);

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

  basefold::FixedProofEncodingOptions options;
  options.challenge_ext_degree = NTL::deg(cfg.challenge_extension_modulus);
  const basefold::Bytes proof_bytes =
      basefold::SerializeBaseFoldPCSEvalProofFixedBytes(proof, options);
  const basefold::BaseFoldPCSEvalProof decoded =
      basefold::DeserializeBaseFoldPCSEvalProofFixedBytes(proof_bytes, options);

  CHECK(decoded.extension.has_extension_payload);
  CHECK(decoded.extension.r_by_level.empty());
  CHECK(decoded.extension.query_multiproofs.size() ==
        proof.extension.query_multiproofs.size());
  CHECK(decoded.extension.base_top_query_multiproof.queried_indices.empty());
  for (const basefold::ExtensionMerkleMultiproof &multiproof :
       decoded.extension.query_multiproofs) {
    CHECK(multiproof.queried_indices.empty());
  }
  CHECK(basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
      commitment_root, z, y, num_queries, decoded, params, cfg));
}

basefold_bench_pcs_artifact::ArtifactMetadata MakeSampleArtifactMetadata() {
  basefold_bench_pcs_artifact::ArtifactMetadata meta;
  meta.context_id = "field-default";
  meta.context_label = "Field";
  meta.mode = "field";
  meta.c = 2;
  meta.k0 = 1;
  meta.d = 2;
  meta.poly_dim = basefold_bench_pcs_artifact::ComputePolyDimOrThrow(meta.k0, meta.d);
  meta.lambda = "128";
  meta.gamma = "auto";
  meta.queries = 3;
  meta.seed = 7;
  meta.use_checked_prover_path = false;
  meta.use_extension_challenges = true;
  meta.scalar_modulus = to_ZZ(2);
  meta.base_prime = ZZ(0);
  meta.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};
  meta.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};
  meta.zeta_source = "explicit";
  meta.challenge_extension_coeffs = {
      {to_ZZ(0), to_ZZ(1)}, {to_ZZ(1)}, {to_ZZ(1)}};
  meta.hash_backend = basefold::SelectedHashBackendName();
  meta.proof_encoding = "basefold_fixed_v1";
  meta.proof_size_bytes = 1234;
  meta.artifact_id = basefold_bench_pcs_artifact::ComputeCanonicalArtifactId(meta);
  meta.display_key = basefold_bench_pcs_artifact::BuildArtifactDisplayKey(meta);
  return meta;
}

void TestPCS_ArtifactMetadataJsonRoundTrip() {
  testutil::PrintInfo("PCS artifact: metadata JSON round-trip preserves Phase 2 schema");

  const auto meta = MakeSampleArtifactMetadata();
  const std::string json =
      basefold_bench_pcs_artifact::SerializeMetadataJson(meta);
  const auto decoded = basefold_bench_pcs_artifact::ParseMetadataJson(json);

  CHECK_EQ(decoded.artifact_id, meta.artifact_id);
  CHECK_EQ(decoded.display_key, meta.display_key);
  CHECK_EQ(decoded.context_id, meta.context_id);
  CHECK_EQ(decoded.mode, meta.mode);
  CHECK_EQ(decoded.poly_dim, meta.poly_dim);
  CHECK_EQ(decoded.lambda, meta.lambda);
  CHECK_EQ(decoded.gamma, meta.gamma);
  CHECK_EQ(decoded.queries, meta.queries);
  CHECK_EQ(decoded.seed, meta.seed);
  CHECK(decoded.use_extension_challenges);
  CHECK_EQ(decoded.scalar_modulus, meta.scalar_modulus);
  CHECK_EQ(decoded.base_prime, meta.base_prime);
  CHECK_EQ(decoded.F_coeffs, meta.F_coeffs);
  CHECK_EQ(decoded.zeta_coeffs, meta.zeta_coeffs);
  CHECK_EQ(decoded.challenge_extension_coeffs, meta.challenge_extension_coeffs);
  CHECK_EQ(decoded.hash_backend, meta.hash_backend);
  CHECK_EQ(decoded.proof_encoding, meta.proof_encoding);
  CHECK_EQ(decoded.proof_size_bytes, meta.proof_size_bytes);
}

void TestPCS_ArtifactManifestJsonlRoundTrip() {
  testutil::PrintInfo("PCS artifact: manifest JSONL round-trip preserves deterministic ids and paths");

  const auto meta = MakeSampleArtifactMetadata();
  const auto entry = basefold_bench_pcs_artifact::ManifestEntryFromMetadata(meta);
  const fs::path root = fs::temp_directory_path() / "basefold_phase2_manifest_test";
  fs::remove_all(root);
  fs::create_directories(root);
  const fs::path manifest_path =
      basefold_bench_pcs_artifact::ArtifactManifestPath(root);

  basefold_bench_pcs_artifact::AppendManifestEntry(manifest_path, entry);
  const auto entries = basefold_bench_pcs_artifact::LoadManifestEntries(manifest_path);

  CHECK_EQ(entries.size(), static_cast<std::size_t>(1));
  CHECK_EQ(entries[0].artifact_id, meta.artifact_id);
  CHECK_EQ(entries[0].display_key, meta.display_key);
  CHECK_EQ(entries[0].object_relpath,
           basefold_bench_pcs_artifact::ManifestObjectRelPath(meta.artifact_id));

  fs::remove_all(root);
}

void TestPCS_ArtifactPublicInputsRoundTrip_GF4() {
  testutil::PrintInfo("PCS artifact: public_inputs binary round-trip preserves commitment, z, and y");

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

  basefold_bench_pcs_artifact::ArtifactPublicInputs inputs;
  inputs.commitment_root = basefold::BaseFoldPCSCommit(f_coeffs, params);
  inputs.z = {alpha, alpha + testutil::ConstZZpE(1)};
  inputs.claimed_y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, inputs.z);

  const basefold::Bytes bytes =
      basefold_bench_pcs_artifact::SerializePublicInputs(inputs);
  const auto decoded =
      basefold_bench_pcs_artifact::DeserializePublicInputs(bytes);

  CHECK_EQ(decoded.commitment_root, inputs.commitment_root);
  CHECK_EQ(decoded.z.size(), inputs.z.size());
  for (std::size_t i = 0; i < decoded.z.size(); ++i) {
    CHECK_EQ(decoded.z[i], inputs.z[i]);
  }
  CHECK_EQ(decoded.claimed_y, inputs.claimed_y);
}

void TestPCS_ArtifactRestoreVerificationContext_GF4() {
  testutil::PrintInfo("PCS artifact: metadata restores verifier context and challenge config");

  const auto meta = MakeSampleArtifactMetadata();
  const auto restored =
      basefold_bench_pcs_artifact::RestoreVerificationContext(meta);

  CHECK_EQ(restored.params.c, meta.c);
  CHECK_EQ(restored.params.k0, meta.k0);
  CHECK_EQ(restored.params.d, meta.d);
  CHECK(restored.challenge_cfg.use_extension_challenges);
  CHECK_EQ(deg(restored.challenge_cfg.challenge_extension_modulus), 2);

  const ZZ p = to_ZZ(2);
  ZZ_pPush p_push(p);
  const ZZ_pX F = basefold_bench_pcs_common::BuildZZpX(meta.F_coeffs);
  ZZ_pEPush e_push(F);
  const ZZ_pE expected_zeta = basefold_bench_pcs_common::BuildZZpE(meta.zeta_coeffs);
  CHECK_EQ(restored.params.zeta, expected_zeta);
}

basefold_bench_pcs_common::ContextSpec MakeDefaultArtifactFieldContext() {
  basefold_bench_pcs_common::ContextSpec field;
  field.label = "Field";
  field.scalar_modulus = to_ZZ(2);
  field.base_prime = ZZ(0);
  field.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};
  field.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};
  return field;
}

basefold_bench_pcs_common::ContextSpec MakeDefaultArtifactRingContext() {
  basefold_bench_pcs_common::ContextSpec ring;
  ring.label = "Ring";
  ring.scalar_modulus = to_ZZ(4);
  ring.base_prime = to_ZZ(2);
  ring.F_coeffs = {to_ZZ(1), to_ZZ(1), to_ZZ(1)};
  ring.zeta_coeffs = {to_ZZ(0), to_ZZ(1)};
  return ring;
}

basefold_bench_pcs_artifact::DumpArtifactRequest MakeArtifactRequest(
    const fs::path &root, const string &artifact_id, const string &mode, long d,
    long queries, std::uint64_t seed) {
  basefold_bench_pcs_artifact::DumpArtifactRequest request;
  request.artifact_root = root;
  request.artifact_id = artifact_id;
  request.mode = mode;
  request.c = 2;
  request.k0 = 1;
  request.d = d;
  request.queries = queries;
  request.seed = seed;
  return request;
}

struct DirectArtifactCase {
  basefold_bench_pcs_artifact::ArtifactPublicInputs public_inputs;
  basefold::BaseFoldPCSEvalProof proof;
  basefold::Bytes proof_bytes;
  std::uint64_t proof_size_bytes = 0;
  long challenge_ext_degree = 0;
};

DirectArtifactCase BuildDirectArtifactCase(
    const basefold_bench_pcs_common::ContextSpec &spec,
    const basefold_bench_pcs_artifact::DumpArtifactRequest &request) {
  using namespace basefold_bench_pcs_common;

  ZZ_pPush mod_push(spec.scalar_modulus);
  const ZZ_pX F = BuildZZpX(spec.F_coeffs);
  ZZ_pEPush e_push(F);

  ZZ_pE zeta;
  if (request.auto_zeta_teich) {
    ZZ p_base;
    long k_base = 0;
    DeduceBasePrimeAndExponent(spec, p_base, k_base);
    zeta = FindTeichmullerGenerator(p_base, k_base, NTL::deg(F), F);
  } else {
    zeta = BuildZZpE(spec.zeta_coeffs);
  }

  const basefold::FoldableCodeParams params =
      (request.k0 == 1)
          ? BuildParams_k0_1(
                request.c, request.d,
                (spec.base_prime > 1) ? spec.base_prime : spec.scalar_modulus,
                zeta)
          : BuildParams_k0_pow2(
                request.c, request.k0, request.d,
                (spec.base_prime > 1) ? spec.base_prime : spec.scalar_modulus,
                zeta);

  const long poly_dim =
      basefold_bench_pcs_artifact::ComputePolyDimOrThrow(request.k0, request.d);
  const vec_ZZ_pE f_coeffs = MakeDeterministicCoefficients(poly_dim, request.seed);
  const long point_dim = request.d + Log2ExactPowerOfTwoLong(request.k0);
  const std::vector<ZZ_pE> z =
      MakeDeterministicPoint(point_dim, request.seed ^ 0xdeadbeefULL);
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  basefold::BaseFoldPCSChallengeConfig challenge_cfg;
  const basefold::BaseFoldPCSChallengeConfig *challenge_cfg_ptr = nullptr;
  long challenge_ext_degree = 0;
  if (request.use_extension_challenges) {
    ZZ_pEX challenge_modulus;
    if (spec.challenge_ext_coeffs.empty()) {
      NTL::clear(challenge_modulus);
      NTL::SetCoeff(challenge_modulus, 0, zeta);
      NTL::SetCoeff(challenge_modulus, 1, ZZ_pE(1));
      NTL::SetCoeff(challenge_modulus, 2, ZZ_pE(1));
    } else {
      challenge_modulus = BuildZZpEX(spec.challenge_ext_coeffs);
    }
    challenge_cfg.use_extension_challenges = true;
    challenge_cfg.challenge_extension_modulus = challenge_modulus;
    challenge_cfg_ptr = &challenge_cfg;
    challenge_ext_degree = NTL::deg(challenge_modulus);
  }

  const basefold::BaseFoldPCSCommitArtifacts commit_artifacts =
      basefold::BaseFoldPCSBuildCommitArtifactsUnchecked(f_coeffs, params);
  const basefold::BaseFoldPCSEvalProof proof =
      (challenge_cfg_ptr != nullptr)
          ? (request.use_checked_prover_path
                 ? basefold::BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle(
                       f_coeffs, z, y, request.queries, params,
                       commit_artifacts, *challenge_cfg_ptr)
                 : basefold::BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked(
                       f_coeffs, z, y, request.queries, params,
                       commit_artifacts, *challenge_cfg_ptr))
          : (request.use_checked_prover_path
                 ? basefold::BaseFoldPCSProveEvalFromCommittedTopOracle(
                       f_coeffs, z, y, request.queries, params, commit_artifacts)
                 : basefold::BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
                       f_coeffs, z, y, request.queries, params,
                       commit_artifacts));

  basefold::FixedProofEncodingOptions options;
  options.include_version_byte = true;
  if (request.use_extension_challenges || proof.extension.has_extension_payload) {
    options.challenge_ext_degree = challenge_ext_degree;
  }

  DirectArtifactCase out;
  out.public_inputs.commitment_root = commit_artifacts.root_d;
  out.public_inputs.z = z;
  out.public_inputs.claimed_y = y;
  out.proof = proof;
  out.proof_bytes =
      basefold::SerializeBaseFoldPCSEvalProofFixedBytes(proof, options);
  out.proof_size_bytes = basefold_bench_pcs_common::ComputeProofSizeBytes(
      proof, request.use_extension_challenges, challenge_ext_degree);
  out.challenge_ext_degree = challenge_ext_degree;
  CHECK_EQ(out.proof_bytes.size(), out.proof_size_bytes);
  return out;
}

basefold::FixedProofEncodingOptions MakeProofEncodingOptions(
    const basefold_bench_pcs_artifact::LoadedArtifactCase &loaded) {
  basefold::FixedProofEncodingOptions options;
  options.include_version_byte = true;
  if (loaded.metadata.use_extension_challenges) {
    options.challenge_ext_degree = loaded.challenge_ext_degree;
  }
  return options;
}

void RewriteProofFile(
    const fs::path &path,
    const basefold_bench_pcs_artifact::LoadedArtifactCase &loaded,
    const basefold::BaseFoldPCSEvalProof &proof) {
  const basefold::Bytes proof_bytes =
      basefold::SerializeBaseFoldPCSEvalProofFixedBytes(
          proof, MakeProofEncodingOptions(loaded));
  basefold_bench_pcs_artifact::WriteBytesToFile(path, proof_bytes);
}

void CheckDumpedArtifactVerifies(
    const basefold_bench_pcs_artifact::DumpArtifactResult &result) {
  const auto meta =
      basefold_bench_pcs_artifact::ReadMetadataJson(result.metadata_path);
  CHECK_EQ(meta.artifact_id, result.metadata.artifact_id);
  CHECK_EQ(meta.proof_size_bytes, result.metadata.proof_size_bytes);

  const auto restored =
      basefold_bench_pcs_artifact::RestoreVerificationContext(meta);
  const auto inputs =
      basefold_bench_pcs_artifact::ReadPublicInputsBinary(result.public_inputs_path);
  const basefold::Bytes proof_bytes =
      basefold_bench_pcs_artifact::ReadBytesFromFile(result.proof_path);

  basefold::FixedProofEncodingOptions options;
  options.include_version_byte = true;
  long challenge_ext_degree = 0;
  if (meta.use_extension_challenges) {
    challenge_ext_degree = NTL::deg(
        restored.challenge_cfg.challenge_extension_modulus);
    options.challenge_ext_degree = challenge_ext_degree;
  }
  const basefold::BaseFoldPCSEvalProof proof =
      basefold::DeserializeBaseFoldPCSEvalProofFixedBytes(proof_bytes, options);
  CHECK_EQ(proof_bytes.size(), meta.proof_size_bytes);

  const bool ok =
      meta.use_extension_challenges
          ? basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
                inputs.commitment_root, inputs.z, inputs.claimed_y, meta.queries,
                proof, restored.params, restored.challenge_cfg)
          : basefold::BaseFoldPCSVerifyEval(inputs.commitment_root, inputs.z,
                                            inputs.claimed_y, meta.queries,
                                            proof, restored.params);
  CHECK(ok);
}

void TestPCS_DumpEvalArtifact_FieldCompleteCase() {
  testutil::PrintInfo("PCS artifact: dump tool writes a complete field artifact case");

  const fs::path root =
      fs::temp_directory_path() / "basefold_phase3_dump_field_test";
  fs::remove_all(root);

  const auto field = MakeDefaultArtifactFieldContext();
  auto request = MakeArtifactRequest(root, "", "field", 3, 2, 11);
  request.lambda = "128";
  request.gamma = "auto";

  const auto result =
      basefold_bench_pcs_artifact::DumpEvalArtifact(field, request);
  CHECK(fs::exists(result.manifest_path));
  CHECK(fs::exists(result.metadata_path));
  CHECK(fs::exists(result.public_inputs_path));
  CHECK(fs::exists(result.proof_path));

  const auto entries =
      basefold_bench_pcs_artifact::LoadManifestEntries(result.manifest_path);
  CHECK_EQ(entries.size(), static_cast<std::size_t>(1));
  CHECK_EQ(entries[0].artifact_id, result.metadata.artifact_id);
  CHECK_EQ(entries[0].object_relpath,
           basefold_bench_pcs_artifact::ManifestObjectRelPath(
               result.metadata.artifact_id));

  CheckDumpedArtifactVerifies(result);
  fs::remove_all(root);
}

void TestPCS_DumpEvalArtifact_RingCompleteCase_ExtensionMode() {
  testutil::PrintInfo("PCS artifact: dump tool writes a complete ring artifact case with extension challenges");

  const fs::path root =
      fs::temp_directory_path() / "basefold_phase3_dump_ring_test";
  fs::remove_all(root);

  const auto ring = MakeDefaultArtifactRingContext();
  auto request = MakeArtifactRequest(root, "", "ring", 2, 2, 19);
  request.use_extension_challenges = true;
  request.use_checked_prover_path = true;
  request.lambda = "128";
  request.gamma = "auto";

  const auto result =
      basefold_bench_pcs_artifact::DumpEvalArtifact(ring, request);
  const auto meta =
      basefold_bench_pcs_artifact::ReadMetadataJson(result.metadata_path);
  CHECK(meta.use_extension_challenges);
  CHECK(!meta.challenge_extension_coeffs.empty());
  CHECK_EQ(meta.proof_size_bytes,
           basefold_bench_pcs_artifact::ReadBytesFromFile(result.proof_path)
               .size());

  CheckDumpedArtifactVerifies(result);
  fs::remove_all(root);
}

void TestPCS_LoadArtifactCaseForVerify_FieldBenchmark() {
  testutil::PrintInfo("PCS artifact: verifier benchmark loads one dumped field case and times verify only");

  const fs::path root =
      fs::temp_directory_path() / "basefold_phase4_verify_field_test";
  fs::remove_all(root);

  const auto field = MakeDefaultArtifactFieldContext();
  auto request = MakeArtifactRequest(root, "phase4_field_case", "field", 3, 2, 23);

  const auto dumped =
      basefold_bench_pcs_artifact::DumpEvalArtifact(field, request);
  const auto loaded = basefold_bench_pcs_artifact::LoadArtifactCaseForVerify(
      root, dumped.metadata.artifact_id);
  CHECK_EQ(loaded.metadata.artifact_id, dumped.metadata.artifact_id);
  CHECK(loaded.load_wall_ms >= 0.0);
  CHECK(loaded.deserialize_wall_ms >= 0.0);
  CHECK(basefold_bench_pcs_artifact::VerifyLoadedArtifactCase(loaded));

  const auto bench =
      basefold_bench_pcs_artifact::RunArtifactVerifyBenchmark(
          loaded, /*enable_profile=*/false, /*warmup=*/0, /*reps=*/2);
  CHECK_EQ(bench.proof_size_bytes, dumped.metadata.proof_size_bytes);
  CHECK(bench.verifier.mean_ms >= 0.0);

  fs::remove_all(root);
}

void TestPCS_LoadArtifactCaseForVerify_ExtensionMetadataTamperFails() {
  testutil::PrintInfo("PCS artifact: verifier fails if extension metadata is tampered after load");

  const fs::path root =
      fs::temp_directory_path() / "basefold_phase4_verify_ring_test";
  fs::remove_all(root);

  const auto ring = MakeDefaultArtifactRingContext();
  auto request =
      MakeArtifactRequest(root, "phase4_ring_ext_case", "ring", 2, 2, 29);
  request.use_extension_challenges = true;

  const auto dumped =
      basefold_bench_pcs_artifact::DumpEvalArtifact(ring, request);
  const auto loaded = basefold_bench_pcs_artifact::LoadArtifactCaseForVerify(
      root, dumped.metadata.artifact_id);
  CHECK(basefold_bench_pcs_artifact::VerifyLoadedArtifactCase(loaded));

  auto tampered = loaded;
  CHECK(tampered.restored.challenge_cfg.use_extension_challenges);
  NTL::SetCoeff(tampered.restored.challenge_cfg.challenge_extension_modulus, 0,
                NTL::coeff(
                    tampered.restored.challenge_cfg.challenge_extension_modulus, 0) +
                    testutil::ConstZZpE(1));
  CHECK(!basefold_bench_pcs_artifact::VerifyLoadedArtifactCase(tampered));

  fs::remove_all(root);
}

void TestPCS_ArtifactCompareDirectAndArtifactProof_Field() {
  testutil::PrintInfo("PCS artifact: direct verify input matches artifact verify input for the same field case");

  const fs::path root =
      fs::temp_directory_path() / "basefold_phase5_compare_field_test";
  fs::remove_all(root);

  const auto field = MakeDefaultArtifactFieldContext();
  const auto request =
      MakeArtifactRequest(root, "phase5_compare_field_case", "field", 3, 2, 31);

  const auto dumped =
      basefold_bench_pcs_artifact::DumpEvalArtifact(field, request);
  const auto loaded = basefold_bench_pcs_artifact::LoadArtifactCaseForVerify(
      root, dumped.metadata.artifact_id);
  const auto direct = BuildDirectArtifactCase(field, request);

  CHECK_EQ(loaded.public_inputs.commitment_root, direct.public_inputs.commitment_root);
  CHECK_EQ(loaded.public_inputs.z.size(), direct.public_inputs.z.size());
  for (std::size_t i = 0; i < loaded.public_inputs.z.size(); ++i) {
    CHECK_EQ(loaded.public_inputs.z[i], direct.public_inputs.z[i]);
  }
  CHECK_EQ(loaded.public_inputs.claimed_y, direct.public_inputs.claimed_y);
  CHECK_EQ(loaded.metadata.proof_size_bytes, direct.proof_size_bytes);
  CHECK_EQ(basefold_bench_pcs_artifact::ReadBytesFromFile(dumped.proof_path),
           direct.proof_bytes);

  fs::remove_all(root);
}

void TestPCS_ArtifactVerifyFailsWhenProofFileTampered() {
  testutil::PrintInfo("PCS artifact: verifier fails when proof.bin is tampered");

  const fs::path root =
      fs::temp_directory_path() / "basefold_phase5_tamper_proof_test";
  fs::remove_all(root);

  const auto field = MakeDefaultArtifactFieldContext();
  const auto request =
      MakeArtifactRequest(root, "phase5_tamper_proof_case", "field", 3, 2, 37);

  const auto dumped =
      basefold_bench_pcs_artifact::DumpEvalArtifact(field, request);
  const auto loaded = basefold_bench_pcs_artifact::LoadArtifactCaseForVerify(
      root, dumped.metadata.artifact_id);
  CHECK(basefold_bench_pcs_artifact::VerifyLoadedArtifactCase(loaded));

  auto proof_tampered = loaded.proof;
  CHECK(!proof_tampered.query_multiproofs.empty());
  CHECK(proof_tampered.query_multiproofs[0].values.length() > 0);
  proof_tampered.query_multiproofs[0].values[0] += testutil::ConstZZpE(1);
  RewriteProofFile(dumped.proof_path, loaded, proof_tampered);

  const auto reloaded = basefold_bench_pcs_artifact::LoadArtifactCaseForVerify(
      root, dumped.metadata.artifact_id);
  CHECK(!basefold_bench_pcs_artifact::VerifyLoadedArtifactCase(reloaded));

  fs::remove_all(root);
}

void TestPCS_ArtifactVerifyFailsWhenCommitmentRootTampered() {
  testutil::PrintInfo("PCS artifact: verifier fails when commitment_root in public_inputs.bin is tampered");

  const fs::path root =
      fs::temp_directory_path() / "basefold_phase5_tamper_root_test";
  fs::remove_all(root);

  const auto field = MakeDefaultArtifactFieldContext();
  const auto request =
      MakeArtifactRequest(root, "phase5_tamper_root_case", "field", 3, 2, 41);

  const auto dumped =
      basefold_bench_pcs_artifact::DumpEvalArtifact(field, request);
  auto inputs =
      basefold_bench_pcs_artifact::ReadPublicInputsBinary(dumped.public_inputs_path);
  inputs.commitment_root[0] ^= static_cast<basefold::Byte>(0x01);
  basefold_bench_pcs_artifact::WritePublicInputsBinary(dumped.public_inputs_path,
                                                       inputs);

  const auto reloaded = basefold_bench_pcs_artifact::LoadArtifactCaseForVerify(
      root, dumped.metadata.artifact_id);
  CHECK(!basefold_bench_pcs_artifact::VerifyLoadedArtifactCase(reloaded));

  fs::remove_all(root);
}

void TestPCS_ArtifactVerifyFailsWhenPointZTampered() {
  testutil::PrintInfo("PCS artifact: verifier fails when point z in public_inputs.bin is tampered");

  const fs::path root =
      fs::temp_directory_path() / "basefold_phase5_tamper_z_test";
  fs::remove_all(root);

  const auto field = MakeDefaultArtifactFieldContext();
  const auto request =
      MakeArtifactRequest(root, "phase5_tamper_z_case", "field", 3, 2, 43);

  const auto dumped =
      basefold_bench_pcs_artifact::DumpEvalArtifact(field, request);
  auto inputs =
      basefold_bench_pcs_artifact::ReadPublicInputsBinary(dumped.public_inputs_path);
  CHECK(!inputs.z.empty());
  inputs.z[0] += testutil::ConstZZpE(1);
  basefold_bench_pcs_artifact::WritePublicInputsBinary(dumped.public_inputs_path,
                                                       inputs);

  const auto reloaded = basefold_bench_pcs_artifact::LoadArtifactCaseForVerify(
      root, dumped.metadata.artifact_id);
  CHECK(!basefold_bench_pcs_artifact::VerifyLoadedArtifactCase(reloaded));

  fs::remove_all(root);
}

void TestPCS_ArtifactVerifyFailsWhenClaimedYTampered() {
  testutil::PrintInfo("PCS artifact: verifier fails when claimed y in public_inputs.bin is tampered");

  const fs::path root =
      fs::temp_directory_path() / "basefold_phase5_tamper_y_test";
  fs::remove_all(root);

  const auto field = MakeDefaultArtifactFieldContext();
  const auto request =
      MakeArtifactRequest(root, "phase5_tamper_y_case", "field", 3, 2, 47);

  const auto dumped =
      basefold_bench_pcs_artifact::DumpEvalArtifact(field, request);
  auto inputs =
      basefold_bench_pcs_artifact::ReadPublicInputsBinary(dumped.public_inputs_path);
  inputs.claimed_y += testutil::ConstZZpE(1);
  basefold_bench_pcs_artifact::WritePublicInputsBinary(dumped.public_inputs_path,
                                                       inputs);

  const auto reloaded = basefold_bench_pcs_artifact::LoadArtifactCaseForVerify(
      root, dumped.metadata.artifact_id);
  CHECK(!basefold_bench_pcs_artifact::VerifyLoadedArtifactCase(reloaded));

  fs::remove_all(root);
}

}  // namespace

int main() {
  try {
    RUN_TEST(TestPCS_EqualityTableFromPoint_MatchesEqPolynomialDefinition);
    RUN_TEST(TestPCS_MultiproofValueCache_ReplaysExpectedPositions);
    RUN_TEST(TestPCS_EvalProof_GF4);
    RUN_TEST(TestPCS_PaperAPI_GF4);
    RUN_TEST(TestPCS_EvalProof_GF4_k0_2);
    RUN_TEST(TestPCS_EvalProof_GR42);
    RUN_TEST(TestPCS_EvalProof_GR42_k0_2);
    RUN_TEST(TestPCS_EvalProof_ExtChallengeConfig_GF4);
    RUN_TEST(TestPCS_EvalProof_ExtChallengeConfig_GR42);
    RUN_TEST(TestPCS_ProofSizeFixedWidth_GF4_HandCheck);
    RUN_TEST(TestPCS_ProofSizeFixedWidth_ExtensionWidthDerivation);
    RUN_TEST(TestPCS_ProofFixedWidthRoundTrip_GF4);
    RUN_TEST(TestPCS_ProofFixedWidthRoundTrip_ExtChallenge_GF4);
    RUN_TEST(TestPCS_ArtifactMetadataJsonRoundTrip);
    RUN_TEST(TestPCS_ArtifactManifestJsonlRoundTrip);
    RUN_TEST(TestPCS_ArtifactPublicInputsRoundTrip_GF4);
    RUN_TEST(TestPCS_ArtifactRestoreVerificationContext_GF4);
    RUN_TEST(TestPCS_DumpEvalArtifact_FieldCompleteCase);
    RUN_TEST(TestPCS_DumpEvalArtifact_RingCompleteCase_ExtensionMode);
    RUN_TEST(TestPCS_LoadArtifactCaseForVerify_FieldBenchmark);
    RUN_TEST(TestPCS_LoadArtifactCaseForVerify_ExtensionMetadataTamperFails);
    RUN_TEST(TestPCS_ArtifactCompareDirectAndArtifactProof_Field);
    RUN_TEST(TestPCS_ArtifactVerifyFailsWhenProofFileTampered);
    RUN_TEST(TestPCS_ArtifactVerifyFailsWhenCommitmentRootTampered);
    RUN_TEST(TestPCS_ArtifactVerifyFailsWhenPointZTampered);
    RUN_TEST(TestPCS_ArtifactVerifyFailsWhenClaimedYTampered);
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
