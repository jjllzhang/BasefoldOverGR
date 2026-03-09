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

#include "PCS/BaseFold/FoldableCode.hpp"
#include "PCS/BaseFold/IOPP.hpp"
#include "tests/test_common.hpp"

using NTL::conv;
using NTL::mat_ZZ_pE;
using NTL::mul;
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
using std::ostringstream;
using std::string;

int g_test_failure_count = 0;

namespace {

ZZ_pE LookupGF4ElementByCode(int idx, const ZZ_pE &alpha) {
  const ZZ_pE one = testutil::ConstZZpE(1);
  if (idx == 0) return ZZ_pE(0);
  if (idx == 1) return one;
  if (idx == 2) return alpha;
  return alpha + one;
}

basefold::FoldableCodeParams BuildParamsGF4(const ZZ &p, const ZZ_pE &alpha) {
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

  mat_ZZ_pE G0;
  G0.SetDims(k0, n0);
  const ZZ_pE one = testutil::ConstZZpE(1);
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

void TestIOPP_CommitAndQuery() {
  testutil::PrintInfo("IOPP: commit/query consistency over GF(2^2)");

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

  const basefold::FoldableCodeParams params = BuildParamsGF4(p, alpha);

  vec_ZZ_pE msg;
  msg.SetLength(basefold::MessageLength(params));
  for (long i = 0; i < msg.length(); ++i) {
    msg[i] = LookupGF4ElementByCode(static_cast<int>(i % 4), alpha);
  }

  vec_ZZ_pE pi_d;
  basefold::EncodeFoldable(pi_d, msg, params);

  basefold::IOPPChallenges challenges;
  challenges.alphas.resize(static_cast<std::size_t>(params.d));
  challenges.alphas[0] = alpha;
  challenges.alphas[1] = alpha + testutil::ConstZZpE(1);

  basefold::IOPPOracles oracles;
  basefold::ProverCommitAll(oracles, pi_d, challenges, params);

  const long mu0 = 6;
  const basefold::IOPPQueryPlan plan = basefold::MakeQueryPlan(mu0, params);

  CHECK(basefold::VerifyQueryFromOracles(plan, challenges, oracles, params));

  basefold::IOPPQueryOpenings openings;
  openings.upper_left_by_level.resize(static_cast<std::size_t>(params.d));
  openings.upper_right_by_level.resize(static_cast<std::size_t>(params.d));
  openings.folded_by_level.resize(static_cast<std::size_t>(params.d));
  for (long i = 0; i < params.d; ++i) {
    const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
    const long ni = basefold::CodewordLengthAtLevel(params, i);
    openings.upper_left_by_level[static_cast<std::size_t>(i)] =
        oracles.oracles_by_level[static_cast<std::size_t>(i + 1)][mu];
    openings.upper_right_by_level[static_cast<std::size_t>(i)] =
        oracles.oracles_by_level[static_cast<std::size_t>(i + 1)][mu + ni];
    openings.folded_by_level[static_cast<std::size_t>(i)] =
        oracles.oracles_by_level[static_cast<std::size_t>(i)][mu];
  }
  openings.pi0_codeword = oracles.oracles_by_level[0];

  CHECK(basefold::VerifyQueryFromOpenings(plan, challenges, openings, params));
}

void TestIOPP_GR42_CommitAndQuery() {
  testutil::PrintInfo("IOPP: commit/query consistency over GR(4,2)");

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

  vec_ZZ_pE msg;
  msg.SetLength(basefold::MessageLength(params));
  msg[0] = testutil::ConstZZpE(0);
  msg[1] = testutil::ConstZZpE(1);
  msg[2] = alpha;
  msg[3] = testutil::ConstZZpE(2);

  vec_ZZ_pE pi_d;
  basefold::EncodeFoldable(pi_d, msg, params);

  basefold::IOPPChallenges challenges;
  challenges.alphas.resize(static_cast<std::size_t>(params.d));
  challenges.alphas[0] = alpha;
  challenges.alphas[1] = alpha + testutil::ConstZZpE(1);

  basefold::IOPPOracles oracles;
  basefold::ProverCommitAll(oracles, pi_d, challenges, params);

  const long mu0 = 3;
  const basefold::IOPPQueryPlan plan = basefold::MakeQueryPlan(mu0, params);
  CHECK(basefold::VerifyQueryFromOracles(plan, challenges, oracles, params));

  basefold::IOPPQueryOpenings openings;
  openings.upper_left_by_level.resize(static_cast<std::size_t>(params.d));
  openings.upper_right_by_level.resize(static_cast<std::size_t>(params.d));
  openings.folded_by_level.resize(static_cast<std::size_t>(params.d));
  for (long i = 0; i < params.d; ++i) {
    const long mu = plan.mu_by_level[static_cast<std::size_t>(i)];
    const long ni = basefold::CodewordLengthAtLevel(params, i);
    openings.upper_left_by_level[static_cast<std::size_t>(i)] =
        oracles.oracles_by_level[static_cast<std::size_t>(i + 1)][mu];
    openings.upper_right_by_level[static_cast<std::size_t>(i)] =
        oracles.oracles_by_level[static_cast<std::size_t>(i + 1)][mu + ni];
    openings.folded_by_level[static_cast<std::size_t>(i)] =
        oracles.oracles_by_level[static_cast<std::size_t>(i)][mu];
  }
  openings.pi0_codeword = oracles.oracles_by_level[0];
  CHECK(basefold::VerifyQueryFromOpenings(plan, challenges, openings, params));
}

void TestIOPP_MerkleMultiproof() {
  testutil::PrintInfo("IOPP: Merkle multiproof verifies and rejects tampering");

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

  basefold::Oracle oracle;
  oracle.SetLength(5);
  oracle[0] = testutil::ConstZZpE(0);
  oracle[1] = testutil::ConstZZpE(1);
  oracle[2] = alpha;
  oracle[3] = alpha + testutil::ConstZZpE(1);
  oracle[4] = testutil::ConstZZpE(1);

  const basefold::MerkleRoot root = basefold::MerkleCommitOracle(oracle);
  const std::vector<long> queried = {4, 1, 1, 3};
  const basefold::MerkleMultiproof proof =
      basefold::MerkleOpenOracleMany(oracle, queried);

  CHECK(proof.queried_indices.size() == 3U);
  CHECK(proof.queried_indices[0] == 1);
  CHECK(proof.queried_indices[1] == 3);
  CHECK(proof.queried_indices[2] == 4);
  CHECK(basefold::MerkleVerifyMultiproof(root, oracle.length(), proof));

  basefold::MerkleMultiproof proof_value_tampered = proof;
  proof_value_tampered.values[0] += testutil::ConstZZpE(1);
  CHECK(!basefold::MerkleVerifyMultiproof(root, oracle.length(),
                                          proof_value_tampered));

  basefold::MerkleMultiproof proof_hash_tampered = proof;
  CHECK(!proof_hash_tampered.sibling_hashes.empty());
  proof_hash_tampered.sibling_hashes[0][0] ^=
      static_cast<basefold::Byte>(0x01);
  CHECK(!basefold::MerkleVerifyMultiproof(root, oracle.length(),
                                          proof_hash_tampered));

  basefold::MerkleMultiproof proof_order_tampered = proof;
  std::swap(proof_order_tampered.queried_indices[0],
            proof_order_tampered.queried_indices[1]);
  CHECK(!basefold::MerkleVerifyMultiproof(root, oracle.length(),
                                          proof_order_tampered));

  CHECK(!basefold::MerkleVerifyMultiproof(root, oracle.length() - 1, proof));
}

void TestIOPP_DecodeC0() {
  testutil::PrintInfo("IOPP: DecodeC0 recovers a witness for pi0");

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

  const basefold::FoldableCodeParams params = BuildParamsGF4(p, alpha);

  vec_ZZ_pE msg;
  msg.SetLength(basefold::MessageLength(params));
  for (long i = 0; i < msg.length(); ++i) {
    msg[i] = LookupGF4ElementByCode(static_cast<int>((i + 2) % 4), alpha);
  }

  vec_ZZ_pE pi_d;
  basefold::EncodeFoldable(pi_d, msg, params);

  basefold::IOPPChallenges challenges;
  challenges.alphas.resize(static_cast<std::size_t>(params.d));
  challenges.alphas[0] = alpha;
  challenges.alphas[1] = alpha + testutil::ConstZZpE(1);

  basefold::IOPPOracles oracles;
  basefold::ProverCommitAll(oracles, pi_d, challenges, params);

  vec_ZZ_pE msg0;
  CHECK(basefold::DecodeC0(msg0, oracles.oracles_by_level[0], params));
  vec_ZZ_pE rec;
  mul(rec, msg0, params.G0);
  CHECK_EQ(rec, oracles.oracles_by_level[0]);
  CHECK(basefold::IsCodewordC0(oracles.oracles_by_level[0], params));
}

void TestIOPP_GR42_DecodeC0() {
  testutil::PrintInfo("IOPP: DecodeC0 over GR(4,2)");

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

  vec_ZZ_pE msg;
  msg.SetLength(basefold::MessageLength(params));
  msg[0] = testutil::ConstZZpE(1);
  msg[1] = testutil::ConstZZpE(2);
  msg[2] = alpha;
  msg[3] = alpha + testutil::ConstZZpE(1);

  vec_ZZ_pE pi_d;
  basefold::EncodeFoldable(pi_d, msg, params);

  basefold::IOPPChallenges challenges;
  challenges.alphas.resize(static_cast<std::size_t>(params.d));
  challenges.alphas[0] = alpha;
  challenges.alphas[1] = alpha + testutil::ConstZZpE(1);

  basefold::IOPPOracles oracles;
  basefold::ProverCommitAll(oracles, pi_d, challenges, params);

  vec_ZZ_pE msg0;
  CHECK(basefold::DecodeC0(msg0, oracles.oracles_by_level[0], params));
  vec_ZZ_pE rec;
  mul(rec, msg0, params.G0);
  CHECK_EQ(rec, oracles.oracles_by_level[0]);
  CHECK(basefold::IsCodewordC0(oracles.oracles_by_level[0], params));
}

}  // namespace

int main() {
  try {
    RUN_TEST(TestIOPP_CommitAndQuery);
    RUN_TEST(TestIOPP_GR42_CommitAndQuery);
    RUN_TEST(TestIOPP_MerkleMultiproof);
    RUN_TEST(TestIOPP_DecodeC0);
    RUN_TEST(TestIOPP_GR42_DecodeC0);
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
