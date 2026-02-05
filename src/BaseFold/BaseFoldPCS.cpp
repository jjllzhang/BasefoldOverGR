#include "BaseFold/BaseFoldPCS.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ_pE.h>

#include <openssl/evp.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "BaseFold/Multilinear.hpp"
#include "BaseFold/Profile.hpp"

using NTL::BytesFromZZ;
using NTL::coeff;
using NTL::LogicError;
using NTL::mat_ZZ_pE;
using NTL::mul;
using NTL::NumBits;
using NTL::NumBytes;
using NTL::rep;
using NTL::vec_ZZ_pE;
using NTL::ZZ;
using NTL::ZZFromBytes;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pX;

namespace basefold {
namespace {

long Pow2Checked(long e) {
  if (e < 0)
    LogicError("Pow2Checked: negative exponent");
  if (e >= static_cast<long>(8 * sizeof(long) - 1)) {
    LogicError("Pow2Checked: exponent too large for long");
  }
  return 1L << e;
}

long CodewordLengthAtLevelNoValidate(const FoldableCodeParams &params,
                                     long level) {
  if (level < 0 || level > params.d) {
    LogicError("CodewordLengthAtLevelNoValidate: level out of range");
  }
  const long pow2 = Pow2Checked(level);
  if (params.k0 <= 0 || params.c <= 0) {
    LogicError("CodewordLengthAtLevelNoValidate: invalid c/k0");
  }
  if (params.k0 > std::numeric_limits<long>::max() / pow2) {
    LogicError("CodewordLengthAtLevelNoValidate: overflow");
  }
  const long k = params.k0 * pow2;
  if (params.c > std::numeric_limits<long>::max() / k) {
    LogicError("CodewordLengthAtLevelNoValidate: overflow");
  }
  return params.c * k;
}

void ProverCommitRoundNoValidate(Oracle &pi_i, const Oracle &pi_ip1,
                                 const FieldElement &alpha_i, long level_i,
                                 const FoldableCodeParams &params) {
  const long n_i = CodewordLengthAtLevelNoValidate(params, level_i);
  pi_i.SetLength(n_i);

  for (long j = 0; j < n_i; ++j) {
    const FieldElement &t = params.diag_T[static_cast<std::size_t>(level_i)][j];
    const FieldElement x1 = t;
    const FieldElement x2 = params.zeta * t;

    const FieldElement &y1 = pi_ip1[j];
    const FieldElement &y2 = pi_ip1[j + n_i];
    pi_i[j] = EvalLineAt(alpha_i, x1, y1, x2, y2);
  }
}

IOPPQueryPlan MakeQueryPlanNoValidate(long initial_mu,
                                     const FoldableCodeParams &params) {
  IOPPQueryPlan plan;
  plan.initial_mu = initial_mu;
  plan.mu_by_level.resize(static_cast<std::size_t>(params.d));

  if (params.d == 0)
    return plan;

  const long n_last = CodewordLengthAtLevelNoValidate(params, params.d - 1);
  if (initial_mu < 0 || initial_mu >= n_last) {
    LogicError("MakeQueryPlanNoValidate: initial_mu out of range");
  }

  long mu = initial_mu;
  for (long i = params.d; i-- > 0;) {
    plan.mu_by_level[static_cast<std::size_t>(i)] = mu;
    if (i > 0) {
      const long n_prev = CodewordLengthAtLevelNoValidate(params, i - 1);
      if (mu >= n_prev)
        mu -= n_prev;
    }
  }
  return plan;
}

void AppendU64(Bytes &out, std::uint64_t v) {
  for (int i = 7; i >= 0; --i) {
    out.push_back(static_cast<Byte>((v >> (8 * i)) & 0xff));
  }
}

Bytes ZZToBytes(const ZZ &x) {
  const long n = NumBytes(x);
  Bytes out;
  out.resize(static_cast<std::size_t>(n));
  if (n > 0) {
    BytesFromZZ(reinterpret_cast<unsigned char *>(out.data()), x, n);
  }
  return out;
}

Bytes SerializeFieldElement(const FieldElement &x) {
  const long r = ZZ_pE::degree();
  if (r <= 0)
    LogicError("SerializeFieldElement: invalid extension degree");

  const ZZ_pX &poly = rep(x);
  Bytes out;
  AppendU64(out, static_cast<std::uint64_t>(r));
  for (long i = 0; i < r; ++i) {
    const ZZ c = rep(coeff(poly, i));
    const Bytes c_bytes = ZZToBytes(c);
    AppendU64(out, static_cast<std::uint64_t>(c_bytes.size()));
    out.insert(out.end(), c_bytes.begin(), c_bytes.end());
  }
  return out;
}

Bytes Sha256(const Byte *data, std::size_t len) {
  Bytes digest;
  digest.resize(32);
  unsigned int out_len = 0;
  const int ok = EVP_Digest(static_cast<const void *>(data), len,
                            reinterpret_cast<unsigned char *>(digest.data()),
                            &out_len, EVP_sha256(), nullptr);
  if (ok != 1) {
    LogicError("Sha256: EVP_Digest failed");
  }
  if (out_len != digest.size()) {
    LogicError("Sha256: unexpected digest size");
  }
  return digest;
}

Bytes Sha256(const Bytes &data) { return Sha256(data.data(), data.size()); }

Bytes TaggedHash(Byte tag, const Bytes &state, const Bytes &payload) {
  Bytes in;
  in.reserve(1 + state.size() + 8 + payload.size());
  in.push_back(tag);
  in.insert(in.end(), state.begin(), state.end());
  AppendU64(in, static_cast<std::uint64_t>(payload.size()));
  in.insert(in.end(), payload.begin(), payload.end());
  return Sha256(in);
}

Bytes TaggedHash(Byte tag, const Bytes &state, const std::string &payload) {
  Bytes p;
  p.reserve(payload.size());
  p.insert(p.end(), payload.begin(), payload.end());
  return TaggedHash(tag, state, p);
}

class Sha256Transcript {
 public:
  Sha256Transcript() {
    const std::string domain = "BaseFoldPCS/v1";
    state_ = TaggedHash(static_cast<Byte>(0x42), Bytes{}, domain);
  }

  void AbsorbBytes(const Bytes &data) {
    state_ = TaggedHash(static_cast<Byte>(0x01), state_, data);
  }

  void AbsorbDigest(const Digest &digest) {
    const Bytes tmp(digest.begin(), digest.end());
    AbsorbBytes(tmp);
  }

  void AbsorbFieldElement(const FieldElement &x) {
    state_ = TaggedHash(static_cast<Byte>(0x02), state_, SerializeFieldElement(x));
  }

  void AbsorbQuadraticPoly(const QuadraticPoly &p) {
    AbsorbFieldElement(p.a0);
    AbsorbFieldElement(p.a1);
    AbsorbFieldElement(p.a2);
  }

  FieldElement ChallengeFieldElement(const std::string &label) const {
    const long r = ZZ_pE::degree();
    if (r <= 0)
      LogicError("ChallengeFieldElement: invalid extension degree");

    const ZZ modulus = ZZ_p::modulus();
    if (modulus <= 1)
      LogicError("ChallengeFieldElement: invalid base modulus");

    ChallengeStream stream(state_, "fe/" + label);
    ZZ_pX poly;
    NTL::clear(poly);
    for (long i = 0; i < r; ++i) {
      const ZZ c = stream.SampleZZLessThan(modulus);
      ZZ_p c_base;
      NTL::conv(c_base, c);
      NTL::SetCoeff(poly, i, c_base);
    }
    FieldElement out;
    NTL::conv(out, poly);
    return out;
  }

  long ChallengeIndex(const std::string &label, long upper_bound) const {
    if (upper_bound <= 0)
      LogicError("ChallengeIndex: upper_bound must be positive");
    if (upper_bound == 1)
      return 0;

    ChallengeStream stream(state_, "idx/" + label);

    std::uint64_t ub = static_cast<std::uint64_t>(upper_bound);
    std::uint64_t t = ub - 1;
    int bits = 0;
    while (t > 0) {
      ++bits;
      t >>= 1;
    }
    if (bits <= 0 || bits > 63)
      LogicError("ChallengeIndex: unsupported range");
    const int byte_len = (bits + 7) / 8;
    const std::uint64_t mask =
        (bits == 64) ? ~0ULL : ((1ULL << bits) - 1ULL);

    while (true) {
      std::uint8_t buf[8];
      std::memset(buf, 0, sizeof(buf));
      stream.ReadBytes(buf, static_cast<std::size_t>(byte_len));
      std::uint64_t x = 0;
      for (int i = 0; i < byte_len; ++i) {
        x = (x << 8) | static_cast<std::uint64_t>(buf[i]);
      }
      x &= mask;
      if (x < ub)
        return static_cast<long>(x);
    }
  }

 private:
  class ChallengeStream {
   public:
    ChallengeStream(const Bytes &state, const std::string &label)
        : state_(state), label_(label) {}

    void ReadBytes(std::uint8_t *out, std::size_t len) {
      std::size_t written = 0;
      while (written < len) {
        if (offset_ == buf_.size()) {
          buf_ = Digest(counter_++);
          offset_ = 0;
        }
        const std::size_t take = std::min(len - written, buf_.size() - offset_);
        std::memcpy(out + written, buf_.data() + offset_, take);
        offset_ += take;
        written += take;
      }
    }

    ZZ SampleZZLessThan(const ZZ &upper_bound) {
      if (upper_bound <= 0)
        LogicError("SampleZZLessThan: upper_bound must be positive");
      if (upper_bound == 1)
        return ZZ(0);

      const ZZ ub_minus_1 = upper_bound - 1;
      const long bits = NumBits(ub_minus_1);
      if (bits <= 0)
        return ZZ(0);
      const long byte_len = (bits + 7) / 8;
      const ZZ two_to_bits = ZZ(1) << bits;

      Bytes tmp;
      tmp.resize(static_cast<std::size_t>(byte_len));
      while (true) {
        ReadBytes(reinterpret_cast<std::uint8_t *>(tmp.data()),
                  static_cast<std::size_t>(byte_len));
        ZZ x = ZZFromBytes(reinterpret_cast<const unsigned char *>(tmp.data()),
                           byte_len);
        x %= two_to_bits;
        if (x < upper_bound)
          return x;
      }
    }

   private:
    Bytes Digest(std::uint64_t ctr) const {
      Bytes payload;
      payload.reserve(8);
      AppendU64(payload, ctr);
      Bytes st = TaggedHash(static_cast<Byte>(0x20), state_, label_);
      return TaggedHash(static_cast<Byte>(0x21), st, payload);
    }

    Bytes state_;
    std::string label_;

    std::uint64_t counter_ = 0;
    Bytes buf_;
    std::size_t offset_ = 0;
  };

  Bytes state_;
};

void ValidateParamsOrThrow(const FoldableCodeParams &params) {
  (void)MessageLength(params);
}

void AbsorbPublicInput(Sha256Transcript &transcript, const MerkleRoot &commitment,
                       const std::vector<FieldElement> &z,
                       const FieldElement &y) {
  transcript.AbsorbDigest(commitment);
  for (const FieldElement &zi : z) {
    transcript.AbsorbFieldElement(zi);
  }
  transcript.AbsorbFieldElement(y);
}

}  // namespace

MerkleRoot BaseFoldPCSCommit(const vec_ZZ_pE &f_coeffs,
                             const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (f_coeffs.length() != MessageLength(params))
    LogicError("BaseFoldPCSCommit: f_coeffs has wrong length");

  Oracle pi_d;
  EncodeFoldable(pi_d, f_coeffs, params);
  return MerkleCommitOracle(pi_d);
}

BaseFoldPCSEvalProof BaseFoldPCSProveEval(const vec_ZZ_pE &f_coeffs,
                                          const std::vector<FieldElement> &z,
                                          const FieldElement &claimed_y,
                                          long num_queries,
                                          const FoldableCodeParams &params) {
  ValidateParamsOrThrow(params);
  if (params.k0 != 1)
    LogicError("BaseFoldPCSProveEval: only supports k0 == 1");
  if (static_cast<long>(z.size()) != params.d)
    LogicError("BaseFoldPCSProveEval: z has wrong dimension");
  if (f_coeffs.length() != MessageLength(params))
    LogicError("BaseFoldPCSProveEval: f_coeffs has wrong length");
  if (num_queries < 0)
    LogicError("BaseFoldPCSProveEval: num_queries must be non-negative");

  if (EvalMultilinearMonomialCoeffs(f_coeffs, z) != claimed_y) {
    LogicError("BaseFoldPCSProveEval: claimed_y != f(z)");
  }

  return BaseFoldPCSProveEvalUnchecked(f_coeffs, z, claimed_y, num_queries,
                                      params);
}

BaseFoldPCSEvalProof BaseFoldPCSProveEvalUnchecked(
    const vec_ZZ_pE &f_coeffs, const std::vector<FieldElement> &z,
    const FieldElement &claimed_y, long num_queries,
    const FoldableCodeParams &params) {
  BaseFoldPCSEvalProof proof;
  proof.commitments.roots_by_level.resize(static_cast<std::size_t>(params.d + 1));
  proof.h_by_level.resize(static_cast<std::size_t>(params.d));

  IOPPOracles oracles;
  oracles.pi.resize(static_cast<std::size_t>(params.d + 1));

  std::vector<MerkleTree> merkle;
  merkle.resize(static_cast<std::size_t>(params.d + 1));

  EncodeFoldableUnchecked(oracles.pi[static_cast<std::size_t>(params.d)], f_coeffs,
                          params);
  merkle[static_cast<std::size_t>(params.d)] =
      MerkleTree::Build(oracles.pi[static_cast<std::size_t>(params.d)]);
  const MerkleRoot root_d =
      merkle[static_cast<std::size_t>(params.d)].Root();
  proof.commitments.roots_by_level[static_cast<std::size_t>(params.d)] = root_d;

  Sha256Transcript transcript;
  AbsorbPublicInput(transcript, root_d, z, claimed_y);

  SumcheckProver sumcheck(f_coeffs, z);

  // h_d
  const QuadraticPoly h_d = sumcheck.CurrentPolynomial();
  proof.h_by_level[static_cast<std::size_t>(params.d - 1)] = h_d;
  transcript.AbsorbQuadraticPoly(h_d);

  std::vector<FieldElement> r_by_level;
  r_by_level.resize(static_cast<std::size_t>(params.d));

  for (long i = params.d; i-- > 0;) {
    const FieldElement r_i =
        transcript.ChallengeFieldElement("r/" + std::to_string(i));
    r_by_level[static_cast<std::size_t>(i)] = r_i;

    ProverCommitRoundNoValidate(oracles.pi[static_cast<std::size_t>(i)],
                                oracles.pi[static_cast<std::size_t>(i + 1)],
                                r_i, i, params);

    merkle[static_cast<std::size_t>(i)] =
        MerkleTree::Build(oracles.pi[static_cast<std::size_t>(i)]);
    const MerkleRoot root_i = merkle[static_cast<std::size_t>(i)].Root();
    proof.commitments.roots_by_level[static_cast<std::size_t>(i)] = root_i;
    transcript.AbsorbDigest(root_i);

    sumcheck.ReceiveChallenge(r_i);
    if (i > 0) {
      const QuadraticPoly h_i = sumcheck.CurrentPolynomial();
      proof.h_by_level[static_cast<std::size_t>(i - 1)] = h_i;
      transcript.AbsorbQuadraticPoly(h_i);
    }
  }

  proof.pi0_full = oracles.pi[0];

  proof.query_proofs.resize(static_cast<std::size_t>(num_queries));
  if (params.d == 0) {
    return proof;
  }

  const long n_last = CodewordLengthAtLevelNoValidate(params, params.d - 1);
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    const IOPPQueryPlan plan = MakeQueryPlanNoValidate(mu, params);

    BaseFoldPCSQueryProof qp;
    qp.left.resize(static_cast<std::size_t>(params.d));
    qp.right.resize(static_cast<std::size_t>(params.d));
    qp.folded.resize(static_cast<std::size_t>(params.d));

    for (long i = 0; i < params.d; ++i) {
      const long mu_i = plan.mu_by_level[static_cast<std::size_t>(i)];
      const long n_i = CodewordLengthAtLevelNoValidate(params, i);
      qp.left[static_cast<std::size_t>(i)] =
          merkle[static_cast<std::size_t>(i + 1)].Open(
              oracles.pi[static_cast<std::size_t>(i + 1)], mu_i);
      qp.right[static_cast<std::size_t>(i)] =
          merkle[static_cast<std::size_t>(i + 1)].Open(
              oracles.pi[static_cast<std::size_t>(i + 1)], mu_i + n_i);
      qp.folded[static_cast<std::size_t>(i)] =
          merkle[static_cast<std::size_t>(i)].Open(
              oracles.pi[static_cast<std::size_t>(i)], mu_i);
    }

    proof.query_proofs[static_cast<std::size_t>(q)] = std::move(qp);
  }

  return proof;
}

bool BaseFoldPCSVerifyEval(const MerkleRoot &commitment_C,
                           const std::vector<FieldElement> &z,
                           const FieldElement &claimed_y,
                           long num_queries,
                           const BaseFoldPCSEvalProof &proof,
                           const FoldableCodeParams &params) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->pcs_verify_ns : nullptr,
                    prof ? &prof->pcs_verify_calls : nullptr);

  ValidateParamsOrThrow(params);
  if (params.k0 != 1)
    return false;
  if (static_cast<long>(z.size()) != params.d)
    return false;
  if (num_queries < 0)
    return false;
  if (static_cast<long>(proof.commitments.roots_by_level.size()) != params.d + 1)
    return false;
  if (static_cast<long>(proof.h_by_level.size()) != params.d)
    return false;
  if (static_cast<long>(proof.query_proofs.size()) != num_queries)
    return false;

  if (proof.commitments.roots_by_level[static_cast<std::size_t>(params.d)] !=
      commitment_C) {
    return false;
  }

  const long n0 = CodewordLengthAtLevel(params, 0);
  if (proof.pi0_full.length() != n0)
    return false;

  Sha256Transcript transcript;
  AbsorbPublicInput(transcript, commitment_C, z, claimed_y);

  // h_d
  transcript.AbsorbQuadraticPoly(
      proof.h_by_level[static_cast<std::size_t>(params.d - 1)]);

  std::vector<FieldElement> r_by_level;
  r_by_level.resize(static_cast<std::size_t>(params.d));

  for (long i = params.d; i-- > 0;) {
    const FieldElement r_i =
        transcript.ChallengeFieldElement("r/" + std::to_string(i));
    r_by_level[static_cast<std::size_t>(i)] = r_i;

    transcript.AbsorbDigest(
        proof.commitments.roots_by_level[static_cast<std::size_t>(i)]);
    if (i > 0) {
      transcript.AbsorbQuadraticPoly(
          proof.h_by_level[static_cast<std::size_t>(i - 1)]);
    }
  }

  if (!CheckSumcheckRelations(proof.h_by_level, r_by_level, claimed_y))
    return false;

  const FieldElement r0 = r_by_level[0];
  const FieldElement h1_r0 = proof.h_by_level[0].Eval(r0);

  std::vector<FieldElement> r_point;
  r_point.resize(static_cast<std::size_t>(params.d));
  for (long i = 0; i < params.d; ++i) {
    r_point[static_cast<std::size_t>(i)] = r_by_level[static_cast<std::size_t>(i)];
  }

  const FieldElement eta = EqPolynomial(z, r_point);

  vec_ZZ_pE msg0;
  msg0.SetLength(1);
  msg0[0] = h1_r0;
  vec_ZZ_pE enc0;
  mul(enc0, msg0, params.G0);

  vec_ZZ_pE rhs = proof.pi0_full;
  for (long i = 0; i < rhs.length(); ++i) {
    rhs[i] *= eta;
  }
  if (enc0 != rhs)
    return false;

  IOPPChallenges challenges;
  challenges.alphas = r_by_level;

  if (params.d == 0) {
    return true;
  }

  const long n_last = CodewordLengthAtLevel(params, params.d - 1);
  for (long q = 0; q < num_queries; ++q) {
    const long mu =
        transcript.ChallengeIndex("mu/" + std::to_string(q), n_last);
    const IOPPQueryPlan plan = MakeQueryPlan(mu, params);

    const BaseFoldPCSQueryProof &qp =
        proof.query_proofs[static_cast<std::size_t>(q)];
    if (static_cast<long>(qp.left.size()) != params.d)
      return false;
    if (static_cast<long>(qp.right.size()) != params.d)
      return false;
    if (static_cast<long>(qp.folded.size()) != params.d)
      return false;

    IOPPQueryMerkleOpenings open;
    open.left = qp.left;
    open.right = qp.right;
    open.folded = qp.folded;
    open.pi0_full = proof.pi0_full;

    if (!VerifyQueryFromMerkleOpenings(plan, challenges, open, proof.commitments,
                                       params)) {
      return false;
    }
  }

  return true;
}

}  // namespace basefold
