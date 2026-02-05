#include "BaseFold/ProofSize.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <cstddef>
#include <cstdint>

namespace basefold {
namespace {

std::size_t SerializedFieldElementSize(const FieldElement &x) {
  const long r = NTL::ZZ_pE::degree();
  if (r <= 0)
    NTL::LogicError("SerializedFieldElementSize: invalid extension degree");

  const NTL::ZZ_pX poly = NTL::rep(x);
  std::size_t total = 8;  // degree r as u64
  for (long i = 0; i < r; ++i) {
    const NTL::ZZ c = NTL::rep(NTL::coeff(poly, i));
    const long n = NTL::NumBytes(c);
    if (n < 0)
      NTL::LogicError("SerializedFieldElementSize: NumBytes returned < 0");
    total += 8 + static_cast<std::size_t>(n);
  }
  return total;
}

std::size_t SerializedMerkleOpeningSize(const MerkleOpening &o) {
  std::size_t total = 0;
  total += 8;  // index as u64
  total += SerializedFieldElementSize(o.value);

  if (o.auth_path.sibling_hashes.size() != o.auth_path.sibling_is_left.size()) {
    NTL::LogicError("SerializedMerkleOpeningSize: auth_path size mismatch");
  }

  for (const Bytes &h : o.auth_path.sibling_hashes) {
    total += h.size();  // hash bytes (sha256 in current implementation)
  }
  total += o.auth_path.sibling_is_left.size();  // 1 byte per direction bit
  return total;
}

}  // namespace

std::uint64_t BaseFoldPCSEvalProofSizeBytes(const BaseFoldPCSEvalProof &p) {
  std::uint64_t total = 0;

  for (const MerkleRoot &r : p.commitments.roots_by_level) {
    total += static_cast<std::uint64_t>(r.size());
  }

  for (const QuadraticPoly &h : p.h_by_level) {
    total += static_cast<std::uint64_t>(SerializedFieldElementSize(h.a0));
    total += static_cast<std::uint64_t>(SerializedFieldElementSize(h.a1));
    total += static_cast<std::uint64_t>(SerializedFieldElementSize(h.a2));
  }

  for (long i = 0; i < p.pi0_full.length(); ++i) {
    total += static_cast<std::uint64_t>(SerializedFieldElementSize(p.pi0_full[i]));
  }

  for (const BaseFoldPCSQueryProof &qp : p.query_proofs) {
    for (const MerkleOpening &o : qp.left) {
      total += static_cast<std::uint64_t>(SerializedMerkleOpeningSize(o));
    }
    for (const MerkleOpening &o : qp.right) {
      total += static_cast<std::uint64_t>(SerializedMerkleOpeningSize(o));
    }
    for (const MerkleOpening &o : qp.folded) {
      total += static_cast<std::uint64_t>(SerializedMerkleOpeningSize(o));
    }
  }

  return total;
}

double BaseFoldPCSEvalProofSizeKB(const BaseFoldPCSEvalProof &proof) {
  const std::uint64_t bytes = BaseFoldPCSEvalProofSizeBytes(proof);
  return static_cast<double>(bytes) / 1024.0;
}

}  // namespace basefold

