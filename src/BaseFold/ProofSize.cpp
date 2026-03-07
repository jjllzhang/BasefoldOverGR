#include "BaseFold/ProofSize.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pX.h>

#include <cstddef>
#include <cstdint>

namespace basefold {
namespace {

std::size_t SerializedFieldElementSize(const FieldElement &x) {
  const long r = NTL::ZZ_pE::degree();
  if (r <= 0)
    NTL::LogicError("SerializedFieldElementSize: invalid extension degree");

  const NTL::ZZ_pX &poly = NTL::rep(x);
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
  static constexpr std::size_t kHashBytes = 32;
  std::size_t total = 0;
  total += 8;  // index as u64
  total += SerializedFieldElementSize(o.value);
  total += o.auth_path.sibling_hashes.size() * kHashBytes;
  return total;
}

std::size_t SerializedMerkleMultiproofSize(const MerkleMultiproof &proof) {
  static constexpr std::size_t kHashBytes = 32;
  std::size_t total = 0;
  total += proof.queried_indices.size() * 8;  // indices as u64
  for (long i = 0; i < proof.values.length(); ++i) {
    total += SerializedFieldElementSize(proof.values[i]);
  }
  total += proof.sibling_hashes.size() * kHashBytes;
  return total;
}

std::size_t SerializedExtensionElementSize(const NTL::ZZ_pEX &x) {
  const long d = NTL::deg(x);
  std::size_t total = 8;  // coeff count as u64
  for (long i = 0; i <= d; ++i) {
    const std::size_t c_size =
        SerializedFieldElementSize(NTL::coeff(x, i));
    total += 8 + c_size;
  }
  return total;
}

std::size_t SerializedExtensionMerkleOpeningSize(
    const ExtensionMerkleOpening &o) {
  static constexpr std::size_t kHashBytes = 32;
  std::size_t total = 0;
  total += 8;  // index as u64
  total += SerializedExtensionElementSize(o.value);
  total += o.auth_path.sibling_hashes.size() * kHashBytes;
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

  for (const MerkleMultiproof &proof : p.query_multiproofs) {
    total += static_cast<std::uint64_t>(SerializedMerkleMultiproofSize(proof));
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

  if (p.extension.enabled) {
    total += 1;  // enabled flag
    for (const MerkleRoot &r : p.extension.roots_by_level) {
      total += static_cast<std::uint64_t>(r.size());
    }
    for (const ExtensionQuadraticPoly &h : p.extension.h_by_level) {
      total += static_cast<std::uint64_t>(SerializedExtensionElementSize(h.a0));
      total += static_cast<std::uint64_t>(SerializedExtensionElementSize(h.a1));
      total += static_cast<std::uint64_t>(SerializedExtensionElementSize(h.a2));
    }
    for (const NTL::ZZ_pEX &r : p.extension.r_by_level) {
      total += static_cast<std::uint64_t>(SerializedExtensionElementSize(r));
    }
    for (const NTL::ZZ_pEX &c : p.extension.msg0_coeffs) {
      total += static_cast<std::uint64_t>(SerializedExtensionElementSize(c));
    }
    for (const NTL::ZZ_pEX &v : p.extension.pi0_full) {
      total += static_cast<std::uint64_t>(SerializedExtensionElementSize(v));
    }
    for (const BaseFoldPCSQueryProofExtension &qp : p.extension.query_proofs) {
      for (const ExtensionMerkleOpening &o : qp.left) {
        total += static_cast<std::uint64_t>(
            SerializedExtensionMerkleOpeningSize(o));
      }
      for (const ExtensionMerkleOpening &o : qp.right) {
        total += static_cast<std::uint64_t>(
            SerializedExtensionMerkleOpeningSize(o));
      }
      for (const ExtensionMerkleOpening &o : qp.folded) {
        total += static_cast<std::uint64_t>(
            SerializedExtensionMerkleOpeningSize(o));
      }
    }
  }

  return total;
}

double BaseFoldPCSEvalProofSizeKB(const BaseFoldPCSEvalProof &proof) {
  const std::uint64_t bytes = BaseFoldPCSEvalProofSizeBytes(proof);
  return static_cast<double>(bytes) / 1024.0;
}

}  // namespace basefold
