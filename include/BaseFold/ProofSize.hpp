#ifndef BASEFOLD_PROOFSIZE_HPP_
#define BASEFOLD_PROOFSIZE_HPP_

#include <cstdint>

#include "BaseFold/BaseFoldPCS.hpp"

namespace basefold {

// Estimates the serialized size (bytes) of a BaseFold PCS evaluation proof.
//
// This is a bench-oriented estimate consistent with the Merkle+FS data actually
// carried by BaseFoldPCSEvalProof:
// - Merkle roots for π_0..π_d
// - Sumcheck messages h_i
// - Full π_0
// - Merkle openings for all queries and levels
std::uint64_t BaseFoldPCSEvalProofSizeBytes(const BaseFoldPCSEvalProof &proof);

// Convenience wrapper returning size in KiB (1024 bytes).
double BaseFoldPCSEvalProofSizeKB(const BaseFoldPCSEvalProof &proof);

}  // namespace basefold

#endif  // BASEFOLD_PROOFSIZE_HPP_

