#ifndef BASEFOLD_TRANSCRIPT_HPP_
#define BASEFOLD_TRANSCRIPT_HPP_

#include <cstddef>
#include <string>

#include "PCS/Common/Merkle.hpp"

namespace basefold {

class FiatShamirTranscript {
 public:
  virtual ~FiatShamirTranscript() = default;

  virtual void AbsorbBytes(const Byte *data, std::size_t len) = 0;
  virtual void AbsorbBytes(const Bytes &data) = 0;
  virtual void AbsorbFieldElement(const FieldElement &x) = 0;

  virtual FieldElement ChallengeFieldElement(const std::string &label) = 0;
  virtual long ChallengeIndex(const std::string &label, long upper_bound) = 0;
};

}  // namespace basefold

#endif  // BASEFOLD_TRANSCRIPT_HPP_
