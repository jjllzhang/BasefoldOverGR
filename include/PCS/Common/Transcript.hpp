#ifndef BASEFOLD_TRANSCRIPT_HPP_
#define BASEFOLD_TRANSCRIPT_HPP_

#include <cstddef>
#include <string>

#include "PCS/Common/Merkle.hpp"
#include "PCS/Common/Sumcheck.hpp"

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

enum class TranscriptByteOrder {
  kBigEndian = 0,
  kLittleEndian = 1,
};

struct HashTranscriptConfig {
  std::string domain_separator;
  TranscriptByteOrder byte_order = TranscriptByteOrder::kBigEndian;
  std::string error_prefix = "HashTranscript";
};

class HashTranscript : public FiatShamirTranscript {
 public:
  explicit HashTranscript(HashTranscriptConfig config);

  void AbsorbBytes(const Byte *data, std::size_t len) override;
  void AbsorbBytes(const Bytes &data) override;
  void AbsorbDigest(const Digest &digest);
  void AbsorbFieldElement(const FieldElement &x) override;

  FieldElement ChallengeFieldElement(const std::string &label) override;
  FieldElement ChallengeFieldElement(const std::string &label) const;

  long ChallengeIndex(const std::string &label, long upper_bound) override;
  long ChallengeIndex(const std::string &label, long upper_bound) const;

 private:
  HashTranscriptConfig config_;
  Bytes state_;
};

void AbsorbQuadraticPoly(HashTranscript &transcript, const QuadraticPoly &poly);

}  // namespace basefold

#endif  // BASEFOLD_TRANSCRIPT_HPP_
