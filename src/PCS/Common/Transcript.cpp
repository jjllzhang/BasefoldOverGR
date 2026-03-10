#include "PCS/Common/Transcript.hpp"

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>

#include "PCS/Common/Profile.hpp"

using NTL::BytesFromZZ;
using NTL::LogicError;
using NTL::NumBits;
using NTL::NumBytes;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pX;
using NTL::ZZFromBytes;

namespace basefold {
namespace {

void AppendU64(Bytes &out, std::uint64_t value, TranscriptByteOrder byte_order) {
  if (byte_order == TranscriptByteOrder::kBigEndian) {
    for (int i = 7; i >= 0; --i) {
      out.push_back(static_cast<Byte>((value >> (8 * i)) & 0xff));
    }
    return;
  }
  for (int i = 0; i < 8; ++i) {
    out.push_back(static_cast<Byte>((value >> (8 * i)) & 0xff));
  }
}

void AppendSerializedFieldElement(Bytes &out, const FieldElement &x,
                                  const char *func_name,
                                  TranscriptByteOrder byte_order) {
  const long r = ZZ_pE::degree();
  if (r <= 0) {
    const std::string msg =
        std::string(func_name) + ": invalid extension degree";
    LogicError(msg.c_str());
  }

  const ZZ_pX &poly = NTL::rep(x);
  AppendU64(out, static_cast<std::uint64_t>(r), byte_order);
  for (long i = 0; i < r; ++i) {
    const ZZ c = NTL::rep(NTL::coeff(poly, i));
    const long n = NumBytes(c);
    AppendU64(out, static_cast<std::uint64_t>(n), byte_order);
    if (n > 0) {
      const std::size_t old_size = out.size();
      out.resize(old_size + static_cast<std::size_t>(n));
      BytesFromZZ(reinterpret_cast<unsigned char *>(out.data() + old_size), c, n);
    }
  }
}

Bytes SerializeFieldElement(const FieldElement &x,
                            TranscriptByteOrder byte_order) {
  Bytes out;
  AppendSerializedFieldElement(out, x, "SerializeFieldElement", byte_order);
  return out;
}

Bytes TaggedHash(Byte tag, const Bytes &state, const Bytes &payload,
                 TranscriptByteOrder byte_order) {
  Bytes in;
  in.reserve(1 + state.size() + 8 + payload.size());
  in.push_back(tag);
  in.insert(in.end(), state.begin(), state.end());
  AppendU64(in, static_cast<std::uint64_t>(payload.size()), byte_order);
  in.insert(in.end(), payload.begin(), payload.end());
  return HashBytes(in);
}

Bytes TaggedHash(Byte tag, const Bytes &state, const std::string &payload,
                 TranscriptByteOrder byte_order) {
  Bytes bytes;
  bytes.reserve(payload.size());
  bytes.insert(bytes.end(), payload.begin(), payload.end());
  return TaggedHash(tag, state, bytes, byte_order);
}

class ChallengeStream {
 public:
  ChallengeStream(const HashTranscriptConfig &config, const Bytes &state,
                  const std::string &label)
      : config_(config), state_(state), label_(label) {}

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
    if (upper_bound <= 0) {
      LogicError((config_.error_prefix + "::SampleZZLessThan: upper_bound must "
                  "be positive")
                     .c_str());
    }
    if (upper_bound == 1) {
      return ZZ(0);
    }

    const ZZ ub_minus_1 = upper_bound - 1;
    const long bits = NumBits(ub_minus_1);
    if (bits <= 0) {
      return ZZ(0);
    }
    const long byte_len = (bits + 7) / 8;
    const ZZ two_to_bits = ZZ(1) << bits;

    Bytes tmp(static_cast<std::size_t>(byte_len));
    while (true) {
      ReadBytes(reinterpret_cast<std::uint8_t *>(tmp.data()),
                static_cast<std::size_t>(byte_len));
      ZZ x =
          ZZFromBytes(reinterpret_cast<const unsigned char *>(tmp.data()), byte_len);
      x %= two_to_bits;
      if (x < upper_bound) {
        return x;
      }
    }
  }

 private:
  Bytes Digest(std::uint64_t ctr) const {
    Bytes payload;
    payload.reserve(8);
    AppendU64(payload, ctr, config_.byte_order);
    const Bytes digest_state =
        TaggedHash(static_cast<Byte>(0x20), state_, label_, config_.byte_order);
    return TaggedHash(static_cast<Byte>(0x21), digest_state, payload,
                      config_.byte_order);
  }

  const HashTranscriptConfig &config_;
  Bytes state_;
  std::string label_;
  std::uint64_t counter_ = 0;
  Bytes buf_;
  std::size_t offset_ = 0;
};

}  // namespace

HashTranscript::HashTranscript(HashTranscriptConfig config)
    : config_(std::move(config)) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->transcript_absorb_ns : nullptr,
                    prof ? &prof->transcript_absorb_calls : nullptr);
  state_ = TaggedHash(static_cast<Byte>(0x42), Bytes{}, config_.domain_separator,
                      config_.byte_order);
}

void HashTranscript::AbsorbBytes(const Byte *data, std::size_t len) {
  Bytes bytes;
  bytes.resize(len);
  if (len > 0) {
    std::memcpy(bytes.data(), data, len);
  }
  AbsorbBytes(bytes);
}

void HashTranscript::AbsorbBytes(const Bytes &data) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->transcript_absorb_ns : nullptr,
                    prof ? &prof->transcript_absorb_calls : nullptr);
  state_ = TaggedHash(static_cast<Byte>(0x01), state_, data, config_.byte_order);
}

void HashTranscript::AbsorbDigest(const Digest &digest) {
  const Bytes bytes(digest.begin(), digest.end());
  AbsorbBytes(bytes);
}

void HashTranscript::AbsorbFieldElement(const FieldElement &x) {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->transcript_absorb_ns : nullptr,
                    prof ? &prof->transcript_absorb_calls : nullptr);
  state_ = TaggedHash(static_cast<Byte>(0x02), state_,
                      SerializeFieldElement(x, config_.byte_order),
                      config_.byte_order);
}

FieldElement HashTranscript::ChallengeFieldElement(const std::string &label) {
  return static_cast<const HashTranscript &>(*this).ChallengeFieldElement(label);
}

FieldElement HashTranscript::ChallengeFieldElement(
    const std::string &label) const {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->transcript_challenge_ns : nullptr,
                    prof ? &prof->transcript_challenge_calls : nullptr);

  const long r = ZZ_pE::degree();
  if (r <= 0) {
    LogicError((config_.error_prefix +
                "::ChallengeFieldElement: invalid extension degree")
                   .c_str());
  }

  const ZZ modulus = ZZ_p::modulus();
  if (modulus <= 1) {
    LogicError((config_.error_prefix +
                "::ChallengeFieldElement: invalid base modulus")
                   .c_str());
  }

  ChallengeStream stream(config_, state_, "fe/" + label);
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

long HashTranscript::ChallengeIndex(const std::string &label, long upper_bound) {
  return static_cast<const HashTranscript &>(*this).ChallengeIndex(label,
                                                                  upper_bound);
}

long HashTranscript::ChallengeIndex(const std::string &label,
                                    long upper_bound) const {
  Profile *prof = ActiveProfile();
  ScopedTimer timer(prof ? &prof->transcript_challenge_ns : nullptr,
                    prof ? &prof->transcript_challenge_calls : nullptr);

  if (upper_bound <= 0) {
    LogicError((config_.error_prefix + "::ChallengeIndex: upper_bound must be "
                "positive")
                   .c_str());
  }
  if (upper_bound == 1) {
    return 0;
  }

  ChallengeStream stream(config_, state_, "idx/" + label);
  std::uint64_t ub = static_cast<std::uint64_t>(upper_bound);
  std::uint64_t t = ub - 1;
  int bits = 0;
  while (t > 0) {
    ++bits;
    t >>= 1;
  }
  if (bits <= 0 || bits > 63) {
    LogicError(
        (config_.error_prefix + "::ChallengeIndex: unsupported range").c_str());
  }
  const int byte_len = (bits + 7) / 8;
  const std::uint64_t mask = (bits == 64) ? ~0ULL : ((1ULL << bits) - 1ULL);

  while (true) {
    std::uint8_t buf[8];
    std::memset(buf, 0, sizeof(buf));
    stream.ReadBytes(buf, static_cast<std::size_t>(byte_len));
    std::uint64_t x = 0;
    for (int i = 0; i < byte_len; ++i) {
      x = (x << 8) | static_cast<std::uint64_t>(buf[i]);
    }
    x &= mask;
    if (x < ub) {
      return static_cast<long>(x);
    }
  }
}

void AbsorbQuadraticPoly(HashTranscript &transcript, const QuadraticPoly &poly) {
  transcript.AbsorbFieldElement(poly.a0);
  transcript.AbsorbFieldElement(poly.a1);
  transcript.AbsorbFieldElement(poly.a2);
}

}  // namespace basefold
