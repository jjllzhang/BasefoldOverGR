#include "BaseFold/Hash.hpp"

#include <NTL/ZZ.h>

#include <algorithm>
#include <cstddef>
#include <string>

#include "blake3.h"

using NTL::LogicError;

namespace basefold {
namespace {

static_assert(BLAKE3_OUT_LEN == 32, "BLAKE3 output size must stay 32 bytes");

Bytes ComputeHashBytes(const Byte *data, std::size_t len, const char *func_name) {
  Bytes out(Digest{}.size(), 0);
  (void)func_name;
  blake3_hasher hasher;
  blake3_hasher_init(&hasher);
  if (len > 0) {
    blake3_hasher_update(&hasher, data, len);
  }
  blake3_hasher_finalize(&hasher, out.data(), out.size());
  return out;
}

}  // namespace

std::string SelectedHashBackendName() { return "blake3"; }

std::size_t DigestBytes() { return Digest{}.size(); }

Bytes HashBytes(const Byte *data, std::size_t len) {
  return ComputeHashBytes(data, len, "HashBytes");
}

Digest HashDigest(const Byte *data, std::size_t len, const char *func_name) {
  const Bytes hashed = ComputeHashBytes(data, len, func_name);
  if (hashed.size() != Digest{}.size()) {
    const std::string msg = std::string(func_name) + ": unexpected digest size";
    LogicError(msg.c_str());
  }
  Digest out{};
  std::copy(hashed.begin(), hashed.end(), out.begin());
  return out;
}

}  // namespace basefold
