#ifndef BASEFOLD_HASH_HPP_
#define BASEFOLD_HASH_HPP_

#include <cstddef>
#include <string>

#include "BaseFold/IOPP.hpp"

namespace basefold {

std::string SelectedHashBackendName();

std::size_t DigestBytes();
Bytes HashBytes(const Byte *data, std::size_t len);

inline Bytes HashBytes(const Bytes &data) {
  return HashBytes(data.data(), data.size());
}

Digest HashDigest(const Byte *data, std::size_t len, const char *func_name);

}  // namespace basefold

#endif  // BASEFOLD_HASH_HPP_
