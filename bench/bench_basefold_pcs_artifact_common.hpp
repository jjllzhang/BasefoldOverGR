#ifndef BASEFOLD_BENCH_PCS_ARTIFACT_COMMON_HPP_
#define BASEFOLD_BENCH_PCS_ARTIFACT_COMMON_HPP_

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <NTL/ZZ_pX.h>

#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "PCS/BaseFold/BaseFoldPCS.hpp"
#include "PCS/BaseFold/ProofDeserialize.hpp"
#include "PCS/BaseFold/ProofSerialize.hpp"
#include "PCS/Common/Hash.hpp"
#include "PCS/Common/Merkle.hpp"
#include "bench_basefold_pcs_common.hpp"

namespace basefold_bench_pcs_artifact {

using NTL::SetCoeff;
using NTL::ZZ;
using NTL::ZZ_p;
using NTL::ZZ_pE;
using NTL::ZZ_pEX;
using NTL::ZZ_pPush;
using NTL::ZZ_pX;
using NTL::conv;
using NTL::deg;

namespace fs = std::filesystem;

struct ArtifactMetadata {
  std::string artifact_id;
  std::string display_key;
  std::string context_id;
  std::string context_label;
  std::string mode;
  long c = 0;
  long k0 = 0;
  long d = 0;
  long poly_dim = 0;
  std::string lambda;
  std::string gamma;
  long queries = 0;
  std::uint64_t seed = 0;
  bool use_checked_prover_path = false;
  bool use_extension_challenges = false;
  ZZ scalar_modulus = ZZ(0);
  ZZ base_prime = ZZ(0);
  std::vector<ZZ> F_coeffs;
  std::vector<ZZ> zeta_coeffs;
  std::string zeta_source;
  std::vector<std::vector<ZZ>> challenge_extension_coeffs;
  std::string hash_backend;
  std::string proof_encoding;
  std::uint64_t proof_size_bytes = 0;
};

struct ArtifactManifestEntry {
  std::string artifact_id;
  std::string display_key;
  std::string context_id;
  std::string context_label;
  std::string mode;
  long c = 0;
  long k0 = 0;
  long d = 0;
  long poly_dim = 0;
  std::string lambda;
  std::string gamma;
  long queries = 0;
  std::uint64_t seed = 0;
  bool use_extension_challenges = false;
  std::string object_relpath;
};

struct ArtifactPublicInputs {
  basefold::MerkleRoot commitment_root{};
  std::vector<basefold::FieldElement> z;
  basefold::FieldElement claimed_y;
};

struct RestoredVerificationContext {
  basefold::FoldableCodeParams params;
  basefold::BaseFoldPCSChallengeConfig challenge_cfg;
};

inline RestoredVerificationContext RestoreVerificationContext(
    const ArtifactMetadata &meta);

inline std::string NormalizeModeLabel(const std::string &mode) {
  if (mode == "field" || mode == "ring") {
    return mode;
  }
  throw std::runtime_error("NormalizeModeLabel: mode must be field or ring");
}

inline std::string EscapeJson(const std::string &s) {
  std::ostringstream out;
  for (unsigned char ch : s) {
    switch (ch) {
      case '\\':
        out << "\\\\";
        break;
      case '"':
        out << "\\\"";
        break;
      case '\b':
        out << "\\b";
        break;
      case '\f':
        out << "\\f";
        break;
      case '\n':
        out << "\\n";
        break;
      case '\r':
        out << "\\r";
        break;
      case '\t':
        out << "\\t";
        break;
      default:
        if (ch < 0x20) {
          out << "\\u" << std::hex << std::setw(4) << std::setfill('0')
              << static_cast<unsigned int>(ch) << std::dec;
        } else {
          out << static_cast<char>(ch);
        }
    }
  }
  return out.str();
}

inline std::string QuoteJson(const std::string &s) {
  return "\"" + EscapeJson(s) + "\"";
}

inline std::string ZZToString(const ZZ &x) {
  std::ostringstream out;
  out << x;
  return out.str();
}

inline std::string JoinStrings(const std::vector<std::string> &parts,
                               const std::string &sep) {
  std::ostringstream out;
  for (std::size_t i = 0; i < parts.size(); ++i) {
    if (i != 0) {
      out << sep;
    }
    out << parts[i];
  }
  return out.str();
}

inline std::string SerializeZZList(const std::vector<ZZ> &xs) {
  std::vector<std::string> parts;
  parts.reserve(xs.size());
  for (const ZZ &x : xs) {
    parts.push_back(QuoteJson(ZZToString(x)));
  }
  return "[" + JoinStrings(parts, ",") + "]";
}

inline std::string SerializeNestedZZList(
    const std::vector<std::vector<ZZ>> &xss) {
  std::vector<std::string> parts;
  parts.reserve(xss.size());
  for (const std::vector<ZZ> &xs : xss) {
    parts.push_back(SerializeZZList(xs));
  }
  return "[" + JoinStrings(parts, ",") + "]";
}

inline void AppendJsonKV(std::ostringstream &out, bool &first,
                         const std::string &key, const std::string &value) {
  if (!first) {
    out << ",\n";
  }
  first = false;
  out << "  " << QuoteJson(key) << ": " << value;
}

inline std::string ToCanonicalArtifactIdInput(const ArtifactMetadata &meta) {
  std::ostringstream out;
  out << "context_id=" << meta.context_id << '\n';
  out << "context_label=" << meta.context_label << '\n';
  out << "mode=" << NormalizeModeLabel(meta.mode) << '\n';
  out << "c=" << meta.c << '\n';
  out << "k0=" << meta.k0 << '\n';
  out << "d=" << meta.d << '\n';
  out << "poly_dim=" << meta.poly_dim << '\n';
  out << "lambda=" << meta.lambda << '\n';
  out << "gamma=" << meta.gamma << '\n';
  out << "queries=" << meta.queries << '\n';
  out << "seed=" << meta.seed << '\n';
  out << "use_extension_challenges="
      << (meta.use_extension_challenges ? 1 : 0) << '\n';
  return out.str();
}

inline std::string HexDigestPrefix(const basefold::Digest &digest,
                                   std::size_t hex_len) {
  std::ostringstream out;
  out << std::hex << std::setfill('0');
  for (basefold::Byte b : digest) {
    out << std::setw(2) << static_cast<unsigned int>(b);
    if (static_cast<std::size_t>(out.tellp()) >= hex_len) {
      break;
    }
  }
  std::string s = out.str();
  if (s.size() > hex_len) {
    s.resize(hex_len);
  }
  return s;
}

inline std::string ComputeCanonicalArtifactId(const ArtifactMetadata &meta) {
  const std::string canonical = ToCanonicalArtifactIdInput(meta);
  const basefold::Digest digest =
      basefold::HashDigest(reinterpret_cast<const basefold::Byte *>(
                               canonical.data()),
                           canonical.size(),
                           "ComputeCanonicalArtifactId");
  return "bfv_" + HexDigestPrefix(digest, 16);
}

inline std::string BuildArtifactDisplayKey(const ArtifactMetadata &meta) {
  std::ostringstream out;
  out << "ctx=" << meta.context_id << ",label=" << meta.context_label
      << ",mode=" << NormalizeModeLabel(meta.mode) << ",d=" << meta.d
      << ",poly_dim=" << meta.poly_dim << ",c=" << meta.c
      << ",k0=" << meta.k0 << ",lambda=" << meta.lambda
      << ",gamma=" << meta.gamma << ",queries=" << meta.queries
      << ",seed=" << meta.seed;
  if (meta.use_extension_challenges) {
    out << ",ext=on";
  }
  return out.str();
}

inline long ComputePolyDimOrThrow(long k0, long d) {
  if (k0 <= 0 || d < 0) {
    throw std::runtime_error("ComputePolyDimOrThrow: invalid k0 or d");
  }
  const long pow2 = basefold_bench_pcs_common::Pow2Checked(d);
  if (k0 > std::numeric_limits<long>::max() / pow2) {
    throw std::runtime_error("ComputePolyDimOrThrow: overflow");
  }
  return k0 * pow2;
}

inline std::string ManifestObjectRelPath(const std::string &artifact_id) {
  return std::string("objects/") + artifact_id;
}

inline void ValidateArtifactIdOrThrow(const std::string &artifact_id) {
  if (artifact_id.empty()) {
    throw std::runtime_error("ValidateArtifactIdOrThrow: artifact_id is empty");
  }
  if (artifact_id == "." || artifact_id == "..") {
    throw std::runtime_error(
        "ValidateArtifactIdOrThrow: artifact_id must not be . or ..");
  }
  for (unsigned char ch : artifact_id) {
    if (std::isalnum(ch) || ch == '_' || ch == '-' || ch == '.') {
      continue;
    }
    throw std::runtime_error(
        "ValidateArtifactIdOrThrow: artifact_id contains unsafe characters");
  }
}

inline fs::path ArtifactObjectDir(const fs::path &root,
                                  const std::string &artifact_id) {
  return root / "objects" / artifact_id;
}

inline fs::path ArtifactManifestPath(const fs::path &root) {
  return root / "manifest.jsonl";
}

inline fs::path ArtifactMetadataPath(const fs::path &root,
                                     const std::string &artifact_id) {
  return ArtifactObjectDir(root, artifact_id) / "meta.json";
}

inline fs::path ArtifactPublicInputsPath(const fs::path &root,
                                         const std::string &artifact_id) {
  return ArtifactObjectDir(root, artifact_id) / "public_inputs.bin";
}

inline fs::path ArtifactProofPath(const fs::path &root,
                                  const std::string &artifact_id) {
  return ArtifactObjectDir(root, artifact_id) / "proof.bin";
}

inline std::string SerializeMetadataJson(const ArtifactMetadata &meta) {
  std::ostringstream out;
  out << "{\n";
  bool first = true;
  AppendJsonKV(out, first, "artifact_id", QuoteJson(meta.artifact_id));
  AppendJsonKV(out, first, "display_key", QuoteJson(meta.display_key));
  AppendJsonKV(out, first, "context_id", QuoteJson(meta.context_id));
  AppendJsonKV(out, first, "context_label", QuoteJson(meta.context_label));
  AppendJsonKV(out, first, "mode", QuoteJson(meta.mode));
  AppendJsonKV(out, first, "c", std::to_string(meta.c));
  AppendJsonKV(out, first, "k0", std::to_string(meta.k0));
  AppendJsonKV(out, first, "d", std::to_string(meta.d));
  AppendJsonKV(out, first, "poly_dim", std::to_string(meta.poly_dim));
  AppendJsonKV(out, first, "lambda", QuoteJson(meta.lambda));
  AppendJsonKV(out, first, "gamma", QuoteJson(meta.gamma));
  AppendJsonKV(out, first, "queries", std::to_string(meta.queries));
  AppendJsonKV(out, first, "seed", std::to_string(meta.seed));
  AppendJsonKV(out, first, "use_checked_prover_path",
               meta.use_checked_prover_path ? "true" : "false");
  AppendJsonKV(out, first, "use_extension_challenges",
               meta.use_extension_challenges ? "true" : "false");
  AppendJsonKV(out, first, "scalar_modulus", QuoteJson(ZZToString(meta.scalar_modulus)));
  AppendJsonKV(out, first, "base_prime", QuoteJson(ZZToString(meta.base_prime)));
  AppendJsonKV(out, first, "F_coeffs", SerializeZZList(meta.F_coeffs));
  AppendJsonKV(out, first, "zeta_coeffs", SerializeZZList(meta.zeta_coeffs));
  AppendJsonKV(out, first, "zeta_source", QuoteJson(meta.zeta_source));
  AppendJsonKV(out, first, "challenge_extension_coeffs",
               SerializeNestedZZList(meta.challenge_extension_coeffs));
  AppendJsonKV(out, first, "hash_backend", QuoteJson(meta.hash_backend));
  AppendJsonKV(out, first, "proof_encoding", QuoteJson(meta.proof_encoding));
  AppendJsonKV(out, first, "proof_size_bytes",
               std::to_string(meta.proof_size_bytes));
  out << "\n}\n";
  return out.str();
}

inline ArtifactManifestEntry ManifestEntryFromMetadata(
    const ArtifactMetadata &meta) {
  ArtifactManifestEntry entry;
  entry.artifact_id = meta.artifact_id;
  entry.display_key = meta.display_key;
  entry.context_id = meta.context_id;
  entry.context_label = meta.context_label;
  entry.mode = meta.mode;
  entry.c = meta.c;
  entry.k0 = meta.k0;
  entry.d = meta.d;
  entry.poly_dim = meta.poly_dim;
  entry.lambda = meta.lambda;
  entry.gamma = meta.gamma;
  entry.queries = meta.queries;
  entry.seed = meta.seed;
  entry.use_extension_challenges = meta.use_extension_challenges;
  entry.object_relpath = ManifestObjectRelPath(meta.artifact_id);
  return entry;
}

inline std::string SerializeManifestEntryJsonLine(
    const ArtifactManifestEntry &entry) {
  std::ostringstream out;
  out << "{";
  out << QuoteJson("artifact_id") << ":" << QuoteJson(entry.artifact_id)
      << ",";
  out << QuoteJson("display_key") << ":" << QuoteJson(entry.display_key)
      << ",";
  out << QuoteJson("context_id") << ":" << QuoteJson(entry.context_id)
      << ",";
  out << QuoteJson("context_label") << ":" << QuoteJson(entry.context_label)
      << ",";
  out << QuoteJson("mode") << ":" << QuoteJson(entry.mode) << ",";
  out << QuoteJson("c") << ":" << entry.c << ",";
  out << QuoteJson("k0") << ":" << entry.k0 << ",";
  out << QuoteJson("d") << ":" << entry.d << ",";
  out << QuoteJson("poly_dim") << ":" << entry.poly_dim << ",";
  out << QuoteJson("lambda") << ":" << QuoteJson(entry.lambda) << ",";
  out << QuoteJson("gamma") << ":" << QuoteJson(entry.gamma) << ",";
  out << QuoteJson("queries") << ":" << entry.queries << ",";
  out << QuoteJson("seed") << ":" << entry.seed << ",";
  out << QuoteJson("use_extension_challenges") << ":"
      << (entry.use_extension_challenges ? "true" : "false") << ",";
  out << QuoteJson("object_relpath") << ":"
      << QuoteJson(entry.object_relpath) << "}";
  return out.str();
}

inline std::string ReadFileToString(const fs::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("ReadFileToString: failed to open " +
                             path.string());
  }
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

inline void WriteStringToFile(const fs::path &path, const std::string &text) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("WriteStringToFile: failed to open " +
                             path.string());
  }
  out << text;
  if (!out) {
    throw std::runtime_error("WriteStringToFile: failed to write " +
                             path.string());
  }
}

inline void WriteBytesToFile(const fs::path &path, const basefold::Bytes &bytes) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("WriteBytesToFile: failed to open " +
                             path.string());
  }
  if (!bytes.empty()) {
    out.write(reinterpret_cast<const char *>(bytes.data()),
              static_cast<std::streamsize>(bytes.size()));
  }
  if (!out) {
    throw std::runtime_error("WriteBytesToFile: failed to write " +
                             path.string());
  }
}

inline basefold::Bytes ReadBytesFromFile(const fs::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("ReadBytesFromFile: failed to open " +
                             path.string());
  }
  std::vector<char> raw((std::istreambuf_iterator<char>(in)),
                        std::istreambuf_iterator<char>());
  basefold::Bytes bytes;
  bytes.reserve(raw.size());
  for (char ch : raw) {
    bytes.push_back(
        static_cast<basefold::Byte>(static_cast<unsigned char>(ch)));
  }
  return bytes;
}

inline std::string TrimAscii(const std::string &s) {
  std::size_t start = 0;
  while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
    ++start;
  }
  std::size_t end = s.size();
  while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
    --end;
  }
  return s.substr(start, end - start);
}

inline std::string ParseJsonStringToken(const std::string &token) {
  if (token.size() < 2 || token.front() != '"' || token.back() != '"') {
    throw std::runtime_error("ParseJsonStringToken: expected JSON string");
  }
  std::ostringstream out;
  for (std::size_t i = 1; i + 1 < token.size(); ++i) {
    char ch = token[i];
    if (ch != '\\') {
      out << ch;
      continue;
    }
    if (i + 1 >= token.size() - 1) {
      throw std::runtime_error("ParseJsonStringToken: truncated escape");
    }
    ++i;
    const char esc = token[i];
    switch (esc) {
      case '\\': out << '\\'; break;
      case '"': out << '"'; break;
      case 'b': out << '\b'; break;
      case 'f': out << '\f'; break;
      case 'n': out << '\n'; break;
      case 'r': out << '\r'; break;
      case 't': out << '\t'; break;
      case 'u':
        throw std::runtime_error(
            "ParseJsonStringToken: unicode escapes are unsupported");
      default:
        throw std::runtime_error("ParseJsonStringToken: unsupported escape");
    }
  }
  return out.str();
}

inline std::string ExtractJsonValueToken(const std::string &json,
                                         const std::string &key) {
  const std::string needle = QuoteJson(key);
  const std::size_t key_pos = json.find(needle);
  if (key_pos == std::string::npos) {
    throw std::runtime_error("ExtractJsonValueToken: missing key " + key);
  }
  const std::size_t colon_pos = json.find(':', key_pos + needle.size());
  if (colon_pos == std::string::npos) {
    throw std::runtime_error("ExtractJsonValueToken: missing colon for key " + key);
  }
  std::size_t pos = colon_pos + 1;
  while (pos < json.size() && std::isspace(static_cast<unsigned char>(json[pos]))) {
    ++pos;
  }
  if (pos >= json.size()) {
    throw std::runtime_error("ExtractJsonValueToken: missing value for key " + key);
  }

  if (json[pos] == '"') {
    std::size_t end = pos + 1;
    bool escape = false;
    for (; end < json.size(); ++end) {
      const char ch = json[end];
      if (escape) {
        escape = false;
        continue;
      }
      if (ch == '\\') {
        escape = true;
        continue;
      }
      if (ch == '"') {
        return json.substr(pos, end - pos + 1);
      }
    }
    throw std::runtime_error("ExtractJsonValueToken: unterminated string for key " + key);
  }

  if (json[pos] == '[') {
    std::size_t end = pos;
    int depth = 0;
    bool in_string = false;
    bool escape = false;
    for (; end < json.size(); ++end) {
      const char ch = json[end];
      if (in_string) {
        if (escape) {
          escape = false;
        } else if (ch == '\\') {
          escape = true;
        } else if (ch == '"') {
          in_string = false;
        }
        continue;
      }
      if (ch == '"') {
        in_string = true;
      } else if (ch == '[') {
        ++depth;
      } else if (ch == ']') {
        --depth;
        if (depth == 0) {
          return json.substr(pos, end - pos + 1);
        }
      }
    }
    throw std::runtime_error("ExtractJsonValueToken: unterminated array for key " + key);
  }

  std::size_t end = pos;
  while (end < json.size() && json[end] != ',' && json[end] != '}' &&
         json[end] != '\n') {
    ++end;
  }
  return TrimAscii(json.substr(pos, end - pos));
}

inline std::string ParseJsonStringField(const std::string &json,
                                        const std::string &key) {
  return ParseJsonStringToken(ExtractJsonValueToken(json, key));
}

inline long ParseJsonLongField(const std::string &json, const std::string &key) {
  const std::string token = ExtractJsonValueToken(json, key);
  std::size_t idx = 0;
  const long value = std::stol(token, &idx, 10);
  if (idx != token.size()) {
    throw std::runtime_error("ParseJsonLongField: invalid integer for key " + key);
  }
  return value;
}

inline std::uint64_t ParseJsonU64Field(const std::string &json,
                                       const std::string &key) {
  const std::string token = ExtractJsonValueToken(json, key);
  std::size_t idx = 0;
  const std::uint64_t value = std::stoull(token, &idx, 10);
  if (idx != token.size()) {
    throw std::runtime_error("ParseJsonU64Field: invalid integer for key " + key);
  }
  return value;
}

inline bool ParseJsonBoolField(const std::string &json, const std::string &key) {
  const std::string token = ExtractJsonValueToken(json, key);
  if (token == "true") {
    return true;
  }
  if (token == "false") {
    return false;
  }
  throw std::runtime_error("ParseJsonBoolField: invalid bool for key " + key);
}

inline std::vector<ZZ> ParseJsonZZListField(const std::string &json,
                                            const std::string &key) {
  const std::string token = ExtractJsonValueToken(json, key);
  if (token == "[]") {
    return {};
  }
  if (token.size() < 2 || token.front() != '[' || token.back() != ']') {
    throw std::runtime_error("ParseJsonZZListField: expected array for key " + key);
  }
  std::vector<ZZ> out;
  std::size_t pos = 1;
  while (pos + 1 < token.size()) {
    while (pos + 1 < token.size() &&
           std::isspace(static_cast<unsigned char>(token[pos]))) {
      ++pos;
    }
    if (pos + 1 >= token.size() || token[pos] == ']') {
      break;
    }
    std::size_t end = pos;
    bool escape = false;
    if (token[pos] != '"') {
      throw std::runtime_error("ParseJsonZZListField: expected string element");
    }
    ++end;
    for (; end < token.size(); ++end) {
      const char ch = token[end];
      if (escape) {
        escape = false;
        continue;
      }
      if (ch == '\\') {
        escape = true;
        continue;
      }
      if (ch == '"') {
        break;
      }
    }
    if (end >= token.size()) {
      throw std::runtime_error("ParseJsonZZListField: unterminated string");
    }
    const std::string elem = ParseJsonStringToken(token.substr(pos, end - pos + 1));
    ZZ value;
    if (!basefold_bench_pcs_common::ParseZZString(elem, value)) {
      throw std::runtime_error("ParseJsonZZListField: invalid integer string");
    }
    out.push_back(value);
    pos = end + 1;
    while (pos < token.size() &&
           (std::isspace(static_cast<unsigned char>(token[pos])) || token[pos] == ',')) {
      ++pos;
    }
  }
  return out;
}

inline std::vector<std::vector<ZZ>> ParseJsonNestedZZListField(
    const std::string &json, const std::string &key) {
  const std::string token = ExtractJsonValueToken(json, key);
  if (token == "[]") {
    return {};
  }
  if (token.size() < 2 || token.front() != '[' || token.back() != ']') {
    throw std::runtime_error(
        "ParseJsonNestedZZListField: expected array for key " + key);
  }
  std::vector<std::vector<ZZ>> out;
  std::size_t pos = 1;
  while (pos + 1 < token.size()) {
    while (pos + 1 < token.size() &&
           (std::isspace(static_cast<unsigned char>(token[pos])) || token[pos] == ',')) {
      ++pos;
    }
    if (pos + 1 >= token.size() || token[pos] == ']') {
      break;
    }
    if (token[pos] != '[') {
      throw std::runtime_error("ParseJsonNestedZZListField: expected nested array");
    }
    std::size_t end = pos;
    int depth = 0;
    bool in_string = false;
    bool escape = false;
    for (; end < token.size(); ++end) {
      const char ch = token[end];
      if (in_string) {
        if (escape) {
          escape = false;
        } else if (ch == '\\') {
          escape = true;
        } else if (ch == '"') {
          in_string = false;
        }
        continue;
      }
      if (ch == '"') {
        in_string = true;
      } else if (ch == '[') {
        ++depth;
      } else if (ch == ']') {
        --depth;
        if (depth == 0) {
          break;
        }
      }
    }
    if (end >= token.size()) {
      throw std::runtime_error("ParseJsonNestedZZListField: unterminated nested array");
    }
    const std::string nested_json = std::string("{\"tmp\":") +
                                    token.substr(pos, end - pos + 1) + "}";
    out.push_back(ParseJsonZZListField(nested_json, "tmp"));
    pos = end + 1;
  }
  return out;
}

inline ArtifactMetadata ParseMetadataJson(const std::string &json) {
  ArtifactMetadata meta;
  meta.artifact_id = ParseJsonStringField(json, "artifact_id");
  meta.display_key = ParseJsonStringField(json, "display_key");
  meta.context_id = ParseJsonStringField(json, "context_id");
  meta.context_label = ParseJsonStringField(json, "context_label");
  meta.mode = ParseJsonStringField(json, "mode");
  meta.c = ParseJsonLongField(json, "c");
  meta.k0 = ParseJsonLongField(json, "k0");
  meta.d = ParseJsonLongField(json, "d");
  meta.poly_dim = ParseJsonLongField(json, "poly_dim");
  meta.lambda = ParseJsonStringField(json, "lambda");
  meta.gamma = ParseJsonStringField(json, "gamma");
  meta.queries = ParseJsonLongField(json, "queries");
  meta.seed = ParseJsonU64Field(json, "seed");
  meta.use_checked_prover_path =
      ParseJsonBoolField(json, "use_checked_prover_path");
  meta.use_extension_challenges =
      ParseJsonBoolField(json, "use_extension_challenges");
  if (!basefold_bench_pcs_common::ParseZZString(
          ParseJsonStringField(json, "scalar_modulus"), meta.scalar_modulus)) {
    throw std::runtime_error("ParseMetadataJson: bad scalar_modulus");
  }
  if (!basefold_bench_pcs_common::ParseZZString(
          ParseJsonStringField(json, "base_prime"), meta.base_prime)) {
    throw std::runtime_error("ParseMetadataJson: bad base_prime");
  }
  meta.F_coeffs = ParseJsonZZListField(json, "F_coeffs");
  meta.zeta_coeffs = ParseJsonZZListField(json, "zeta_coeffs");
  meta.zeta_source = ParseJsonStringField(json, "zeta_source");
  meta.challenge_extension_coeffs =
      ParseJsonNestedZZListField(json, "challenge_extension_coeffs");
  meta.hash_backend = ParseJsonStringField(json, "hash_backend");
  meta.proof_encoding = ParseJsonStringField(json, "proof_encoding");
  meta.proof_size_bytes = ParseJsonU64Field(json, "proof_size_bytes");
  return meta;
}

inline void WriteMetadataJson(const fs::path &path,
                              const ArtifactMetadata &meta) {
  WriteStringToFile(path, SerializeMetadataJson(meta));
}

inline ArtifactMetadata ReadMetadataJson(const fs::path &path) {
  return ParseMetadataJson(ReadFileToString(path));
}

inline void AppendManifestEntry(const fs::path &manifest_path,
                                const ArtifactManifestEntry &entry) {
  std::ofstream out(manifest_path, std::ios::binary | std::ios::app);
  if (!out) {
    throw std::runtime_error("AppendManifestEntry: failed to open " +
                             manifest_path.string());
  }
  out << SerializeManifestEntryJsonLine(entry) << '\n';
  if (!out) {
    throw std::runtime_error("AppendManifestEntry: failed to write " +
                             manifest_path.string());
  }
}

inline std::vector<ArtifactManifestEntry> LoadManifestEntries(
    const fs::path &manifest_path) {
  std::vector<ArtifactManifestEntry> out;
  std::ifstream in(manifest_path, std::ios::binary);
  if (!in) {
    throw std::runtime_error("LoadManifestEntries: failed to open " +
                             manifest_path.string());
  }
  std::string line;
  while (std::getline(in, line)) {
    line = TrimAscii(line);
    if (line.empty()) {
      continue;
    }
    ArtifactManifestEntry entry;
    entry.artifact_id = ParseJsonStringField(line, "artifact_id");
    entry.display_key = ParseJsonStringField(line, "display_key");
    entry.context_id = ParseJsonStringField(line, "context_id");
    entry.context_label = ParseJsonStringField(line, "context_label");
    entry.mode = ParseJsonStringField(line, "mode");
    entry.c = ParseJsonLongField(line, "c");
    entry.k0 = ParseJsonLongField(line, "k0");
    entry.d = ParseJsonLongField(line, "d");
    entry.poly_dim = ParseJsonLongField(line, "poly_dim");
    entry.lambda = ParseJsonStringField(line, "lambda");
    entry.gamma = ParseJsonStringField(line, "gamma");
    entry.queries = ParseJsonLongField(line, "queries");
    entry.seed = ParseJsonU64Field(line, "seed");
    entry.use_extension_challenges =
        ParseJsonBoolField(line, "use_extension_challenges");
    entry.object_relpath = ParseJsonStringField(line, "object_relpath");
    out.push_back(entry);
  }
  return out;
}

inline std::vector<basefold::Byte> DigestToBytes(const basefold::Digest &digest) {
  return std::vector<basefold::Byte>(digest.begin(), digest.end());
}

inline void AppendBytes(basefold::Bytes &out, const basefold::Byte *data,
                        std::size_t size) {
  if (size == 0) {
    return;
  }
  out.insert(out.end(), data, data + size);
}

inline void AppendU64LE(basefold::Bytes &out, std::uint64_t value) {
  for (int i = 0; i < 8; ++i) {
    out.push_back(static_cast<basefold::Byte>((value >> (8 * i)) & 0xff));
  }
}

inline std::uint64_t ReadU64LE(const basefold::Bytes &bytes, std::size_t &pos,
                               const char *what) {
  if (pos + 8 > bytes.size()) {
    throw std::runtime_error(std::string(what) + ": truncated u64");
  }
  std::uint64_t value = 0;
  for (int i = 0; i < 8; ++i) {
    value |= (static_cast<std::uint64_t>(bytes[pos + i]) << (8 * i));
  }
  pos += 8;
  return value;
}

inline basefold::FixedProofEncodingContext BuildPublicInputsEncodingContext() {
  basefold::BaseFoldPCSEvalProof dummy;
  return basefold::BuildFixedProofEncodingContext(dummy, {});
}

inline void SerializeFieldElementFixed(basefold::Bytes &out,
                                       const basefold::FieldElement &value,
                                       const basefold::FixedProofEncodingContext &ctx) {
  const NTL::ZZ_pX poly = NTL::rep(value);
  for (long i = 0; i < ctx.base_ext_degree; ++i) {
    const ZZ coeff = NTL::rep(NTL::coeff(poly, i));
    if (coeff < 0) {
      throw std::runtime_error("SerializeFieldElementFixed: negative coeff");
    }
    std::vector<unsigned char> tmp(static_cast<std::size_t>(ctx.coeff_bytes), 0);
    const long n = NTL::NumBytes(coeff);
    if (static_cast<std::uint64_t>(n) > ctx.coeff_bytes) {
      throw std::runtime_error("SerializeFieldElementFixed: coeff width overflow");
    }
    if (n > 0) {
      NTL::BytesFromZZ(tmp.data(), coeff, n);
    }
    AppendBytes(out, reinterpret_cast<const basefold::Byte *>(tmp.data()),
                tmp.size());
  }
}

inline basefold::FieldElement DeserializeFieldElementFixed(
    const basefold::Bytes &bytes, std::size_t &pos,
    const basefold::FixedProofEncodingContext &ctx, const char *what) {
  ZZ_pX poly;
  NTL::clear(poly);
  for (long i = 0; i < ctx.base_ext_degree; ++i) {
    if (pos + ctx.coeff_bytes > bytes.size()) {
      throw std::runtime_error(std::string(what) + ": truncated field element");
    }
    ZZ coeff = NTL::ZZFromBytes(bytes.data() + pos,
                                static_cast<long>(ctx.coeff_bytes));
    pos += ctx.coeff_bytes;
    if (coeff != 0) {
      SetCoeff(poly, i, conv<ZZ_p>(coeff));
    }
  }
  basefold::FieldElement out;
  conv(out, poly);
  return out;
}

inline basefold::Bytes SerializePublicInputs(const ArtifactPublicInputs &inputs) {
  const basefold::FixedProofEncodingContext ctx =
      BuildPublicInputsEncodingContext();
  basefold::Bytes out;
  out.reserve(ctx.hash_bytes + 8 +
              static_cast<std::size_t>(ctx.field_elem_bytes) *
                  static_cast<std::size_t>(inputs.z.size() + 1));
  AppendBytes(out, inputs.commitment_root.data(), inputs.commitment_root.size());
  AppendU64LE(out, static_cast<std::uint64_t>(inputs.z.size()));
  for (const basefold::FieldElement &zi : inputs.z) {
    SerializeFieldElementFixed(out, zi, ctx);
  }
  SerializeFieldElementFixed(out, inputs.claimed_y, ctx);
  return out;
}

inline ArtifactPublicInputs DeserializePublicInputs(const basefold::Bytes &bytes) {
  const basefold::FixedProofEncodingContext ctx =
      BuildPublicInputsEncodingContext();
  ArtifactPublicInputs out;
  std::size_t pos = 0;
  if (ctx.hash_bytes != out.commitment_root.size()) {
    throw std::runtime_error(
        "DeserializePublicInputs: digest width mismatch with MerkleRoot");
  }
  if (pos + ctx.hash_bytes > bytes.size()) {
    throw std::runtime_error("DeserializePublicInputs: truncated commitment_root");
  }
  for (std::size_t i = 0; i < out.commitment_root.size(); ++i) {
    out.commitment_root[i] = bytes[pos + i];
  }
  pos += ctx.hash_bytes;
  const std::uint64_t z_count =
      ReadU64LE(bytes, pos, "DeserializePublicInputs");
  out.z.reserve(static_cast<std::size_t>(z_count));
  for (std::uint64_t i = 0; i < z_count; ++i) {
    out.z.push_back(
        DeserializeFieldElementFixed(bytes, pos, ctx, "DeserializePublicInputs"));
  }
  out.claimed_y =
      DeserializeFieldElementFixed(bytes, pos, ctx, "DeserializePublicInputs");
  if (pos != bytes.size()) {
    throw std::runtime_error("DeserializePublicInputs: trailing bytes remain");
  }
  return out;
}

inline void WritePublicInputsBinary(const fs::path &path,
                                    const ArtifactPublicInputs &inputs) {
  const basefold::Bytes bytes = SerializePublicInputs(inputs);
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("WritePublicInputsBinary: failed to open " +
                             path.string());
  }
  out.write(reinterpret_cast<const char *>(bytes.data()),
            static_cast<std::streamsize>(bytes.size()));
  if (!out) {
    throw std::runtime_error("WritePublicInputsBinary: failed to write " +
                             path.string());
  }
}

inline ArtifactPublicInputs ReadPublicInputsBinary(const fs::path &path) {
  return DeserializePublicInputs(ReadBytesFromFile(path));
}

inline std::vector<ZZ> FieldElementToCoeffVector(const basefold::FieldElement &x) {
  std::vector<ZZ> coeffs;
  const ZZ_pX poly = NTL::rep(x);
  const long degree = NTL::deg(poly);
  if (degree < 0) {
    coeffs.push_back(ZZ(0));
    return coeffs;
  }
  coeffs.reserve(static_cast<std::size_t>(degree + 1));
  for (long i = 0; i <= degree; ++i) {
    coeffs.push_back(NTL::rep(NTL::coeff(poly, i)));
  }
  return coeffs;
}

inline std::vector<std::vector<ZZ>> ExtensionPolyToCoeffVectors(
    const ZZ_pEX &poly) {
  std::vector<std::vector<ZZ>> coeffs;
  const long degree = NTL::deg(poly);
  if (degree < 0) {
    return coeffs;
  }
  coeffs.reserve(static_cast<std::size_t>(degree + 1));
  for (long i = 0; i <= degree; ++i) {
    coeffs.push_back(FieldElementToCoeffVector(NTL::coeff(poly, i)));
  }
  return coeffs;
}

inline std::string CoeffVectorToken(const std::vector<ZZ> &coeffs) {
  std::ostringstream out;
  for (std::size_t i = 0; i < coeffs.size(); ++i) {
    if (i != 0) {
      out << '_';
    }
    out << coeffs[i];
  }
  return out.str();
}

inline std::string NestedCoeffVectorToken(
    const std::vector<std::vector<ZZ>> &coeffs) {
  std::ostringstream out;
  for (std::size_t i = 0; i < coeffs.size(); ++i) {
    if (i != 0) {
      out << "__";
    }
    out << CoeffVectorToken(coeffs[i]);
  }
  return out.str();
}

inline std::string BuildDefaultContextId(
    const std::string &mode, const basefold_bench_pcs_common::ContextSpec &spec,
    const std::vector<ZZ> &actual_zeta_coeffs, const std::string &zeta_source,
    const std::vector<std::vector<ZZ>> &actual_challenge_ext_coeffs) {
  std::ostringstream out;
  out << NormalizeModeLabel(mode) << "-mod" << spec.scalar_modulus;
  if (spec.base_prime > 1) {
    out << "-p" << spec.base_prime;
  }
  out << "-F" << CoeffVectorToken(spec.F_coeffs) << "-zeta"
      << CoeffVectorToken(actual_zeta_coeffs) << "-zsrc-" << zeta_source;
  if (!actual_challenge_ext_coeffs.empty()) {
    out << "-E" << NestedCoeffVectorToken(actual_challenge_ext_coeffs);
  }
  return out.str();
}

struct DumpArtifactRequest {
  fs::path artifact_root;
  std::string artifact_id;
  std::string context_id;
  std::string context_label;
  std::string mode;
  std::string lambda = "unspecified";
  std::string gamma = "unspecified";
  long c = 2;
  long k0 = 1;
  long d = 16;
  long queries = 4;
  std::uint64_t seed = 0;
  bool use_checked_prover_path = false;
  bool use_extension_challenges = false;
  bool auto_zeta_teich = false;
};

struct DumpArtifactResult {
  ArtifactMetadata metadata;
  fs::path manifest_path;
  fs::path object_dir;
  fs::path metadata_path;
  fs::path public_inputs_path;
  fs::path proof_path;
};

struct LoadedArtifactCase {
  ArtifactManifestEntry manifest_entry;
  ArtifactMetadata metadata;
  RestoredVerificationContext restored;
  ArtifactPublicInputs public_inputs;
  basefold::BaseFoldPCSEvalProof proof;
  long challenge_ext_degree = 0;
  double load_wall_ms = 0.0;
  double deserialize_wall_ms = 0.0;
};

struct ArtifactVerifyBenchResult {
  basefold_bench_pcs_common::Stats verifier;
  std::uint64_t anti_opt_checksum = 0;
  std::uint64_t proof_size_bytes = 0;
  double proof_size_kb = 0.0;
  basefold::Profile verifier_profile;
  bool has_profile = false;
};

inline DumpArtifactResult DumpEvalArtifact(
    const basefold_bench_pcs_common::ContextSpec &spec,
    const DumpArtifactRequest &request) {
  using namespace basefold_bench_pcs_common;

  if (request.artifact_root.empty()) {
    throw std::runtime_error("DumpEvalArtifact: artifact_root is required");
  }
  if (request.c <= 0) {
    throw std::runtime_error("DumpEvalArtifact: c must be positive");
  }
  if (request.k0 <= 0 || !IsPowerOfTwoLong(request.k0)) {
    throw std::runtime_error(
        "DumpEvalArtifact: k0 must be a positive power of two");
  }
  if (request.d < 0) {
    throw std::runtime_error("DumpEvalArtifact: d must be non-negative");
  }
  if (request.queries < 0) {
    throw std::runtime_error("DumpEvalArtifact: queries must be non-negative");
  }
  if (spec.scalar_modulus <= 1) {
    throw std::runtime_error("DumpEvalArtifact: scalar_modulus must be > 1");
  }

  const std::string mode = NormalizeModeLabel(request.mode);
  ZZ_pPush mod_push(spec.scalar_modulus);

  ValidateMonic(spec.F_coeffs, spec.scalar_modulus, "F");
  const ZZ_pX F = BuildZZpX(spec.F_coeffs);
  ZZ_pEPush e_push(F);

  ZZ_pE zeta;
  std::string zeta_source = "explicit";
  if (request.auto_zeta_teich) {
    ZZ p_base;
    long k_base = 0;
    DeduceBasePrimeAndExponent(spec, p_base, k_base);
    zeta = FindTeichmullerGenerator(p_base, k_base, NTL::deg(F), F);
    zeta_source = "auto_teich";
  } else {
    zeta = BuildZZpE(spec.zeta_coeffs);
  }
  const std::vector<ZZ> actual_zeta_coeffs = FieldElementToCoeffVector(zeta);

  const basefold::FoldableCodeParams params =
      (request.k0 == 1)
          ? BuildParams_k0_1(
                request.c, request.d,
                (spec.base_prime > 1) ? spec.base_prime : spec.scalar_modulus,
                zeta)
          : BuildParams_k0_pow2(
                request.c, request.k0, request.d,
                (spec.base_prime > 1) ? spec.base_prime : spec.scalar_modulus,
                zeta);

  const long poly_dim = ComputePolyDimOrThrow(request.k0, request.d);
  const vec_ZZ_pE f_coeffs = MakeDeterministicCoefficients(poly_dim, request.seed);
  const long point_dim = request.d + Log2ExactPowerOfTwoLong(request.k0);
  const std::vector<ZZ_pE> z =
      MakeDeterministicPoint(point_dim, request.seed ^ 0xdeadbeefULL);
  const ZZ_pE y = basefold::EvalMultilinearMonomialCoeffs(f_coeffs, z);

  basefold::BaseFoldPCSChallengeConfig challenge_cfg;
  ZZ_pEX challenge_modulus;
  long challenge_ext_degree = 0;
  const basefold::BaseFoldPCSChallengeConfig *challenge_cfg_ptr = nullptr;
  std::vector<std::vector<ZZ>> actual_challenge_ext_coeffs;
  if (request.use_extension_challenges) {
    if (spec.challenge_ext_coeffs.empty()) {
      NTL::clear(challenge_modulus);
      NTL::SetCoeff(challenge_modulus, 0, zeta);
      NTL::SetCoeff(challenge_modulus, 1, ZZ_pE(1));
      NTL::SetCoeff(challenge_modulus, 2, ZZ_pE(1));
    } else {
      challenge_modulus = BuildZZpEX(spec.challenge_ext_coeffs);
    }
    challenge_cfg.use_extension_challenges = true;
    challenge_cfg.challenge_extension_modulus = challenge_modulus;
    challenge_cfg_ptr = &challenge_cfg;
    challenge_ext_degree = NTL::deg(challenge_modulus);
    actual_challenge_ext_coeffs = ExtensionPolyToCoeffVectors(challenge_modulus);
  }

  const basefold::BaseFoldPCSCommitArtifacts commit_artifacts =
      basefold::BaseFoldPCSBuildCommitArtifactsUnchecked(f_coeffs, params);
  const basefold::BaseFoldPCSEvalProof proof =
      (challenge_cfg_ptr != nullptr)
          ? (request.use_checked_prover_path
                 ? basefold::BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracle(
                       f_coeffs, z, y, request.queries, params,
                       commit_artifacts, *challenge_cfg_ptr)
                 : basefold::BaseFoldPCSProveEvalWithChallengeConfigFromCommittedTopOracleUnchecked(
                       f_coeffs, z, y, request.queries, params,
                       commit_artifacts, *challenge_cfg_ptr))
          : (request.use_checked_prover_path
                 ? basefold::BaseFoldPCSProveEvalFromCommittedTopOracle(
                       f_coeffs, z, y, request.queries, params,
                       commit_artifacts)
                 : basefold::BaseFoldPCSProveEvalFromCommittedTopOracleUnchecked(
                       f_coeffs, z, y, request.queries, params,
                       commit_artifacts));

  ArtifactPublicInputs public_inputs;
  public_inputs.commitment_root = commit_artifacts.root_d;
  public_inputs.z = z;
  public_inputs.claimed_y = y;

  basefold::FixedProofEncodingOptions encoding_options;
  encoding_options.include_version_byte = true;
  if (request.use_extension_challenges || proof.extension.has_extension_payload) {
    encoding_options.challenge_ext_degree = challenge_ext_degree;
  }
  const basefold::Bytes proof_bytes =
      basefold::SerializeBaseFoldPCSEvalProofFixedBytes(proof, encoding_options);
  const std::uint64_t expected_proof_size_bytes = ComputeProofSizeBytes(
      proof, request.use_extension_challenges, challenge_ext_degree);
  if (proof_bytes.size() != expected_proof_size_bytes) {
    throw std::runtime_error(
        "DumpEvalArtifact: serializer byte size does not match proof size accounting");
  }

  ArtifactMetadata metadata;
  metadata.context_label =
      request.context_label.empty() ? spec.label : request.context_label;
  metadata.context_id =
      request.context_id.empty()
          ? BuildDefaultContextId(mode, spec, actual_zeta_coeffs, zeta_source,
                                  actual_challenge_ext_coeffs)
          : request.context_id;
  metadata.mode = mode;
  metadata.c = request.c;
  metadata.k0 = request.k0;
  metadata.d = request.d;
  metadata.poly_dim = poly_dim;
  metadata.lambda = request.lambda;
  metadata.gamma = request.gamma;
  metadata.queries = request.queries;
  metadata.seed = request.seed;
  metadata.use_checked_prover_path = request.use_checked_prover_path;
  metadata.use_extension_challenges = request.use_extension_challenges;
  metadata.scalar_modulus = spec.scalar_modulus;
  metadata.base_prime = spec.base_prime;
  metadata.F_coeffs = spec.F_coeffs;
  metadata.zeta_coeffs = actual_zeta_coeffs;
  metadata.zeta_source = zeta_source;
  metadata.challenge_extension_coeffs = actual_challenge_ext_coeffs;
  metadata.hash_backend = basefold::SelectedHashBackendName();
  metadata.proof_encoding = "basefold_fixed_v1";
  metadata.proof_size_bytes = expected_proof_size_bytes;
  metadata.artifact_id = request.artifact_id.empty()
                             ? ComputeCanonicalArtifactId(metadata)
                             : request.artifact_id;
  ValidateArtifactIdOrThrow(metadata.artifact_id);
  metadata.display_key = BuildArtifactDisplayKey(metadata);

  const fs::path manifest_path = ArtifactManifestPath(request.artifact_root);
  const fs::path object_dir =
      ArtifactObjectDir(request.artifact_root, metadata.artifact_id);
  const fs::path metadata_path =
      ArtifactMetadataPath(request.artifact_root, metadata.artifact_id);
  const fs::path public_inputs_path =
      ArtifactPublicInputsPath(request.artifact_root, metadata.artifact_id);
  const fs::path proof_path =
      ArtifactProofPath(request.artifact_root, metadata.artifact_id);

  if (fs::exists(object_dir)) {
    throw std::runtime_error("DumpEvalArtifact: artifact already exists at " +
                             object_dir.string());
  }
  if (fs::exists(manifest_path)) {
    const std::vector<ArtifactManifestEntry> entries =
        LoadManifestEntries(manifest_path);
    for (const ArtifactManifestEntry &entry : entries) {
      if (entry.artifact_id == metadata.artifact_id) {
        throw std::runtime_error(
            "DumpEvalArtifact: artifact_id already exists in manifest");
      }
    }
  }

  fs::create_directories(request.artifact_root / "objects");
  fs::create_directories(object_dir);
  WriteMetadataJson(metadata_path, metadata);
  WritePublicInputsBinary(public_inputs_path, public_inputs);
  WriteBytesToFile(proof_path, proof_bytes);
  AppendManifestEntry(manifest_path, ManifestEntryFromMetadata(metadata));

  DumpArtifactResult result;
  result.metadata = metadata;
  result.manifest_path = manifest_path;
  result.object_dir = object_dir;
  result.metadata_path = metadata_path;
  result.public_inputs_path = public_inputs_path;
  result.proof_path = proof_path;
  return result;
}

inline ArtifactManifestEntry FindManifestEntryOrThrow(
    const std::vector<ArtifactManifestEntry> &entries,
    const std::string &artifact_id) {
  const ArtifactManifestEntry *match = nullptr;
  for (const ArtifactManifestEntry &entry : entries) {
    if (entry.artifact_id != artifact_id) {
      continue;
    }
    if (match != nullptr) {
      throw std::runtime_error(
          "FindManifestEntryOrThrow: duplicate artifact_id in manifest");
    }
    match = &entry;
  }
  if (match == nullptr) {
    throw std::runtime_error(
        "FindManifestEntryOrThrow: artifact_id not found in manifest");
  }
  return *match;
}

inline void ValidateManifestEntryMatchesMetadata(
    const ArtifactManifestEntry &entry, const ArtifactMetadata &meta) {
  if (entry.artifact_id != meta.artifact_id ||
      entry.display_key != meta.display_key || entry.context_id != meta.context_id ||
      entry.context_label != meta.context_label || entry.mode != meta.mode ||
      entry.c != meta.c || entry.k0 != meta.k0 || entry.d != meta.d ||
      entry.poly_dim != meta.poly_dim || entry.lambda != meta.lambda ||
      entry.gamma != meta.gamma || entry.queries != meta.queries ||
      entry.seed != meta.seed ||
      entry.use_extension_challenges != meta.use_extension_challenges) {
    throw std::runtime_error(
        "ValidateManifestEntryMatchesMetadata: manifest/meta mismatch");
  }
  if (entry.object_relpath != ManifestObjectRelPath(meta.artifact_id)) {
    throw std::runtime_error(
        "ValidateManifestEntryMatchesMetadata: unexpected object_relpath");
  }
}

inline bool VerifyLoadedArtifactCase(const LoadedArtifactCase &artifact) {
  return artifact.metadata.use_extension_challenges
             ? basefold::BaseFoldPCSVerifyEvalWithChallengeConfig(
                   artifact.public_inputs.commitment_root, artifact.public_inputs.z,
                   artifact.public_inputs.claimed_y, artifact.metadata.queries,
                   artifact.proof, artifact.restored.params,
                   artifact.restored.challenge_cfg)
             : basefold::BaseFoldPCSVerifyEval(
                   artifact.public_inputs.commitment_root, artifact.public_inputs.z,
                   artifact.public_inputs.claimed_y, artifact.metadata.queries,
                   artifact.proof, artifact.restored.params);
}

inline LoadedArtifactCase LoadArtifactCaseForVerify(const fs::path &artifact_root,
                                                    const std::string &artifact_id) {
  using namespace basefold_bench_pcs_common;

  ValidateArtifactIdOrThrow(artifact_id);
  if (artifact_root.empty()) {
    throw std::runtime_error(
        "LoadArtifactCaseForVerify: artifact_root is required");
  }

  const auto load_t0 = std::chrono::steady_clock::now();
  const std::vector<ArtifactManifestEntry> entries =
      LoadManifestEntries(ArtifactManifestPath(artifact_root));
  const ArtifactManifestEntry entry = FindManifestEntryOrThrow(entries, artifact_id);
  const fs::path object_dir = ArtifactObjectDir(artifact_root, artifact_id);
  const std::string metadata_json = ReadFileToString(object_dir / "meta.json");
  const basefold::Bytes public_inputs_bytes =
      ReadBytesFromFile(object_dir / "public_inputs.bin");
  const basefold::Bytes proof_bytes = ReadBytesFromFile(object_dir / "proof.bin");
  const auto load_t1 = std::chrono::steady_clock::now();

  const auto deserialize_t0 = std::chrono::steady_clock::now();
  LoadedArtifactCase out;
  out.manifest_entry = entry;
  out.metadata = ParseMetadataJson(metadata_json);
  if (out.metadata.artifact_id != artifact_id) {
    throw std::runtime_error(
        "LoadArtifactCaseForVerify: meta.json artifact_id mismatch");
  }
  ValidateManifestEntryMatchesMetadata(out.manifest_entry, out.metadata);
  if (out.metadata.hash_backend != basefold::SelectedHashBackendName()) {
    throw std::runtime_error(
        "LoadArtifactCaseForVerify: artifact hash_backend does not match current build");
  }
  if (out.metadata.proof_encoding != "basefold_fixed_v1") {
    throw std::runtime_error(
        "LoadArtifactCaseForVerify: unsupported proof_encoding");
  }
  out.restored = RestoreVerificationContext(out.metadata);
  out.public_inputs = DeserializePublicInputs(public_inputs_bytes);
  basefold::FixedProofEncodingOptions options;
  options.include_version_byte = true;
  if (out.metadata.use_extension_challenges) {
    out.challenge_ext_degree =
        NTL::deg(out.restored.challenge_cfg.challenge_extension_modulus);
    options.challenge_ext_degree = out.challenge_ext_degree;
  }
  out.proof =
      basefold::DeserializeBaseFoldPCSEvalProofFixedBytes(proof_bytes, options);
  if (proof_bytes.size() != out.metadata.proof_size_bytes) {
    throw std::runtime_error(
        "LoadArtifactCaseForVerify: proof.bin length does not match meta.json");
  }
  const auto deserialize_t1 = std::chrono::steady_clock::now();

  out.load_wall_ms = MsSince(load_t0, load_t1);
  out.deserialize_wall_ms = MsSince(deserialize_t0, deserialize_t1);
  return out;
}

inline ArtifactVerifyBenchResult RunArtifactVerifyBenchmark(
    const LoadedArtifactCase &artifact, bool enable_profile, int warmup,
    int reps) {
  using namespace basefold_bench_pcs_common;

  if (warmup < 0) {
    throw std::runtime_error(
        "RunArtifactVerifyBenchmark: warmup must be >= 0");
  }
  if (reps <= 0) {
    throw std::runtime_error("RunArtifactVerifyBenchmark: reps must be > 0");
  }

  std::vector<double> verifier_ms;
  verifier_ms.reserve(static_cast<std::size_t>(reps));

  basefold::Profile verifier_prof;
  basefold::ResetProfile(verifier_prof);

  std::uint64_t anti_opt_checksum = 0;
  for (int iter = -warmup; iter < reps; ++iter) {
    const auto t0 = std::chrono::steady_clock::now();
    const bool ok = [&] {
      if (enable_profile && iter >= 0) {
        basefold::ProfileGuard guard(&verifier_prof);
        return VerifyLoadedArtifactCase(artifact);
      }
      return VerifyLoadedArtifactCase(artifact);
    }();
    const auto t1 = std::chrono::steady_clock::now();

    if (!ok) {
      throw std::runtime_error(
          "RunArtifactVerifyBenchmark: verification failed");
    }

    if (!artifact.public_inputs.commitment_root.empty()) {
      anti_opt_checksum ^=
          static_cast<std::uint64_t>(artifact.public_inputs.commitment_root[0]);
    }
    anti_opt_checksum ^= static_cast<std::uint64_t>(ok);
    anti_opt_checksum ^= artifact.metadata.proof_size_bytes;

    if (iter >= 0) {
      verifier_ms.push_back(MsSince(t0, t1));
    }
  }

  ArtifactVerifyBenchResult out;
  out.verifier = ComputeStats(verifier_ms);
  out.anti_opt_checksum = anti_opt_checksum;
  out.proof_size_bytes = artifact.metadata.proof_size_bytes;
  out.proof_size_kb =
      static_cast<double>(artifact.metadata.proof_size_bytes) / 1024.0;
  out.verifier_profile = verifier_prof;
  out.has_profile = enable_profile;
  return out;
}

inline RestoredVerificationContext RestoreVerificationContext(
    const ArtifactMetadata &meta) {
  RestoredVerificationContext out;
  const std::string mode = NormalizeModeLabel(meta.mode);
  ZZ_p::init(meta.scalar_modulus);
  const ZZ_pX F = basefold_bench_pcs_common::BuildZZpX(meta.F_coeffs);
  ZZ_pE::init(F);

  const ZZ_pE zeta = basefold_bench_pcs_common::BuildZZpE(meta.zeta_coeffs);
  if (mode == "ring") {
    ZZ p_base;
    long k_base = 0;
    basefold_bench_pcs_common::DeduceBasePrimeAndExponent(
        {meta.context_label, meta.scalar_modulus, meta.base_prime,
         meta.F_coeffs, meta.zeta_coeffs, meta.challenge_extension_coeffs},
        p_base, k_base);
    (void)k_base;
    out.params = (meta.k0 == 1)
                     ? basefold_bench_pcs_common::BuildParams_k0_1(
                           meta.c, meta.d, p_base, zeta)
                     : basefold_bench_pcs_common::BuildParams_k0_pow2(
                           meta.c, meta.k0, meta.d, p_base, zeta);
  } else {
    out.params = (meta.k0 == 1)
                     ? basefold_bench_pcs_common::BuildParams_k0_1(
                           meta.c, meta.d, meta.scalar_modulus, zeta)
                     : basefold_bench_pcs_common::BuildParams_k0_pow2(
                           meta.c, meta.k0, meta.d, meta.scalar_modulus, zeta);
  }

  if (meta.use_extension_challenges) {
    out.challenge_cfg.use_extension_challenges = true;
    out.challenge_cfg.challenge_extension_modulus =
        basefold_bench_pcs_common::BuildZZpEX(meta.challenge_extension_coeffs);
  }
  return out;
}

inline std::string DescribeArtifactSummary(const ArtifactMetadata &meta) {
  std::ostringstream out;
  out << "artifact_id=" << meta.artifact_id << '\n'
      << "display_key=" << meta.display_key << '\n'
      << "mode=" << meta.mode << " c=" << meta.c << " k0=" << meta.k0
      << " d=" << meta.d << " poly_dim=" << meta.poly_dim
      << " queries=" << meta.queries << '\n'
      << "seed=" << meta.seed
      << " ext_challenges=" << (meta.use_extension_challenges ? "on" : "off")
      << " proof_size_bytes=" << meta.proof_size_bytes;
  return out.str();
}

}  // namespace basefold_bench_pcs_artifact

#endif  // BASEFOLD_BENCH_PCS_ARTIFACT_COMMON_HPP_
