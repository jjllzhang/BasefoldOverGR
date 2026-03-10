#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "tests/test_common.hpp"

int g_test_failure_count = 0;

namespace {

namespace fs = std::filesystem;

struct ForbiddenPattern {
  std::string needle;
  std::string message;
};

fs::path FindRepoRoot(fs::path start) {
  start = fs::absolute(start);
  while (!start.empty()) {
    if (fs::exists(start / "CMakeLists.txt") && fs::exists(start / "include") &&
        fs::exists(start / "src")) {
      return start;
    }
    const fs::path parent = start.parent_path();
    if (parent == start) {
      break;
    }
    start = parent;
  }
  return fs::current_path();
}

bool HasAllowedExtension(const fs::path &path) {
  static const std::vector<std::string> kAllowedExtensions = {
      ".c",   ".cc",  ".cmake", ".cpp", ".h",   ".hpp",
      ".md",  ".py",  ".sh",    ".txt", ".yml", ".yaml",
  };
  const std::string ext = path.extension().string();
  for (const std::string &allowed : kAllowedExtensions) {
    if (ext == allowed) {
      return true;
    }
  }
  return false;
}

std::string ReadFileOrEmpty(const fs::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) {
    return std::string();
  }
  return std::string((std::istreambuf_iterator<char>(in)),
                     std::istreambuf_iterator<char>());
}

void CheckFileForForbiddenPatterns(
    const fs::path &path, const std::vector<ForbiddenPattern> &patterns) {
  const std::string content = ReadFileOrEmpty(path);
  CHECK_MSG(!content.empty(),
            std::string("Failed to read text file: ") + path.string());

  for (const ForbiddenPattern &pattern : patterns) {
    if (content.find(pattern.needle) != std::string::npos) {
      CHECK_MSG(false,
                path.string() + " contains forbidden legacy path pattern '" +
                    pattern.needle + "' (" + pattern.message + ")");
    }
  }
}

void VisitTreeAndCheck(const fs::path &root,
                       const std::vector<ForbiddenPattern> &patterns) {
  if (!fs::exists(root)) {
    return;
  }
  for (const fs::directory_entry &entry :
       fs::recursive_directory_iterator(root)) {
    if (!entry.is_regular_file()) {
      continue;
    }
    const fs::path path = entry.path();
    if (path.filename() == "test_path_hygiene.cpp") {
      continue;
    }
    if (!HasAllowedExtension(path)) {
      continue;
    }
    CheckFileForForbiddenPatterns(path, patterns);
  }
}

void TestNoLegacyBaseFoldPathNaming() {
  testutil::PrintInfo("Path hygiene: reject legacy include/src path naming");

  const std::vector<ForbiddenPattern> patterns = {
      {std::string("#include \"") + "Base" + "Fold/",
       "legacy header include path"},
      {std::string("#include <") + "Base" + "Fold/",
       "legacy header include path"},
      {std::string("include/") + "Base" + "Fold/",
       "legacy include directory reference"},
      {std::string("src/") + "Base" + "Fold/",
       "legacy source directory reference"},
  };

  const fs::path repo_root = FindRepoRoot(fs::current_path());
  const std::vector<fs::path> roots = {
      repo_root / "include",
      repo_root / "src",
      repo_root / "tests",
      repo_root / "bench",
      repo_root / "scripts",
  };
  for (const fs::path &root : roots) {
    VisitTreeAndCheck(root, patterns);
  }

  const std::vector<fs::path> top_level_files = {
      repo_root / "CMakeLists.txt",
      repo_root / "README.md",
      repo_root / "README_zh.md",
  };
  for (const fs::path &path : top_level_files) {
    CheckFileForForbiddenPatterns(path, patterns);
  }
}

}  // namespace

int main() {
  RUN_TEST(TestNoLegacyBaseFoldPathNaming);

  if (g_test_failure_count != 0) {
    std::cerr << "\n" << g_test_failure_count << " test assertion(s) failed.\n";
    return 1;
  }

  std::cout << "\nAll tests passed.\n";
  return 0;
}
