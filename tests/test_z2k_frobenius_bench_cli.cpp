#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "tests/test_common.hpp"

int g_test_failure_count = 0;

namespace {

namespace fs = std::filesystem;

fs::path g_executable_dir;
constexpr const char *kRingF64 =
    "1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,"
    "1,1,0,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,0,1,1,0,0,0,0,0,1,0,1";

struct CommandResult {
  int exit_code = -1;
  bool exited = false;
  bool signaled = false;
  std::string output;
};

fs::path DetectExecutableDir(const char *argv0) {
  if (argv0 == nullptr || std::string(argv0).empty()) {
    return fs::current_path();
  }
  return fs::absolute(fs::path(argv0)).parent_path();
}

CommandResult RunCommandCapture(const fs::path &exe,
                                const std::vector<std::string> &args) {
  CommandResult result;
#if defined(__unix__) || defined(__APPLE__)
  int pipe_fds[2];
  CHECK_MSG(pipe(pipe_fds) == 0,
            std::string("RunCommandCapture: pipe() failed for ") + exe.string());
  if (g_test_failure_count != 0) {
    return result;
  }

  std::cout.flush();
  std::cerr.flush();
  const pid_t pid = fork();
  CHECK_MSG(pid >= 0,
            std::string("RunCommandCapture: fork() failed for ") + exe.string());
  if (g_test_failure_count != 0) {
    close(pipe_fds[0]);
    close(pipe_fds[1]);
    return result;
  }

  if (pid == 0) {
    close(pipe_fds[0]);
    dup2(pipe_fds[1], STDOUT_FILENO);
    dup2(pipe_fds[1], STDERR_FILENO);
    close(pipe_fds[1]);

    std::vector<std::string> argv_storage;
    argv_storage.reserve(args.size() + 1);
    argv_storage.push_back(exe.string());
    for (const std::string &arg : args) {
      argv_storage.push_back(arg);
    }

    std::vector<char *> argv_ptrs;
    argv_ptrs.reserve(argv_storage.size() + 1);
    for (std::string &arg : argv_storage) {
      argv_ptrs.push_back(arg.data());
    }
    argv_ptrs.push_back(nullptr);

    execv(exe.c_str(), argv_ptrs.data());
    std::perror("execv");
    _exit(127);
  }

  close(pipe_fds[1]);
  char buffer[256];
  while (true) {
    const ssize_t bytes_read = read(pipe_fds[0], buffer, sizeof(buffer));
    if (bytes_read <= 0) {
      break;
    }
    result.output.append(buffer, static_cast<std::size_t>(bytes_read));
  }
  close(pipe_fds[0]);

  int status = 0;
  CHECK_MSG(waitpid(pid, &status, 0) == pid,
            std::string("RunCommandCapture: waitpid() failed for ") + exe.string());
  if (g_test_failure_count != 0) {
    return result;
  }

  result.exited = WIFEXITED(status);
  result.signaled = WIFSIGNALED(status);
  if (result.exited) {
    result.exit_code = WEXITSTATUS(status);
  }
  return result;
#else
  (void)exe;
  (void)args;
  testutil::PrintInfo(
      "RunCommandCapture: skipped subprocess execution on this platform");
  return result;
#endif
}

void ExpectCommandSuccessContains(const fs::path &exe,
                                  const std::vector<std::string> &args,
                                  const std::vector<std::string> &needles,
                                  const std::string &label) {
#if defined(__unix__) || defined(__APPLE__)
  const CommandResult result = RunCommandCapture(exe, args);
  CHECK_MSG(result.exited && result.exit_code == 0,
            label + ": command failed with output:\n" + result.output);
  for (const std::string &needle : needles) {
    CHECK_MSG(result.output.find(needle) != std::string::npos,
              label + ": missing output needle '" + needle + "' in:\n" +
                  result.output);
  }
#else
  (void)exe;
  (void)args;
  (void)needles;
  testutil::PrintInfo(label + ": skipped subprocess assertion on this platform");
#endif
}

std::uint64_t ExtractByteCountFromOutput(const std::string &output,
                                         const std::string &label) {
  const std::size_t label_pos = output.find(label);
  CHECK_MSG(label_pos != std::string::npos,
            "ExtractByteCountFromOutput: missing label " + label);
  if (g_test_failure_count != 0) {
    return 0;
  }
  const std::size_t lparen = output.find('(', label_pos);
  const std::size_t bytes_suffix = output.find(" B)", lparen);
  CHECK_MSG(lparen != std::string::npos && bytes_suffix != std::string::npos,
            "ExtractByteCountFromOutput: malformed byte field for " + label);
  if (g_test_failure_count != 0) {
    return 0;
  }

  const std::string digits =
      output.substr(lparen + 1, bytes_suffix - (lparen + 1));
  try {
    return static_cast<std::uint64_t>(std::stoull(digits));
  } catch (...) {
    CHECK_MSG(false, "ExtractByteCountFromOutput: invalid byte count for " + label);
    return 0;
  }
}

void TestFrobeniusBenchCommit_HelpTextIsStable() {
  testutil::PrintInfo("Frobenius bench CLI: commit --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_commit";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_frobenius_commit", "--auto-zeta teich", "Headline commit time"},
      "TestFrobeniusBenchCommit_HelpTextIsStable");
}

void TestFrobeniusBenchEval_HelpTextIsStable() {
  testutil::PrintInfo("Frobenius bench CLI: eval --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_eval";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_frobenius_eval", "--queries <int>",
       "serializer-backed bytes"},
      "TestFrobeniusBenchEval_HelpTextIsStable");
}

void TestFrobeniusBenchOuterCommit_HelpTextIsStable() {
  testutil::PrintInfo(
      "Frobenius bench CLI: outer commit --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_outer_commit";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_frobenius_outer_commit", "--auto-zeta teich",
       "No backend PCS commit"},
      "TestFrobeniusBenchOuterCommit_HelpTextIsStable");
}

void TestFrobeniusBenchOuterProve_HelpTextIsStable() {
  testutil::PrintInfo(
      "Frobenius bench CLI: outer prove --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_outer_prove";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_frobenius_outer_prove", "--queries <int>",
       "outer proof generation"},
      "TestFrobeniusBenchOuterProve_HelpTextIsStable");
}

void TestFrobeniusBenchOuterVerify_HelpTextIsStable() {
  testutil::PrintInfo(
      "Frobenius bench CLI: outer verify --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_outer_verify";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_frobenius_outer_verify", "--queries <int>",
       "outer verification"},
      "TestFrobeniusBenchOuterVerify_HelpTextIsStable");
}

void TestFrobeniusBenchCommit_SmokeRunPrintsStableFields() {
  testutil::PrintInfo("Frobenius bench CLI: commit smoke run prints stable result fields");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_commit";
  ExpectCommandSuccessContains(
      exe, {"--warmup", "0", "--reps", "1", "--seed", "7"},
      {"[frobenius commit]", "hash backend",
       "basis mode bench-fixed preferred-normal",
       "basis search distinct candidates", "basis score", "packing mean",
       "commit  mean", "anti-opt checksum"},
      "TestFrobeniusBenchCommit_SmokeRunPrintsStableFields");
}

void TestFrobeniusBenchCommit_HardcodedPresetForGR2P16_64() {
  testutil::PrintInfo(
      "Frobenius bench CLI: commit GR(2^16,64) hits the hardcoded commit-like preset");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_commit";
  ExpectCommandSuccessContains(
      exe,
      {"--ell", "7", "--kappa", "6", "--warmup", "0", "--reps", "1",
       "--seed", "23", "--ring-mod", "65536", "--ring-p", "2", "--ring-F",
       kRingF64, "--ring-zeta", "0,1"},
      {"[frobenius commit]", "basis mode bench-hardcoded preferred-normal",
       "basis preset commit-like (a0=1, a1=1, exponent=14)",
       "anti-opt checksum"},
      "TestFrobeniusBenchCommit_HardcodedPresetForGR2P16_64");
}

void TestFrobeniusBenchOuterCommit_SmokeRunPrintsStableFields() {
  testutil::PrintInfo(
      "Frobenius bench CLI: outer commit smoke run prints stable result fields");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_outer_commit";
  ExpectCommandSuccessContains(
      exe, {"--warmup", "0", "--reps", "1", "--seed", "13"},
      {"[frobenius outer commit]", "hash backend",
       "basis mode bench-fixed preferred-normal",
       "basis search distinct candidates", "basis score", "packing mean",
       "commit  mean", "backend input size", "anti-opt checksum"},
      "TestFrobeniusBenchOuterCommit_SmokeRunPrintsStableFields");
}

void TestFrobeniusBenchEval_SmokeRunPrintsStableFieldsAndSizes() {
  testutil::PrintInfo("Frobenius bench CLI: eval smoke run prints stable fields and sane proof sizes");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_eval";
#if defined(__unix__) || defined(__APPLE__)
  const CommandResult result =
      RunCommandCapture(exe, {"--warmup", "0", "--reps", "1", "--queries", "2",
                              "--seed", "11"});
  CHECK_MSG(result.exited && result.exit_code == 0,
            "TestFrobeniusBenchEval_SmokeRunPrintsStableFieldsAndSizes: command failed with output:\n" +
                result.output);
  if (g_test_failure_count != 0) {
    return;
  }

  const std::vector<std::string> needles = {
      "[frobenius eval]", "hash backend",
      "basis mode bench-fixed preferred-normal",
      "basis search distinct candidates", "basis score", "prove-phase mean",
      "outer prover mean", "backend prover mean", "verifier mean",
      "outer verifier mean", "backend verifier mean", "outer proof size",
      "proof size", "anti-opt checksum"};
  for (const std::string &needle : needles) {
    CHECK_MSG(result.output.find(needle) != std::string::npos,
              "TestFrobeniusBenchEval_SmokeRunPrintsStableFieldsAndSizes: missing output needle '" +
                  needle + "' in:\n" + result.output);
  }
  if (g_test_failure_count != 0) {
    return;
  }

  const std::uint64_t outer_bytes =
      ExtractByteCountFromOutput(result.output, "outer proof size");
  const std::uint64_t proof_bytes =
      ExtractByteCountFromOutput(result.output, "proof size");
  CHECK_GT(outer_bytes, 0U);
  CHECK_GE(proof_bytes, outer_bytes);
#else
  testutil::PrintInfo(
      "TestFrobeniusBenchEval_SmokeRunPrintsStableFieldsAndSizes: skipped subprocess assertion on this platform");
#endif
}

void TestFrobeniusBenchEval_HardcodedPresetForGR2P32_64() {
  testutil::PrintInfo(
      "Frobenius bench CLI: eval GR(2^32,64) hits the hardcoded eval-like preset");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_eval";
#if defined(__unix__) || defined(__APPLE__)
  const CommandResult result =
      RunCommandCapture(exe, {"--ell", "7", "--kappa", "6", "--warmup", "0",
                              "--reps", "1", "--queries", "2", "--seed",
                              "29", "--ring-mod", "4294967296", "--ring-p",
                              "2", "--ring-F", kRingF64, "--ring-zeta",
                              "0,1"});
  CHECK_MSG(result.exited && result.exit_code == 0,
            "TestFrobeniusBenchEval_HardcodedPresetForGR2P32_64: command failed with output:\n" +
                result.output);
  if (g_test_failure_count != 0) {
    return;
  }

  const std::vector<std::string> needles = {
      "[frobenius eval]", "basis mode bench-hardcoded preferred-normal",
      "basis preset eval-like (a0=1, a1=1, exponent=8)", "proof size",
      "anti-opt checksum"};
  for (const std::string &needle : needles) {
    CHECK_MSG(result.output.find(needle) != std::string::npos,
              "TestFrobeniusBenchEval_HardcodedPresetForGR2P32_64: missing output needle '" +
                  needle + "' in:\n" + result.output);
  }
#else
  testutil::PrintInfo(
      "TestFrobeniusBenchEval_HardcodedPresetForGR2P32_64: skipped subprocess assertion on this platform");
#endif
}

void TestFrobeniusBenchOuterProve_SmokeRunPrintsStableFieldsAndSizes() {
  testutil::PrintInfo(
      "Frobenius bench CLI: outer prove smoke run prints stable fields and sane proof sizes");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_outer_prove";
#if defined(__unix__) || defined(__APPLE__)
  const CommandResult result =
      RunCommandCapture(exe, {"--warmup", "0", "--reps", "1", "--queries", "2",
                              "--seed", "17"});
  CHECK_MSG(result.exited && result.exit_code == 0,
            "TestFrobeniusBenchOuterProve_SmokeRunPrintsStableFieldsAndSizes: command failed with output:\n" +
                result.output);
  if (g_test_failure_count != 0) {
    return;
  }

  const std::vector<std::string> needles = {
      "[frobenius outer prove]", "hash backend",
      "basis mode bench-fixed preferred-normal",
      "basis search distinct candidates", "basis score", "prove-phase mean",
      "proof size", "anti-opt checksum"};
  for (const std::string &needle : needles) {
    CHECK_MSG(result.output.find(needle) != std::string::npos,
              "TestFrobeniusBenchOuterProve_SmokeRunPrintsStableFieldsAndSizes: missing output needle '" +
                  needle + "' in:\n" + result.output);
  }
  if (g_test_failure_count != 0) {
    return;
  }

  const std::uint64_t proof_bytes =
      ExtractByteCountFromOutput(result.output, "proof size");
  CHECK_GT(proof_bytes, 0U);
#else
  testutil::PrintInfo(
      "TestFrobeniusBenchOuterProve_SmokeRunPrintsStableFieldsAndSizes: skipped subprocess assertion on this platform");
#endif
}

void TestFrobeniusBenchOuterVerify_SmokeRunPrintsStableFieldsAndSizes() {
  testutil::PrintInfo(
      "Frobenius bench CLI: outer verify smoke run prints stable fields and sane proof sizes");

  const fs::path exe = g_executable_dir / "bench_z2k_frobenius_outer_verify";
#if defined(__unix__) || defined(__APPLE__)
  const CommandResult result =
      RunCommandCapture(exe, {"--warmup", "0", "--reps", "1", "--queries", "2",
                              "--seed", "19"});
  CHECK_MSG(result.exited && result.exit_code == 0,
            "TestFrobeniusBenchOuterVerify_SmokeRunPrintsStableFieldsAndSizes: command failed with output:\n" +
                result.output);
  if (g_test_failure_count != 0) {
    return;
  }

  const std::vector<std::string> needles = {
      "[frobenius outer verify]", "hash backend",
      "basis mode bench-fixed preferred-normal",
      "basis search distinct candidates", "basis score", "verifier mean",
      "input proof size", "anti-opt checksum"};
  for (const std::string &needle : needles) {
    CHECK_MSG(result.output.find(needle) != std::string::npos,
              "TestFrobeniusBenchOuterVerify_SmokeRunPrintsStableFieldsAndSizes: missing output needle '" +
                  needle + "' in:\n" + result.output);
  }
  if (g_test_failure_count != 0) {
    return;
  }

  const std::uint64_t proof_bytes =
      ExtractByteCountFromOutput(result.output, "input proof size");
  CHECK_GT(proof_bytes, 0U);
#else
  testutil::PrintInfo(
      "TestFrobeniusBenchOuterVerify_SmokeRunPrintsStableFieldsAndSizes: skipped subprocess assertion on this platform");
#endif
}

}  // namespace

int main(int argc, char **argv) {
  g_executable_dir = DetectExecutableDir((argc > 0) ? argv[0] : nullptr);

  RUN_TEST(TestFrobeniusBenchCommit_HelpTextIsStable);
  RUN_TEST(TestFrobeniusBenchEval_HelpTextIsStable);
  RUN_TEST(TestFrobeniusBenchOuterCommit_HelpTextIsStable);
  RUN_TEST(TestFrobeniusBenchOuterProve_HelpTextIsStable);
  RUN_TEST(TestFrobeniusBenchOuterVerify_HelpTextIsStable);
  RUN_TEST(TestFrobeniusBenchCommit_SmokeRunPrintsStableFields);
  RUN_TEST(TestFrobeniusBenchCommit_HardcodedPresetForGR2P16_64);
  RUN_TEST(TestFrobeniusBenchOuterCommit_SmokeRunPrintsStableFields);
  RUN_TEST(TestFrobeniusBenchEval_SmokeRunPrintsStableFieldsAndSizes);
  RUN_TEST(TestFrobeniusBenchEval_HardcodedPresetForGR2P32_64);
  RUN_TEST(TestFrobeniusBenchOuterProve_SmokeRunPrintsStableFieldsAndSizes);
  RUN_TEST(TestFrobeniusBenchOuterVerify_SmokeRunPrintsStableFieldsAndSizes);

  if (g_test_failure_count != 0) {
    std::cerr << "\n" << g_test_failure_count << " test assertion(s) failed.\n";
    return 1;
  }

  std::cout << "\nAll tests passed.\n";
  return 0;
}
