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

constexpr const char *kRingF64 =
    "1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,1,0,"
    "0,0,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,0,1,1,0,0,0,0,0,1,0,1";
constexpr const char *kRingF128 =
    "1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
    "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
    "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
    "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
    "1";

fs::path g_executable_dir;

struct CommandResult {
  int exit_code = -1;
  bool exited = false;
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
            std::string("RunCommandCapture: pipe() failed for ") +
                exe.string());
  if (g_test_failure_count != 0) {
    return result;
  }

  std::cout.flush();
  std::cerr.flush();
  const pid_t pid = fork();
  CHECK_MSG(pid >= 0, std::string("RunCommandCapture: fork() failed for ") +
                          exe.string());
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
            std::string("RunCommandCapture: waitpid() failed for ") +
                exe.string());
  if (g_test_failure_count != 0) {
    return result;
  }

  result.exited = WIFEXITED(status);
  if (result.exited) {
    result.exit_code = WEXITSTATUS(status);
  }
#else
  (void)exe;
  (void)args;
#endif
  return result;
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
  testutil::PrintInfo(label +
                      ": skipped subprocess assertion on this platform");
#endif
}

void TestRingSwitchBenchCommit_HelpTextIsStable() {
  testutil::PrintInfo(
      "Ring-switch bench CLI: commit --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_commit";
  ExpectCommandSuccessContains(exe, {"--help"},
                               {"bench_z2k_ring_switch_commit",
                                "--basis-mode <default|provided>",
                                "hardcoded bench presets"},
                               "TestRingSwitchBenchCommit_HelpTextIsStable");
}

void TestRingSwitchBenchEval_HelpTextIsStable() {
  testutil::PrintInfo(
      "Ring-switch bench CLI: eval --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_eval";
  ExpectCommandSuccessContains(exe, {"--help"},
                               {"bench_z2k_ring_switch_eval", "--queries <int>",
                                "--checked", "--profiled-backend-verify",
                                "unchecked prover hot path",
                                "hardcoded bench presets"},
                               "TestRingSwitchBenchEval_HelpTextIsStable");
}

void TestRingSwitchBenchBackendEval_HelpTextIsStable() {
  testutil::PrintInfo(
      "Ring-switch bench CLI: backend eval --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_backend_eval";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_ring_switch_backend_eval", "--queries <int>",
       "packed backend PCS"},
      "TestRingSwitchBenchBackendEval_HelpTextIsStable");
}

void TestRingSwitchBenchOuterCommit_HelpTextIsStable() {
  testutil::PrintInfo(
      "Ring-switch bench CLI: outer commit --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_outer_commit";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_ring_switch_outer_commit", "--basis-mode <default|provided>",
       "Headline time measures", "No backend PCS commit"},
      "TestRingSwitchBenchOuterCommit_HelpTextIsStable");
}

void TestRingSwitchBenchOuterProve_HelpTextIsStable() {
  testutil::PrintInfo(
      "Ring-switch bench CLI: outer prove --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_outer_prove";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_ring_switch_outer_prove", "--queries <int>", "--checked",
       "unchecked prover hot path", "compiler-side outer proof generation"},
      "TestRingSwitchBenchOuterProve_HelpTextIsStable");
}

void TestRingSwitchBenchCommit_DefaultPresetForGR2p16r64() {
  testutil::PrintInfo("Ring-switch bench CLI: commit smoke run picks the "
                      "GR(2^16,64) commit-like preset");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_commit";
  ExpectCommandSuccessContains(
      exe,
      {"--ell", "7", "--kappa", "6", "--warmup", "0", "--reps", "1",
       "--ring-mod", "65536", "--ring-p", "2", "--ring-F", kRingF64,
       "--ring-zeta", "0,1"},
      {"[ring-switch commit]",
       "basis mode default (bench-fixed commit-like preset)",
       "default alpha preset poly", "default beta preset lower_16",
       "packing mean", "commit  mean"},
      "TestRingSwitchBenchCommit_DefaultPresetForGR2p16r64");
}

void TestRingSwitchBenchOuterCommit_DefaultPresetForGR2p16r64() {
  testutil::PrintInfo("Ring-switch bench CLI: outer commit smoke run prints "
                      "stable result fields");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_outer_commit";
  ExpectCommandSuccessContains(
      exe,
      {"--ell", "7", "--kappa", "6", "--warmup", "0", "--reps", "1",
       "--ring-mod", "65536", "--ring-p", "2", "--ring-F", kRingF64,
       "--ring-zeta", "0,1"},
      {"[ring-switch outer commit]",
       "basis mode default (bench-fixed commit-like preset)",
       "default alpha preset poly", "default beta preset lower_16",
       "outer commit mean", "backend input size"},
      "TestRingSwitchBenchOuterCommit_DefaultPresetForGR2p16r64");
}

void TestRingSwitchBenchCommit_DefaultPresetForGR2p16r128() {
  testutil::PrintInfo("Ring-switch bench CLI: commit smoke run picks the "
                      "GR(2^16,128) commit-like preset");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_commit";
  ExpectCommandSuccessContains(
      exe,
      {"--ell", "8", "--kappa", "7", "--warmup", "0", "--reps", "1",
       "--ring-mod", "65536", "--ring-p", "2", "--ring-F", kRingF128,
       "--ring-zeta", "0,1"},
      {"[ring-switch commit]",
       "basis mode default (bench-fixed commit-like preset)",
       "default alpha preset poly", "default beta preset lower_64",
       "packing mean", "commit  mean"},
      "TestRingSwitchBenchCommit_DefaultPresetForGR2p16r128");
}

void TestRingSwitchBenchEval_DefaultPresetForGR2p16r128() {
  testutil::PrintInfo("Ring-switch bench CLI: eval smoke run picks the "
                      "GR(2^16,128) eval-like preset");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_eval";
  ExpectCommandSuccessContains(
      exe,
      {"--ell", "8", "--kappa", "7", "--queries", "2", "--warmup", "0",
       "--reps", "1", "--ring-mod", "65536", "--ring-p", "2", "--ring-F",
       kRingF128, "--ring-zeta", "0,1"},
      {"[ring-switch eval]",
       "basis mode default (bench-fixed eval-like preset)",
       "default alpha preset poly", "default beta preset lower_64",
       "outer commit mean", "backend commit mean", "commit total mean",
       "prove-phase mean", "verifier mean"},
      "TestRingSwitchBenchEval_DefaultPresetForGR2p16r128");
}

void TestRingSwitchBenchEval_DefaultPresetForGR2p32r64() {
  testutil::PrintInfo("Ring-switch bench CLI: eval smoke run picks the "
                      "GR(2^32,64) eval-like preset");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_eval";
  ExpectCommandSuccessContains(
      exe,
      {"--ell", "7", "--kappa", "6", "--queries", "2", "--warmup", "0",
       "--reps", "1", "--ring-mod", "4294967296", "--ring-p", "2", "--ring-F",
       kRingF64, "--ring-zeta", "0,1"},
      {"[ring-switch eval]",
       "basis mode default (bench-fixed eval-like preset)",
       "default alpha preset poly_dual", "default beta preset lower_16",
       "outer commit mean", "backend commit mean", "commit total mean",
       "prove-phase mean", "verifier mean"},
      "TestRingSwitchBenchEval_DefaultPresetForGR2p32r64");
}

void TestRingSwitchBenchBackendEval_DefaultPresetForGR2p32r64() {
  testutil::PrintInfo("Ring-switch bench CLI: backend eval smoke run prints "
                      "stable result fields");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_backend_eval";
  ExpectCommandSuccessContains(
      exe,
      {"--ell", "7", "--kappa", "6", "--queries", "2", "--warmup", "0",
       "--reps", "1", "--ring-mod", "4294967296", "--ring-p", "2", "--ring-F",
       kRingF64, "--ring-zeta", "0,1"},
      {"[ring-switch backend eval]",
       "basis mode default (bench-fixed eval-like preset)",
       "backend commit mean", "prove-phase mean", "verifier mean",
       "proof size"},
      "TestRingSwitchBenchBackendEval_DefaultPresetForGR2p32r64");
}

void TestRingSwitchBenchEval_CheckedFlagAccepted() {
  testutil::PrintInfo("Ring-switch bench CLI: eval accepts --checked and still "
                      "prints stable fields");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_eval";
  ExpectCommandSuccessContains(
      exe, {"--checked", "--warmup", "0", "--reps", "1", "--queries", "2"},
      {"[ring-switch eval]", "prove-phase mean", "verifier mean"},
      "TestRingSwitchBenchEval_CheckedFlagAccepted");
}

void TestRingSwitchBenchEval_ProfiledBackendVerifyFlagAccepted() {
  testutil::PrintInfo(
      "Ring-switch bench CLI: eval accepts --profiled-backend-verify and still "
      "prints stable fields");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_eval";
  ExpectCommandSuccessContains(
      exe,
      {"--profiled-backend-verify", "--warmup", "0", "--reps", "1", "--queries",
       "2"},
      {"[ring-switch eval]", "backend verify timing profiled-subcall",
       "backend verifier mean"},
      "TestRingSwitchBenchEval_ProfiledBackendVerifyFlagAccepted");
}

void TestRingSwitchBenchOuterProve_CheckedFlagAccepted() {
  testutil::PrintInfo("Ring-switch bench CLI: outer prove accepts --checked "
                      "and still prints stable fields");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_outer_prove";
  ExpectCommandSuccessContains(
      exe, {"--checked", "--warmup", "0", "--reps", "1", "--queries", "2"},
      {"[ring-switch outer prove]", "prove-phase mean", "proof size"},
      "TestRingSwitchBenchOuterProve_CheckedFlagAccepted");
}

void TestRingSwitchBenchOuterProve_ProfileFlagAccepted() {
  testutil::PrintInfo("Ring-switch bench CLI: outer prove accepts --profile "
                      "and prints detailed profile buckets");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_outer_prove";
  ExpectCommandSuccessContains(
      exe, {"--profile", "--warmup", "0", "--reps", "1", "--queries", "2"},
      {"[ring-switch outer prove]", "[profile-ring-switch-outer-prover]",
       "BuildComponentTensor", "BuildSuffixEqTable", "RecoverAlphaCoords",
       "LiftAlphaRows", "ComputeSByU", "RecoverPartialsEq1", "BatchPrepPrefix",
       "TranscriptAbsorb", "PrefixChallengesRPrime", "BuildEqPrefixTable",
       "InitialClaim", "BuildBatchedGTable", "SumcheckSuffixRounds",
       "ProductSumcheckInit", "CurrentPolynomial", "ReceiveChallenge",
       "FinalTstarGstarEq3"},
      "TestRingSwitchBenchOuterProve_ProfileFlagAccepted");
}

void TestRingSwitchBenchOuterVerify_ProfileFlagAccepted() {
  testutil::PrintInfo("Ring-switch bench CLI: outer verify accepts --profile "
                      "and prints detailed profile buckets");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_outer_verify";
  ExpectCommandSuccessContains(
      exe, {"--profile", "--warmup", "0", "--reps", "1", "--queries", "2"},
      {"[ring-switch outer verify]", "[profile-ring-switch-outer-verifier]",
       "BuildComponentTensor", "BuildSuffixEqTable", "RecoverAlphaCoords",
       "LiftAlphaRows", "RecoverPartialsEq1", "ReplayPrefixBatching",
       "TranscriptAbsorb", "PrefixChallengesRPrime", "BuildEqPrefixTable",
       "InitialClaim", "ReplaySumcheckChain", "FinalGstarEq3",
       "FoldSuffixPoint", "EvaluateGstar"},
      "TestRingSwitchBenchOuterVerify_ProfileFlagAccepted");
}

} // namespace

int main(int argc, char **argv) {
  g_executable_dir = DetectExecutableDir(argc > 0 ? argv[0] : nullptr);
  TestRingSwitchBenchCommit_HelpTextIsStable();
  TestRingSwitchBenchEval_HelpTextIsStable();
  TestRingSwitchBenchBackendEval_HelpTextIsStable();
  TestRingSwitchBenchOuterCommit_HelpTextIsStable();
  TestRingSwitchBenchOuterProve_HelpTextIsStable();
  TestRingSwitchBenchCommit_DefaultPresetForGR2p16r64();
  TestRingSwitchBenchOuterCommit_DefaultPresetForGR2p16r64();
  TestRingSwitchBenchCommit_DefaultPresetForGR2p16r128();
  TestRingSwitchBenchEval_DefaultPresetForGR2p32r64();
  TestRingSwitchBenchEval_DefaultPresetForGR2p16r128();
  TestRingSwitchBenchBackendEval_DefaultPresetForGR2p32r64();
  TestRingSwitchBenchEval_CheckedFlagAccepted();
  TestRingSwitchBenchEval_ProfiledBackendVerifyFlagAccepted();
  TestRingSwitchBenchOuterProve_CheckedFlagAccepted();
  TestRingSwitchBenchOuterProve_ProfileFlagAccepted();
  TestRingSwitchBenchOuterVerify_ProfileFlagAccepted();

  if (g_test_failure_count != 0) {
    std::cerr << g_test_failure_count
              << " ring-switch bench CLI test(s) failed\n";
    return 1;
  }
  std::cout << "All ring-switch bench CLI tests passed\n";
  return 0;
}
