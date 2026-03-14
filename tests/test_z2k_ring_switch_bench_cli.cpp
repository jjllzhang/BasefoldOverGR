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
    "1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0,1,0,0,1,1,1,0,1,0,1,0,0,1,1,0,0,0,0,0,1,0,1";

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
  testutil::PrintInfo(label + ": skipped subprocess assertion on this platform");
#endif
}

void TestRingSwitchBenchCommit_HelpTextIsStable() {
  testutil::PrintInfo("Ring-switch bench CLI: commit --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_commit";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_ring_switch_commit", "--basis-mode <default|provided>",
       "hardcoded bench presets"},
      "TestRingSwitchBenchCommit_HelpTextIsStable");
}

void TestRingSwitchBenchEval_HelpTextIsStable() {
  testutil::PrintInfo("Ring-switch bench CLI: eval --help prints stable usage anchors");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_eval";
  ExpectCommandSuccessContains(
      exe, {"--help"},
      {"bench_z2k_ring_switch_eval", "--queries <int>",
       "hardcoded bench presets"},
      "TestRingSwitchBenchEval_HelpTextIsStable");
}

void TestRingSwitchBenchCommit_DefaultPresetForGR2p16r64() {
  testutil::PrintInfo(
      "Ring-switch bench CLI: commit smoke run picks the GR(2^16,64) commit-like preset");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_commit";
  ExpectCommandSuccessContains(
      exe,
      {"--ell", "7", "--kappa", "6", "--warmup", "0", "--reps", "1",
       "--ring-mod", "65536", "--ring-p", "2", "--ring-F", kRingF64,
       "--ring-zeta", "0,1"},
      {"[ring-switch commit]", "basis mode default (bench-fixed commit-like preset)",
       "default alpha preset poly", "default beta preset lower_16",
       "packing mean", "commit  mean"},
      "TestRingSwitchBenchCommit_DefaultPresetForGR2p16r64");
}

void TestRingSwitchBenchEval_DefaultPresetForGR2p32r64() {
  testutil::PrintInfo(
      "Ring-switch bench CLI: eval smoke run picks the GR(2^32,64) eval-like preset");

  const fs::path exe = g_executable_dir / "bench_z2k_ring_switch_eval";
  ExpectCommandSuccessContains(
      exe,
      {"--ell", "7", "--kappa", "6", "--queries", "2", "--warmup", "0",
       "--reps", "1", "--ring-mod", "4294967296", "--ring-p", "2",
       "--ring-F", kRingF64, "--ring-zeta", "0,1"},
      {"[ring-switch eval]", "basis mode default (bench-fixed eval-like preset)",
       "default alpha preset poly_dual", "default beta preset lower_16",
       "prove-phase mean", "verifier mean"},
      "TestRingSwitchBenchEval_DefaultPresetForGR2p32r64");
}

}  // namespace

int main(int argc, char **argv) {
  g_executable_dir = DetectExecutableDir(argc > 0 ? argv[0] : nullptr);
  TestRingSwitchBenchCommit_HelpTextIsStable();
  TestRingSwitchBenchEval_HelpTextIsStable();
  TestRingSwitchBenchCommit_DefaultPresetForGR2p16r64();
  TestRingSwitchBenchEval_DefaultPresetForGR2p32r64();

  if (g_test_failure_count != 0) {
    std::cerr << g_test_failure_count << " ring-switch bench CLI test(s) failed\n";
    return 1;
  }
  std::cout << "All ring-switch bench CLI tests passed\n";
  return 0;
}
