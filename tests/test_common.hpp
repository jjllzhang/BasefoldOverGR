#ifndef TESTS_TEST_COMMON_HPP_
#define TESTS_TEST_COMMON_HPP_

#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pX.h>

#include <iostream>
#include <string>

extern int g_test_failure_count;

inline void Check(bool condition, const std::string &message, const char *file,
                  int line) {
  if (condition)
    return;
  std::cerr << file << ":" << line << " FAIL: " << message << "\n";
  g_test_failure_count++;
}

#define CHECK(cond) Check((cond), #cond, __FILE__, __LINE__)
#define CHECK_MSG(cond, msg) Check((cond), (msg), __FILE__, __LINE__)
#define CHECK_EQ(a, b)                                                         \
  Check(((a) == (b)), std::string(#a) + " == " + std::string(#b), __FILE__,    \
        __LINE__)
#define CHECK_LE(a, b)                                                         \
  CHECK_MSG(((a) <= (b)), std::string(#a) + " <= " + std::string(#b))
#define CHECK_LT(a, b)                                                         \
  CHECK_MSG(((a) < (b)), std::string(#a) + " < " + std::string(#b))
#define CHECK_GE(a, b)                                                         \
  CHECK_MSG(((a) >= (b)), std::string(#a) + " >= " + std::string(#b))
#define CHECK_GT(a, b)                                                         \
  CHECK_MSG(((a) > (b)), std::string(#a) + " > " + std::string(#b))

namespace testutil {

inline void PrintTestStart(const char *name) {
  std::cout << "\n[ RUN      ] " << name << "\n";
}

inline void PrintTestResult(const char *name, int new_failures) {
  if (new_failures == 0) {
    std::cout << "[       OK ] " << name << "\n";
  } else {
    std::cout << "[  FAILED  ] " << name << " (" << new_failures
              << " failure(s))\n";
  }
}

inline void PrintInfo(const std::string &message) {
  std::cout << "  " << message << "\n";
}

inline NTL::ZZ_pE ConstZZpE(long value) {
  NTL::ZZ_pX poly;
  SetCoeff(poly, 0, NTL::to_ZZ_p(value));
  NTL::ZZ_pE out;
  conv(out, poly);
  return out;
}

} // namespace testutil

#define RUN_TEST(fn)                                                           \
  do {                                                                         \
    testutil::PrintTestStart(#fn);                                             \
    const int _before = g_test_failure_count;                                  \
    fn();                                                                      \
    const int _after = g_test_failure_count;                                   \
    testutil::PrintTestResult(#fn, _after - _before);                          \
  } while (0)

#endif // TESTS_TEST_COMMON_HPP_
