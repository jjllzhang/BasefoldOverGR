#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct CliArgs {
  long long d = -1;                   // multilinear variable dimension
  long long c = -1;                   // inverse code rate
  long long k0 = 1;                   // base message dimension
  long long lambda = -1;              // target security bits
  std::optional<long double> q;       // randomness domain size
  std::optional<long double> p;       // prime base for q = p^(r*m)
  std::optional<long long> r;         // residue degree
  long long m = 1;                    // extension multiplier
  std::optional<long double> gamma;   // slack parameter
  bool auto_gamma = false;
  std::optional<long double> gamma_min;
  std::optional<long double> gamma_max;
  long long gamma_steps = 4000;
  bool show_levels = false;
  bool help = false;
};

struct LevelInfo {
  long long i = 0;
  long double n_prev = 0;
  long double ell_i = 0;
  long double t_i = 0;
};

struct CalcResult {
  long double q = 0;
  long double gamma = 0;
  bool gamma_auto = false;
  bool gamma_auto_feasible = true;
  bool base_valid = true;
  bool delta_johnson_constraint_valid = true;
  bool delta_cd_constraint_valid = true;
  long double gamma_search_min = 0;
  long double gamma_search_max = 0;
  long long gamma_search_steps = 0;  // user-provided search budget
  long long gamma_search_evals = 0;  // actual gamma evaluations
  long double k_d = 0;
  long double n_d = 0;
  long double t_d = 0;
  long double delta_code_lower = 0;
  long double delta_upper_johnson = 0;
  long double delta = 0;
  long double base = 0;
  bool has_finite_l_iopp = false;
  bool has_finite_l_pcs = false;
  std::uint64_t l_min_iopp = 0;
  std::uint64_t l_min_pcs = 0;  // recommended: satisfies PCS and thus IOPP
  long double target_2_minus_lambda = 0;
  long double iopp_budget = 0;   // 2^-lambda - 2d/(gamma^3*q)
  long double pcs_budget = 0;    // 2^-lambda - 2d/(gamma^3*q) - 2d/q
  long double second_term_iopp = 0;
  long double second_term_pcs = 0;
  long double iopp_first_term = 0;
  long double sumcheck_term = 0;
  long double iopp_bound_at_l_pcs = 0;
  long double pcs_bound_at_l_pcs = 0;
  std::vector<LevelInfo> levels;
};

struct LSolveResult {
  bool has_finite_l = false;
  std::uint64_t l_min = 0;
  long double base_pow_l = 0;
};

struct PrecomputedCodeParams {
  long double k_d = 0;
  long double n_d = 0;
  long double t_d = 0;
  long double delta_code_lower = 0;
  std::vector<LevelInfo> levels;
};

void PrintUsage(const char* argv0) {
  std::cout
      << "Usage:\n"
      << "  " << argv0 << " --d <int> --c <int> --lambda <int> [options]\n\n"
      << "Required:\n"
      << "  --d <int>           Multilinear polynomial dimension (variables)\n"
      << "  --c <int>           Inverse code rate c (c >= 2)\n"
      << "  --lambda <int>      Security parameter lambda (bits)\n\n"
      << "Randomness domain q (choose one mode):\n"
      << "  --q <real>          Directly set q\n"
      << "  --p <real> --r <int> [--m <int>]\n"
      << "                      Set q = p^(r*m) (ring/extension style)\n\n"
      << "Optional:\n"
      << "  --k0 <int>          Base message dimension k0 (default: 1)\n"
      << "  --gamma <real>      Slack gamma in (0,1) (default: 1/(10*max(1,d)))\n"
      << "  --auto-gamma        Search gamma automatically to minimize l_min_for_PCS\n"
      << "  --gamma-min <real>  Auto-search lower bound (default: 1e-9)\n"
      << "  --gamma-max <real>  Auto-search upper bound (default: 0.99)\n"
      << "  --gamma-steps <int> Auto-search budget (default: 4000)\n"
      << "  --show-levels       Print per-level (n_{i-1}, ell_i, t_i)\n"
      << "  --help              Show this message\n\n"
      << "Formulas (from Basefold_over_GR.pdf, parameter selection):\n"
      << "  t0 = k0\n"
      << "  t_i = 2*t_{i-1} + ell_i\n"
      << "  ell_i = [ (2*log2(q/(q-1)) + 2)*t_{i-1} + log2(9/4)*n_{i-1} + lambda ] / (log2(q-1)-1)\n"
      << "  Delta_Cd >= 1 - t_d / n_d\n"
      << "  J_gamma(x) = 1 - sqrt(1 - x*(1-gamma))\n"
      << "  delta < J_gamma(J_gamma(Delta_Cd))\n"
      << "  base = 1 - delta + gamma*d\n"
      << "  validity requires: delta < J_gamma(J_gamma(Delta_Cd)), 0 < base < 1, and 3*delta - gamma*d < Delta_Cd\n"
      << "  l_iopp from: 2d/(gamma^3*q) + base^l <= 2^-lambda\n"
      << "  l_pcs  from: 2d/q + 2d/(gamma^3*q) + base^l <= 2^-lambda\n"
      << "  recommended l = l_pcs\n";
}

long long ParseInt(const std::string& s, const char* name) {
  try {
    std::size_t idx = 0;
    const long long v = std::stoll(s, &idx);
    if (idx != s.size()) {
      throw std::runtime_error("");
    }
    return v;
  } catch (...) {
    throw std::runtime_error(std::string("Invalid integer for ") + name + ": " + s);
  }
}

long double ParseReal(const std::string& s, const char* name) {
  try {
    std::size_t idx = 0;
    const long double v = std::stold(s, &idx);
    if (idx != s.size()) {
      throw std::runtime_error("");
    }
    return v;
  } catch (...) {
    throw std::runtime_error(std::string("Invalid real for ") + name + ": " + s);
  }
}

long double Pow2LD(long long exp) {
  if (exp < 0) {
    throw std::runtime_error("Exponent must be non-negative");
  }
  long double v = 1.0L;
  for (long long i = 0; i < exp; ++i) {
    v *= 2.0L;
    if (!std::isfinite(v)) {
      throw std::runtime_error("Overflow while computing 2^d");
    }
  }
  return v;
}

long double PowIntLD(long double base, long long exp) {
  if (exp < 0) {
    throw std::runtime_error("Exponent must be non-negative");
  }
  long double v = 1.0L;
  for (long long i = 0; i < exp; ++i) {
    v *= base;
    if (!std::isfinite(v)) {
      throw std::runtime_error("Overflow while computing q = p^(r*m)");
    }
  }
  return v;
}

long double Clamp01(long double x) {
  if (x < 0.0L) return 0.0L;
  if (x > 1.0L) return 1.0L;
  return x;
}

long double Johnson(long double x, long double gamma) {
  const long double inner = 1.0L - x * (1.0L - gamma);
  if (inner < 0.0L) {
    throw std::runtime_error("Johnson bound sqrt argument is negative; check inputs");
  }
  return 1.0L - std::sqrt(inner);
}

CliArgs ParseArgs(int argc, char** argv) {
  CliArgs args;
  for (int i = 1; i < argc; ++i) {
    const std::string a(argv[i]);
    auto need_value = [&](const char* name) -> std::string {
      if (i + 1 >= argc) {
        throw std::runtime_error(std::string("Missing value for ") + name);
      }
      ++i;
      return std::string(argv[i]);
    };

    if (a == "--help" || a == "-h") {
      args.help = true;
    } else if (a == "--d") {
      args.d = ParseInt(need_value("--d"), "--d");
    } else if (a == "--c") {
      args.c = ParseInt(need_value("--c"), "--c");
    } else if (a == "--k0") {
      args.k0 = ParseInt(need_value("--k0"), "--k0");
    } else if (a == "--lambda") {
      args.lambda = ParseInt(need_value("--lambda"), "--lambda");
    } else if (a == "--gamma") {
      args.gamma = ParseReal(need_value("--gamma"), "--gamma");
    } else if (a == "--auto-gamma") {
      args.auto_gamma = true;
    } else if (a == "--gamma-min") {
      args.gamma_min = ParseReal(need_value("--gamma-min"), "--gamma-min");
    } else if (a == "--gamma-max") {
      args.gamma_max = ParseReal(need_value("--gamma-max"), "--gamma-max");
    } else if (a == "--gamma-steps") {
      args.gamma_steps = ParseInt(need_value("--gamma-steps"), "--gamma-steps");
    } else if (a == "--q") {
      args.q = ParseReal(need_value("--q"), "--q");
    } else if (a == "--p") {
      args.p = ParseReal(need_value("--p"), "--p");
    } else if (a == "--r") {
      args.r = ParseInt(need_value("--r"), "--r");
    } else if (a == "--m") {
      args.m = ParseInt(need_value("--m"), "--m");
    } else if (a == "--show-levels") {
      args.show_levels = true;
    } else {
      throw std::runtime_error("Unknown argument: " + a);
    }
  }
  return args;
}

LSolveResult SolveMinimalL(long double base, long double residual_budget) {
  LSolveResult out;
  if (!(residual_budget > 0.0L)) {
    return out;  // no finite l
  }
  if (base <= 0.0L) {
    out.has_finite_l = true;
    out.l_min = 1;
    out.base_pow_l = 0.0L;
    return out;
  }
  if (base >= 1.0L) {
    return out;  // no finite l because base^l does not decay
  }
  if (residual_budget >= 1.0L) {
    out.has_finite_l = true;
    out.l_min = 1;
    out.base_pow_l = base;
    return out;
  }

  const long double log_base = std::log2(base);  // < 0
  const long double l_real = std::ceil(std::log2(residual_budget) / log_base);
  if (!std::isfinite(l_real) || l_real > static_cast<long double>(
                                       std::numeric_limits<std::uint64_t>::max())) {
    throw std::runtime_error("Computed l is out of uint64 range");
  }
  std::uint64_t l = (l_real < 1.0L) ? 1 : static_cast<std::uint64_t>(l_real);

  long double term = std::pow(base, static_cast<long double>(l));
  while (term > residual_budget &&
         l < std::numeric_limits<std::uint64_t>::max()) {
    ++l;
    term = std::pow(base, static_cast<long double>(l));
  }
  while (l > 1) {
    const long double prev = std::pow(base, static_cast<long double>(l - 1));
    if (prev <= residual_budget) {
      --l;
      term = prev;
    } else {
      break;
    }
  }

  out.has_finite_l = true;
  out.l_min = l;
  out.base_pow_l = term;
  return out;
}

void ValidateBasicArgs(const CliArgs& args) {
  if (args.d < 0) throw std::runtime_error("--d must be >= 0");
  if (args.c < 2) throw std::runtime_error("--c must be >= 2");
  if (args.k0 <= 0) throw std::runtime_error("--k0 must be > 0");
  if (args.lambda <= 0) throw std::runtime_error("--lambda must be > 0");
  if (args.m <= 0) throw std::runtime_error("--m must be > 0");
  if (args.gamma_steps <= 1) {
    throw std::runtime_error("--gamma-steps must be > 1");
  }
  if (args.auto_gamma && args.gamma.has_value()) {
    throw std::runtime_error("Use either --gamma or --auto-gamma, not both");
  }
}

long double ResolveQ(const CliArgs& args) {
  if (args.q.has_value() && (args.p.has_value() || args.r.has_value())) {
    throw std::runtime_error("Use either --q OR (--p,--r,[--m]), not both");
  }

  if (args.q.has_value()) {
    if (!(*args.q > 3.0L)) {
      throw std::runtime_error("q must satisfy q > 3 (equivalently |F*| > 2)");
    }
    return *args.q;
  } else {
    if (!args.p.has_value() || !args.r.has_value()) {
      throw std::runtime_error(
          "Missing q: provide --q or provide --p and --r (optionally --m)");
    }
    if (*args.p <= 1.0L) throw std::runtime_error("--p must be > 1");
    if (*args.r <= 0) throw std::runtime_error("--r must be > 0");
    const long long exp = (*args.r) * args.m;
    if (exp <= 0) throw std::runtime_error("r*m must be > 0");
    const long double q = PowIntLD(*args.p, exp);
    if (!(q > 3.0L)) {
      throw std::runtime_error("q must satisfy q > 3 (equivalently |F*| > 2)");
    }
    return q;
  }
}

std::vector<long double> BuildLogGrid(long double g_min, long double g_max,
                                      long long points) {
  if (points < 2) {
    throw std::runtime_error("BuildLogGrid requires at least 2 points");
  }
  std::vector<long double> grid;
  grid.reserve(static_cast<std::size_t>(points));
  const long double log_min = std::log(g_min);
  const long double log_max = std::log(g_max);
  for (long long i = 0; i < points; ++i) {
    const long double t = static_cast<long double>(i) /
                          static_cast<long double>(points - 1);
    grid.push_back(std::exp(log_min + t * (log_max - log_min)));
  }
  return grid;
}

PrecomputedCodeParams PrecomputeCode(const CliArgs& args, long double q,
                                     bool keep_levels) {
  PrecomputedCodeParams out;
  const long double c_ld = static_cast<long double>(args.c);
  const long double k0_ld = static_cast<long double>(args.k0);
  const long double lambda_ld = static_cast<long double>(args.lambda);

  const long double log_qm1 = std::log2(q - 1.0L);
  const long double denom = log_qm1 - 1.0L;
  if (!(denom > 0.0L)) {
    throw std::runtime_error("Need log2(q-1) - 1 > 0; this requires q > 3");
  }

  const long double coeff_t = 2.0L * std::log2(q / (q - 1.0L)) + 2.0L;
  const long double coeff_n = std::log2(9.0L / 4.0L);

  long double t = k0_ld;
  long double n_prev = c_ld * k0_ld;  // n_0
  if (!std::isfinite(t) || !std::isfinite(n_prev)) {
    throw std::runtime_error("Numeric overflow before recurrence");
  }

  out.levels.clear();
  if (keep_levels) out.levels.reserve(static_cast<std::size_t>(args.d));
  for (long long i = 1; i <= args.d; ++i) {
    const long double ell_i =
        ((coeff_t * t) + (coeff_n * n_prev) + lambda_ld) / denom;
    t = 2.0L * t + ell_i;
    if (!std::isfinite(ell_i) || !std::isfinite(t)) {
      throw std::runtime_error("Numeric overflow in recurrence");
    }
    if (keep_levels) out.levels.push_back(LevelInfo{i, n_prev, ell_i, t});
    n_prev *= 2.0L;
  }

  out.k_d = k0_ld * Pow2LD(args.d);
  out.n_d = c_ld * out.k_d;
  out.t_d = t;
  out.delta_code_lower = Clamp01(1.0L - out.t_d / out.n_d);
  return out;
}

bool IsValidPoint(const CalcResult& out) {
  return out.base_valid && out.delta_johnson_constraint_valid &&
         out.delta_cd_constraint_valid;
}

bool BetterAny(const CalcResult& lhs, long double lhs_gamma,
               const CalcResult& rhs, long double rhs_gamma) {
  if (lhs.pcs_budget != rhs.pcs_budget) return lhs.pcs_budget > rhs.pcs_budget;
  if (lhs.base != rhs.base) return lhs.base < rhs.base;
  return lhs_gamma < rhs_gamma;
}

bool BetterFeasible(const CalcResult& lhs, long double lhs_gamma,
                    const CalcResult& rhs, long double rhs_gamma) {
  if (lhs.l_min_pcs != rhs.l_min_pcs) return lhs.l_min_pcs < rhs.l_min_pcs;
  if (lhs.pcs_bound_at_l_pcs != rhs.pcs_bound_at_l_pcs) {
    return lhs.pcs_bound_at_l_pcs < rhs.pcs_bound_at_l_pcs;
  }
  return lhs_gamma < rhs_gamma;
}

CalcResult ComputeAtGamma(const CliArgs& args, const PrecomputedCodeParams& pre,
                          long double q, long double gamma, bool keep_levels,
                          bool enforce_validity = true) {
  if (!(gamma > 0.0L && gamma < 1.0L)) {
    throw std::runtime_error("gamma must be in (0,1)");
  }
  CalcResult out;
  out.q = q;
  out.gamma = gamma;

  const long double d_ld = static_cast<long double>(args.d);
  const long double lambda_ld = static_cast<long double>(args.lambda);
  out.k_d = pre.k_d;
  out.n_d = pre.n_d;
  out.t_d = pre.t_d;
  out.delta_code_lower = pre.delta_code_lower;
  out.levels = keep_levels ? pre.levels : std::vector<LevelInfo>{};
  const long double j1 = Johnson(out.delta_code_lower, out.gamma);
  out.delta_upper_johnson = Clamp01(Johnson(j1, out.gamma));

  const long double delta_lb_from_base = out.gamma * d_ld;  // strict: base < 1
  const long double delta_ub_from_delta_cd =
      (out.delta_code_lower + out.gamma * d_ld) /
      3.0L;  // strict: 3*delta - gamma*d < Delta_Cd
  const long double delta_ub_from_base = 1.0L + out.gamma * d_ld;  // strict: base > 0
  const long double delta_ub =
      std::min(std::min(out.delta_upper_johnson, delta_ub_from_delta_cd),
               std::min(delta_ub_from_base, 1.0L));

  if (!(delta_ub > delta_lb_from_base)) {
    out.delta = delta_lb_from_base;
    out.base = 1.0L - out.delta + out.gamma * d_ld;
    out.base_valid = false;
    out.delta_johnson_constraint_valid = false;
    out.delta_cd_constraint_valid = false;
    if (enforce_validity) {
      throw std::runtime_error(
          "Invalid parameters: empty feasible interval for delta under "
          "delta < J_gamma(J_gamma(Delta_Cd)), 0 < base < 1, and "
          "3*delta - gamma*d < Delta_Cd");
    }
    return out;
  }

  long double chosen_delta = std::nextafter(
      delta_ub, -std::numeric_limits<long double>::infinity());
  if (!(chosen_delta > delta_lb_from_base)) {
    chosen_delta = (delta_lb_from_base + delta_ub) / 2.0L;
  }
  out.delta = chosen_delta;

  out.base = 1.0L - out.delta + out.gamma * d_ld;
  out.base_valid = (out.base > 0.0L && out.base < 1.0L);
  out.delta_johnson_constraint_valid = (out.delta < out.delta_upper_johnson);
  out.delta_cd_constraint_valid =
      ((3.0L * out.delta - out.gamma * d_ld) < out.delta_code_lower);
  if (!out.base_valid || !out.delta_johnson_constraint_valid ||
      !out.delta_cd_constraint_valid) {
    if (enforce_validity) {
      if (!out.base_valid && !out.delta_johnson_constraint_valid &&
          !out.delta_cd_constraint_valid) {
        throw std::runtime_error(
            "Invalid parameters: require delta < J_gamma(J_gamma(Delta_Cd)), "
            "0 < base < 1, and 3*delta - gamma*d < Delta_Cd");
      }
      if (!out.base_valid) {
        throw std::runtime_error(
            "Invalid parameters: base = 1 - delta + gamma*d must be in (0,1)");
      }
      if (!out.delta_johnson_constraint_valid) {
        throw std::runtime_error(
            "Invalid parameters: require delta < J_gamma(J_gamma(Delta_Cd))");
      }
      throw std::runtime_error(
          "Invalid parameters: require 3*delta - gamma*d < Delta_Cd");
    }
    return out;
  }
  out.target_2_minus_lambda = std::exp2(-lambda_ld);
  out.iopp_first_term =
      (2.0L * d_ld) / (std::pow(out.gamma, 3.0L) * q);
  out.sumcheck_term = (2.0L * d_ld) / q;

  out.iopp_budget = out.target_2_minus_lambda - out.iopp_first_term;
  out.pcs_budget = out.target_2_minus_lambda - out.iopp_first_term - out.sumcheck_term;

  const LSolveResult l_iopp = SolveMinimalL(out.base, out.iopp_budget);
  out.has_finite_l_iopp = l_iopp.has_finite_l;
  out.l_min_iopp = l_iopp.l_min;
  out.second_term_iopp = l_iopp.base_pow_l;

  const LSolveResult l_pcs = SolveMinimalL(out.base, out.pcs_budget);
  out.has_finite_l_pcs = l_pcs.has_finite_l;
  out.l_min_pcs = l_pcs.l_min;
  out.second_term_pcs = l_pcs.base_pow_l;

  if (out.has_finite_l_pcs) {
    out.iopp_bound_at_l_pcs = out.iopp_first_term + out.second_term_pcs;
    out.pcs_bound_at_l_pcs = out.sumcheck_term + out.iopp_bound_at_l_pcs;
  } else {
    out.iopp_bound_at_l_pcs = std::numeric_limits<long double>::quiet_NaN();
    out.pcs_bound_at_l_pcs = std::numeric_limits<long double>::quiet_NaN();
  }
  return out;
}

CalcResult Compute(const CliArgs& args) {
  ValidateBasicArgs(args);
  const long double q = ResolveQ(args);
  const PrecomputedCodeParams pre = PrecomputeCode(args, q, args.show_levels);

  if (!args.auto_gamma) {
    const long double gamma = args.gamma.value_or(
        1.0L / (10.0L * static_cast<long double>((args.d > 0) ? args.d : 1)));
    CalcResult out = ComputeAtGamma(args, pre, q, gamma, args.show_levels);
    out.gamma_auto = false;
    return out;
  }

  long double g_min = args.gamma_min.value_or(1e-9L);
  long double g_max = args.gamma_max.value_or(0.99L);
  if (!(g_min > 0.0L && g_min < 1.0L)) {
    throw std::runtime_error("--gamma-min must be in (0,1)");
  }
  if (!(g_max > 0.0L && g_max < 1.0L)) {
    throw std::runtime_error("--gamma-max must be in (0,1)");
  }
  if (!(g_min < g_max)) {
    throw std::runtime_error("Need gamma_min < gamma_max");
  }

  const long double d_ld = static_cast<long double>(args.d);
  const long double target = std::exp2(-static_cast<long double>(args.lambda));
  const long double sumcheck_term = (2.0L * d_ld) / q;
  const long double residual = target - sumcheck_term;
  if (residual > 0.0L) {
    const long double lb = std::cbrt((2.0L * d_ld) / (q * residual)) * 1.0000000001L;
    if (std::isfinite(lb) && lb > g_min && lb < 1.0L) {
      g_min = lb;
    }
  }
  if (!(g_min < g_max)) {
    g_min = args.gamma_min.value_or(1e-9L);
    g_max = args.gamma_max.value_or(0.99L);
  }

  bool have_feasible = false;
  CalcResult best_feasible;
  long double best_feasible_gamma = 0.0L;
  bool have_any = false;
  CalcResult best_any;
  long double best_any_gamma = 0.0L;
  long long eval_count = 0;

  auto consider = [&](const CalcResult& cur, long double gamma) {
    if (!IsValidPoint(cur)) return;
    if (!have_any || BetterAny(cur, gamma, best_any, best_any_gamma)) {
      best_any = cur;
      best_any_gamma = gamma;
      have_any = true;
    }
    if (cur.has_finite_l_pcs &&
        (!have_feasible ||
         BetterFeasible(cur, gamma, best_feasible, best_feasible_gamma))) {
      best_feasible = cur;
      best_feasible_gamma = gamma;
      have_feasible = true;
    }
  };

  struct GammaPoint {
    long double gamma = 0;
    CalcResult result;
  };

  const long long eval_budget = args.gamma_steps;
  long long coarse_points = static_cast<long long>(
      std::floor(std::sqrt(static_cast<long double>(eval_budget)) * 4.0L));
  coarse_points = std::max<long long>(17, coarse_points);
  coarse_points = std::min<long long>(coarse_points, 401);
  coarse_points = std::min<long long>(coarse_points, eval_budget);
  if (coarse_points < 2) coarse_points = 2;

  const std::vector<long double> coarse_grid =
      BuildLogGrid(g_min, g_max, coarse_points);
  std::vector<GammaPoint> coarse_samples;
  coarse_samples.reserve(static_cast<std::size_t>(coarse_points));
  for (const long double gamma : coarse_grid) {
    CalcResult cur = ComputeAtGamma(args, pre, q, gamma, false, false);
    ++eval_count;
    coarse_samples.push_back(GammaPoint{gamma, cur});
    consider(cur, gamma);
  }

  if (!have_any) {
    throw std::runtime_error(
        "auto-gamma search failed: no gamma in search range satisfies validity constraints");
  }

  std::vector<long long> ranked_indices;
  ranked_indices.reserve(static_cast<std::size_t>(coarse_samples.size()));
  const bool rank_feasible = have_feasible;
  for (long long i = 0; i < static_cast<long long>(coarse_samples.size()); ++i) {
    const CalcResult& cur = coarse_samples[static_cast<std::size_t>(i)].result;
    if (!IsValidPoint(cur)) continue;
    if (rank_feasible && !cur.has_finite_l_pcs) continue;
    ranked_indices.push_back(i);
  }
  if (ranked_indices.empty()) {
    for (long long i = 0; i < static_cast<long long>(coarse_samples.size()); ++i) {
      const CalcResult& cur = coarse_samples[static_cast<std::size_t>(i)].result;
      if (IsValidPoint(cur)) ranked_indices.push_back(i);
    }
  }

  std::sort(ranked_indices.begin(), ranked_indices.end(),
            [&](long long lhs_idx, long long rhs_idx) {
              const GammaPoint& lhs = coarse_samples[static_cast<std::size_t>(lhs_idx)];
              const GammaPoint& rhs = coarse_samples[static_cast<std::size_t>(rhs_idx)];
              if (rank_feasible) {
                return BetterFeasible(lhs.result, lhs.gamma, rhs.result, rhs.gamma);
              }
              return BetterAny(lhs.result, lhs.gamma, rhs.result, rhs.gamma);
            });

  const long long remaining_budget = std::max<long long>(0, eval_budget - eval_count);
  const long long refine_rounds = (remaining_budget >= 120) ? 4 : 3;
  if (remaining_budget > 0) {
    long long candidate_count =
        std::min<long long>(6, static_cast<long long>(ranked_indices.size()));
    const long long max_candidates_by_budget =
        std::max<long long>(1, remaining_budget / (refine_rounds * 3));
    candidate_count = std::min<long long>(candidate_count, max_candidates_by_budget);
    candidate_count = std::max<long long>(1, candidate_count);

    std::vector<long long> candidate_indices;
    candidate_indices.reserve(static_cast<std::size_t>(candidate_count));
    for (long long idx : ranked_indices) {
      bool too_close = false;
      for (long long picked : candidate_indices) {
        if (std::llabs(idx - picked) <= 1) {
          too_close = true;
          break;
        }
      }
      if (too_close) continue;
      candidate_indices.push_back(idx);
      if (static_cast<long long>(candidate_indices.size()) >= candidate_count) break;
    }
    if (candidate_indices.empty()) {
      candidate_indices.push_back(ranked_indices.front());
    }

    const long long denom = std::max<long long>(
        1, static_cast<long long>(candidate_indices.size()) * refine_rounds);
    long long refine_points = remaining_budget / denom;
    if (refine_points < 3) {
      candidate_indices.clear();
    } else {
      refine_points = std::max<long long>(3, refine_points);
      refine_points = std::min<long long>(refine_points, 61);
      if ((refine_points % 2) == 0) ++refine_points;
      refine_points = std::max<long long>(3, refine_points);
    }

    for (long long center_idx : candidate_indices) {
      long double left = (center_idx > 0)
                             ? coarse_samples[static_cast<std::size_t>(center_idx - 1)].gamma
                             : g_min;
      long double right =
          (center_idx + 1 < static_cast<long long>(coarse_samples.size()))
              ? coarse_samples[static_cast<std::size_t>(center_idx + 1)].gamma
              : g_max;
      if (!(left < right)) continue;

      for (long long round = 0; round < refine_rounds; ++round) {
        if (!(left < right)) break;
        const std::vector<long double> grid = BuildLogGrid(left, right, refine_points);

        bool local_have_any = false;
        bool local_have_feasible = false;
        long long best_local_idx = -1;
        CalcResult best_local;
        long double best_local_gamma = 0.0L;

        for (long long i = 0; i < static_cast<long long>(grid.size()); ++i) {
          const long double gamma = grid[static_cast<std::size_t>(i)];
          CalcResult cur = ComputeAtGamma(args, pre, q, gamma, false, false);
          ++eval_count;
          consider(cur, gamma);
          if (!IsValidPoint(cur)) continue;

          if (cur.has_finite_l_pcs) {
            if (!local_have_feasible ||
                BetterFeasible(cur, gamma, best_local, best_local_gamma)) {
              best_local = cur;
              best_local_gamma = gamma;
              best_local_idx = i;
              local_have_feasible = true;
            }
            local_have_any = true;
            continue;
          }

          if (!local_have_feasible &&
              (!local_have_any || BetterAny(cur, gamma, best_local, best_local_gamma))) {
            best_local = cur;
            best_local_gamma = gamma;
            best_local_idx = i;
            local_have_any = true;
          }
        }

        if (!local_have_any || best_local_idx < 0) break;
        if (round + 1 >= refine_rounds) break;

        long double next_left = left;
        long double next_right = right;
        if (best_local_idx > 0) {
          next_left = grid[static_cast<std::size_t>(best_local_idx - 1)];
        }
        if (best_local_idx + 1 < static_cast<long long>(grid.size())) {
          next_right = grid[static_cast<std::size_t>(best_local_idx + 1)];
        }
        if (!(next_left < next_right) || (next_left == left && next_right == right)) {
          break;
        }
        left = next_left;
        right = next_right;
      }
    }
  }

  CalcResult out = have_feasible
                       ? ComputeAtGamma(args, pre, q, best_feasible_gamma, args.show_levels)
                       : ComputeAtGamma(args, pre, q, best_any_gamma, args.show_levels);
  out.gamma_auto = true;
  out.gamma_auto_feasible = have_feasible;
  out.gamma_search_min = g_min;
  out.gamma_search_max = g_max;
  out.gamma_search_steps = args.gamma_steps;
  out.gamma_search_evals = eval_count;
  return out;
}

void PrintResult(const CliArgs& args, const CalcResult& out) {
  std::cout << std::setprecision(18);
  std::cout << "Input parameters:\n";
  std::cout << "  d        = " << args.d << "\n";
  std::cout << "  c        = " << args.c << "\n";
  std::cout << "  k0       = " << args.k0 << "\n";
  std::cout << "  lambda   = " << args.lambda << "\n";
  std::cout << "  q        = " << out.q << "\n";
  if (out.gamma_auto) {
    std::cout << "  gamma    = " << out.gamma << " (auto-selected)\n";
    std::cout << "  gamma-search range = [" << out.gamma_search_min << ", "
              << out.gamma_search_max << "], budget=" << out.gamma_search_steps
              << ", evals=" << out.gamma_search_evals
              << "\n\n";
  } else {
    std::cout << "  gamma    = " << out.gamma
              << (args.gamma.has_value() ? " (user)"
                                         : " (default=1/(10*max(1,d)))")
              << "\n\n";
  }

  std::cout << "Code parameters:\n";
  std::cout << "  k_d      = " << out.k_d << "\n";
  std::cout << "  n_d      = " << out.n_d << "\n";
  std::cout << "  t_d      = " << out.t_d << "\n";
  std::cout << "  Delta_Cd lower bound = " << out.delta_code_lower << "\n\n";

  if (args.show_levels) {
    std::cout << "Per-level recurrence (i = 1..d):\n";
    for (const LevelInfo& lv : out.levels) {
      std::cout << "  i=" << lv.i << "  n_{i-1}=" << lv.n_prev
                << "  ell_i=" << lv.ell_i << "  t_i=" << lv.t_i << "\n";
    }
    std::cout << "\n";
  }

  std::cout << "IOPP parameters:\n";
  std::cout << "  delta_ub = J_gamma(J_gamma(Delta_Cd)) = " << out.delta_upper_johnson
            << "\n";
  std::cout << "  delta    < delta_ub (chosen near upper bound) = " << out.delta
            << "\n";
  std::cout << "  base     = 1 - delta + gamma*d        = " << out.base << "\n";
  std::cout << "  target   = 2^-lambda = " << out.target_2_minus_lambda << "\n";
  std::cout << "  budget(iopp): 2^-lambda - 2d/(gamma^3*q)          = " << out.iopp_budget
            << "\n";
  std::cout << "  budget(pcs):  2^-lambda - 2d/(gamma^3*q) - 2d/q   = " << out.pcs_budget
            << "\n";

  if (out.has_finite_l_iopp) {
    std::cout << "  l_min_iopp_only = " << out.l_min_iopp
              << "  (base^l=" << out.second_term_iopp << ")\n";
  } else {
    std::cout << "  l_min_iopp_only = N/A (no finite l)\n";
  }

  if (out.has_finite_l_pcs) {
    std::cout << "  l_min_for_PCS   = " << out.l_min_pcs
              << "  (recommended; also ensures IOPP)\n";
    std::cout << "  base^l(PCS-safe) = " << out.second_term_pcs << "\n\n";
    std::cout << "Bounds at l = l_min_for_PCS (ignoring negl(lambda)):\n";
    std::cout << "  IOPP <= 2d/(gamma^3*q) + base^l\n";
    std::cout << "       = " << out.iopp_first_term << " + " << out.second_term_pcs
              << " = " << out.iopp_bound_at_l_pcs << "\n";
    std::cout << "  PCS  <= 2d/q + 2d/(gamma^3*q) + base^l\n";
    std::cout << "       = " << out.sumcheck_term << " + " << out.iopp_first_term
              << " + " << out.second_term_pcs << " = " << out.pcs_bound_at_l_pcs
              << "\n";
  } else {
    std::cout << "  l_min_for_PCS   = N/A (no finite l)\n\n";
    std::cout << "Reason: budget(pcs) <= 0, so no l can force PCS bound below 2^-lambda.\n";
    if (out.gamma_auto && !out.gamma_auto_feasible) {
      std::cout << "Auto-search result: no feasible gamma found in the search range.\n";
    }
    std::cout << "Adjust gamma / c / q (or extension-challenge domain) to increase budget(pcs).\n";
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    const CliArgs args = ParseArgs(argc, argv);
    if (args.help) {
      PrintUsage(argv[0]);
      return 0;
    }
    const CalcResult out = Compute(args);
    PrintResult(args, out);
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << "\n\n";
    PrintUsage(argv[0]);
    return 1;
  }
}
