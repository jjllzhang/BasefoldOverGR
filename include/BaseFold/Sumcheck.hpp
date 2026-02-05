#ifndef BASEFOLD_SUMCHECK_HPP_
#define BASEFOLD_SUMCHECK_HPP_

#include <NTL/ZZ_pE.h>
#include <NTL/vec_ZZ_pE.h>

#include <vector>

namespace basefold {

using FieldElement = NTL::ZZ_pE;
using FieldVec = NTL::vec_ZZ_pE;

// A degree-2 univariate polynomial a0 + a1*X + a2*X^2.
struct QuadraticPoly {
  FieldElement a0;
  FieldElement a1;
  FieldElement a2;

  FieldElement Eval(const FieldElement &x) const;
};

// Computes the sumcheck messages for g(X)=f(X)*eq_z(X), where f is multilinear
// (given in monomial-basis coefficients).
//
// This prover never divides; it works over both finite fields and Galois rings.
class SumcheckProver {
 public:
  // Preconditions:
  // - f_coeffs.length() == 2^d for some d>=0
  // - z.size() == d
  SumcheckProver(const FieldVec &f_coeffs,
                 const std::vector<FieldElement> &z);

  long Dimension() const { return d_; }
  long RemainingVars() const { return cur_k_; }

  // Returns h_{cur_k_}(X), where cur_k_ counts the remaining variables.
  QuadraticPoly CurrentPolynomial() const;

  // Applies the verifier challenge r_{cur_k_-1} for the current variable X_{cur_k_},
  // updating internal state to the next round (cur_k_ := cur_k_-1).
  void ReceiveChallenge(const FieldElement &r_kminus1);

 private:
  long d_ = 0;
  std::vector<FieldElement> z_;

  long cur_k_ = 0;
  FieldVec coeffs_rem_;

  // Product of eq_z factors for variables already fixed to challenges
  // (i.e., variables in {cur_k_+1, ..., d_}).
  FieldElement suffix_eq_prod_;
};

// Verifier-side consistency checks for the sumcheck part of BaseFoldPCS.
//
// Conventions:
// - h_by_level[i] stores h_{i+1}(X) (so h_by_level.back() is h_d).
// - r[i] stores r_i, the challenge for variable X_{i+1}.
bool CheckSumcheckRelations(const std::vector<QuadraticPoly> &h_by_level,
                            const std::vector<FieldElement> &r,
                            const FieldElement &claimed_y);

}  // namespace basefold

#endif  // BASEFOLD_SUMCHECK_HPP_

