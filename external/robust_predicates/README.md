# Robust Geometric Predicates

This directory contains an implementation of robust geometric predicates
using the generic floating-point filter framework by Bartels, Fisikopoulos 
and Weiser.


## Third-party dependencies

The `includes/` directory contains header-only libraries copied from the 
repository of:

    Bartels, Tinko, Vissarion Fisikopoulos, and Martin Weiser. "Fast 
    floating-point filters for robust predicates." BIT Numerical Mathematics
    63.2 (2023): 31.

These files are provided under the Boost Software License 1.0. 
See `includes/LICENSE.txt` for details.


## Printing Error Bounds

Compile with `cmake -DBUILD_PRINT_EXPRESSIONS=ON` to build executable to
print error bounds.


### Summary of "Robust Geometric Predicates" (done with Claude)

  - *expression_tree.hpp*: the foundation. Defines argument<N>, sum, difference, product, abs, the _1.._20 placeholders, and the overloaded +, -, * operators that build expression trees. Everything else depends on this.
  - *expressions.hpp*: pre-built expression trees for standard predicates: orient2d, orient3d, incircle, insphere and some variants. Uses the generic det<> template to construct them.
  - *expression_eval.hpp*: evaluates an expression tree at runtime using floating point arithmetic. Takes an expression type and an array of input values, returns a double. Used by the filter to compute the approximate determinant.
  - *expansion_arithmetic.hpp*: implements Shewchuk's expansion arithmetic primitives: two_sum_tail, two_product_tail, grow_expansion, expansion_sum, scale_expansion etc. The building blocks for exact arithmetic.
  - *expansion_eval.hpp*: evaluates an expression tree using expansion arithmetic instead of floating point. Uses expansion_arithmetic.hpp to compute exact results. Used by stage D.
  - *stage_a_error_bound.hpp*: computes Shewchuk's original stage A error bound for an expression tree at compile time. Produces a coefficient c such that the error is bounded by c·\varepsilon·M.
  - *forward_error_bound.hpp*: computes the improved error bound from the paper using robust_rules. More refined than stage A — produces tighter bounds by tracking error propagation more carefully. Defines forward_error_semi_static and forward_error_bound_expression. This is the paper's main contribution.
  - *semi_static_filter.hpp*: implements the semi-static filter. Takes an expression and an error bound expression, evaluates both at runtime, and returns +1, -1, 0, or sign_uncertain. Used by both stage A and the forward error bound filter.
  - *static_filter.hpp*: like the semi-static filter but the error bound is computed once at construction time from global input bounds, not per call. Faster but requires knowing input magnitude bounds in advance.
  - *almost_static_filter.hpp*: a stateful filter that starts with a static filter and updates it as new inputs arrive, tightening the bound over time.
  -tage_a.hpp — convenience aliases that combine stage_a_error_bound with semi_static_filter, static_filter, or almost_static_filter. Gives you stage_a_semi_static<expr>, stage_a_static<expr>, stage_a_almost_static<expr>.
  - *stage_b.hpp*: implements stage B: checks if all leaf differences are exactly representable, then runs simplified expansion arithmetic if so. Falls back to sign_uncertain if any leaf difference is not exact.
  - *stage_d.hpp*: implements stage D: full expansion arithmetic that always produces the exact sign. The guaranteed fallback. Never returns sign_uncertain.
  - *staged_predicate.hpp*: chains multiple stages together. Tries each in order, returns the first definitive result. The apply method you call on pred{}.
  - *simple_orient2d.hpp*: a hand-crafted semi-static and static filter specifically for the 2D orientation predicate, based on Ozaki et al. Slightly faster than the generic approach for this specific predicate. Also defines the phi and theta constants used in the error bound.
  - *fpg_error_bound.hpp*: implements an alternative error bound derivation inspired by the FPG (Floating Point Generator) approach by Meyer and Pion. Groups translation-invariant input variables together for a tighter bound. An alternative to forward_error_bound.
  - *interval_error_bound.hpp*: transforms a semi-static error bound expression into a static one using interval arithmetic. Takes upper and lower bounds on inputs and produces a single conservative error bound constant.
  - *predicate_approximation.hpp*: evaluates just the floating point approximation of an expression, returning its sign. The cheapest possible evaluation — used as the first step inside the filter before comparing against the error bound.
  - *signs_only_filter.hpp*: a filter that tries to determine the sign purely from the signs of subexpressions, without any error bound computation. Very cheap but only succeeds in simple cases.
  - *zero_pattern.hpp*: a filter that checks whether the result is provably zero based on the structure of the inputs. Cheap first check before doing any arithmetic.
  - *scaled_expansion_arithmetic.hpp*: an experimental extension of expansion arithmetic that handles very large magnitude differences between inputs. Work in progress based on the comments in the file.
  - *stage_a_error_bound.hpp*: already covered above.

