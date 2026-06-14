# Robust Geometric Predicates

This directory contains an implementation of robust geometric predicates
using the generic floating-point filter framework by Bartels, Fisikopoulos 
and Weiser. This MD file was edited with help from Claude.

## TODO / REVIEW
Currently, if stage d is reached for orthow_n3_k4 the required memory is huge 
(around 500 mb), so the predicate would be very slow and inefficient. Perhaps
make this optional? Or print a warning every time it is called? There might be
a better way to handle things.


## Third-party dependencies
The `includes/` directory contains header-only libraries copied from the 
repository of:

    Bartels, Tinko, Vissarion Fisikopoulos, and Martin Weiser. "Fast 
    floating-point filters for robust predicates." BIT Numerical Mathematics
    63.2 (2023): 31.

These files are provided under the Boost Software License 1.0. 
See `includes/LICENSE.txt` for details.

### Local patches

`forward_error_bound.hpp` :
    Replaced `std::abs` with a local `constexpr_abs_local` helper in the
    `abs_constant` structs inside `inexact_leaves` and `exact_leaves`.
    `std::abs` is not marked `constexpr` in AppleClang 17, causing compilation
    failures when using compile-time constants in predicate expressions.

`expression_tree.hpp` :
    Added `constexpr argument<21> _21` to support predicates with 21 input
    arguments (e.g. powertest_n3_k4 which takes 4 weighted 3D points plus
    an alpha parameter).

`expansion_arithmetic.hpp`:
    `scale_expansion` was patched to skip the non-overlapping assertion
    when `b == 0`, since the result is trivially all-zeros regardless of
    the input expansion's internal structure. The assert is a false precondition
    for this degenerate case.

    `expansion_plus` (scalar + scalar overload) had a bug where
    `two_sum_tail(e, f, *(h_begin + 1))` read uninitialized memory. Fixed to
    pass the already-computed sum `x` instead: `two_sum_tail(e, f, x)`. This
    caused overlapping expansions (assertion failure in debug, wrong results in
    release) whenever two bare argument values were added directly.

    Added `scale_expansion_pow2`: a variant of `scale_expansion` for
    power-of-2 constants. Multiplying a floating-point expansion by a power
    of 2 is exact (`two_product_tail == 0` always), so the result has the same
    number of components as the input rather than twice as many.

`expansion_eval.hpp`:
    All leaf specializations of `eval_expansion_impl` were incorrectly
    extracting runtime values from leaves using `input[right::argn - 1]`,
    which is undefined behaviour when `argn == 0` (compile-time constants).
    Fixed to use the pre-existing `get_arg_or_const<Expression>::apply(input)`
    helper, which handles `argn == 0` by returning `Expression::value` instead.

    Added `has_zero_argn<T>` void_t trait and `is_pow2_const_v<Expr>()` helper
    to detect power-of-2 compile-time constants at compile time. A simple
    `if constexpr (Expr::argn == 0)` is ill-formed on Apple Clang even with
    short-circuit `&&`, because `Expr::argn` must be syntactically valid for
    all `Expr`; the void_t trait solves this.

    `expansion_size_impl` for `operator_types::product` extended to handle
    power-of-2 constant operands: `product<E, pow2_const>` maps to `left_size`
    and `product<pow2_const, E>` maps to `right_size`, rather than the general
    `2 * left_size * right_size`. The three `eval_expansion_impl` product
    specializations route to `scale_expansion_pow2` when the constant operand
    is detected as a power of 2.

`stage_d.hpp`:
    Type aliases (`stack`, `evals`, `sizes`, `accumulated_sizes`) and
    `results_size` promoted from inside `apply()` to class scope, making
    the buffer size publicly inspectable at compile time as
    `some_pred_impl::exact::results_size`.

    Fixed stack overflow for large predicates: replaced `std::array<Real,
    results_size>` with a hybrid stack/heap strategy. Predicates with
    `results_size <= 4096` doubles use a stack-allocated `std::array`.
    Larger predicates use a `static thread_local std::vector<Real>` that
    is allocated once per thread on first use and reused for all subsequent
    calls, avoiding repeated malloc/free. Uses `results.data()` to obtain
    a raw `Real*` pointer, avoiding the `std::__wrap_iter` vs `double*`
    iterator type mismatch that occurs with `std::vector::begin()`.

`stage_b.hpp`:
    Type aliases (`stack`, `evals`, `sizes`, `accumulated_sizes`) and
    `results_size` promoted from inside `apply()` to class scope, making
    the buffer size publicly inspectable at compile time. Applied the same
    hybrid stack/heap strategy as `stage_d.hpp` for large predicates.
    Uses `Real*` throughout to avoid iterator type mismatches.


### Bugs found

**Bare leaf node crash in `stage_d` on Apple Clang / ARM64**

`stage_d` crashes with `EXC_BAD_ACCESS` inside `eval_expansion_impl` when
any argument placeholder (e.g. `grp::_5`) appears as a **bare leaf** in the
expression tree — that is, as the direct child of the root node rather than
wrapped inside a `difference<argument<i>, argument<j>>` pair.

All standard predicates (orient2d, orient3d, incircle, insphere) are unaffected
because every argument in those expressions is wrapped in a leaf-difference.
The bug only surfaces for custom predicates that add a runtime threshold or
weight directly, such as `dx*dx + wa - wb - alpha` where `alpha` is passed
as a bare `argument<N>`.

*Root cause.* When `eval_expansion_impl` is instantiated with
`LeftLeaf=true, RightLeaf=true, MostSigOnly=true`, it constructs
`std::array<Real, sizeof...(Reals)> input {{ static_cast<Real>(args)... }}`
from the forwarded variadic `const Reals&... args` pack. On Apple Clang 17 /
ARM64, these references are null by the time this line executes. Whether this
is a compiler bug or undefined behaviour in the library's variadic forwarding
chain has not been conclusively determined.

*Workaround (used in this codebase).* Express every runtime scalar as a
difference of two arguments, with the second argument always passed as `0.0`:

```cpp
// Instead of: constexpr auto expr = dx*dx + wa - wb - alpha;
//                     called as: pred{}.apply(xa, wa, xb, wb, alpha);
constexpr auto alpha_diff = grp::_5 - grp::_6;   // alpha - 0.0 = alpha exactly
constexpr auto expr = dx*dx + wa - wb - alpha_diff;
//                     called as: pred{}.apply(xa, wa, xb, wb, alpha, 0.0);
```

`two_difference_tail(alpha, 0.0, alpha) == 0` exactly, so `alpha - 0.0`
produces the two-component expansion `[0, alpha]` with zero rounding error.
This forces the root node into the `LeftLeaf=false, RightLeaf=false`
specialisation, which reads directly from the pre-computed expansion buffer
and avoids the problematic code path entirely. No library files are modified.


## Power-of-2 constant optimization

Multiplying a floating-point expansion by a power-of-2 constant (1, 2, 4,
8, ...) is exact: `two_product_tail(x, 2^k, x * 2^k) == 0` always. The
result expansion therefore has the same number of components as the input,
not twice as many.

This matters when an expression naturally involves multiplication by small
integer constants, for example `a1 * 4` where `a1` is already a two-component
expansion. Without the optimization the expansion size doubles at each such
multiplication; with it the size is unchanged.

To exploit this, define compile-time constants using `static_constant_interface`:

```cpp
template <int N>
struct int_const : grp::static_constant_interface<double> {
    static constexpr double value = static_cast<double>(N);
    static constexpr bool non_negative = (N >= 0);
};

constexpr auto a1_4 = a1 * int_const<4>{};  // size unchanged vs a1
```


## Inspecting predicate buffer sizes at compile time

After all predicate namespaces in a `.cpp` file, add:

```cpp
template <std::size_t N>
struct [[deprecated("results_size — see template argument")]] show_stage_d_size {};

using _size_pred1_d = show_stage_d_size<pred1_impl::exact::results_size>*;
using _size_pred1_b = show_stage_d_size<pred1_impl::stageb::results_size>*;
```

Each emits a compiler warning such as:
`warning: 'show_stage_d_size<953>' is deprecated`

Convert to MB: `doubles × 8 / 1,048,576`.


| Predicate        | stage_d (size) | stage_b (size) | Notes                    |
|------------------|----------------|----------------|--------------------------|
| orthow_nD_k1     | —              | —              | trivial                  |
| powertest_n1_k1  | —              | —              | fast                     |
| powertest_n1_k2  | —              | —              | fast                     |
| orthow_n1_k2     | —              | —              | fast                     |
| powertest_n2_k2  | 25.2 KB        | 2.7 KB         | manageable               |
| powertest_n2_k3  | 23.9 KB        | 3.0 KB         | manageable               |
| orthow_n2_k3     | 602 KB         | 26.7 KB        | stage_b much smaller     |
| powertest_n3_k2  | 263 KB         | 14.3 KB        | stage_b much smaller     |
| powertest_n3_k3  | 5.1 MB         | 132 KB         | slow if exact triggered  |
| powertest_n3_k4  | 657 KB         | 38.6 KB        | stage_b practical        |
| orthow_n3_k2     | —              | —              | fast                     |
| orthow_n3_k3     | 1.8 MB         | 72.3 KB        | slow if exact triggered  |
| orthow_n3_k4     | 464 MB         | 4.7 MB         | stage_d impractical      |


## Printing Error Bounds

Compile with `cmake -DBUILD_PRINT_EXPRESSIONS=ON` to build executable to
print error bounds.


## Summary of header files

  - *expression_tree.hpp*: the foundation. Defines argument<N>, sum,
    difference, product, abs, the _1.._21 placeholders, and the overloaded
    +, -, * operators that build expression trees. Everything else depends
    on this.
  - *expressions.hpp*: pre-built expression trees for standard predicates:
    orient2d, orient3d, incircle, insphere and some variants. Uses the
    generic det<> template to construct them.
  - *expression_eval.hpp*: evaluates an expression tree at runtime using
    floating point arithmetic. Takes an expression type and an array of
    input values, returns a double. Used by the filter to compute the
    approximate determinant.
  - *expansion_arithmetic.hpp*: implements Shewchuk's expansion arithmetic
    primitives: two_sum_tail, two_product_tail, grow_expansion,
    expansion_sum, scale_expansion, scale_expansion_pow2 etc. The building
    blocks for exact arithmetic.
  - *expansion_eval.hpp*: evaluates an expression tree using expansion
    arithmetic instead of floating point. Uses expansion_arithmetic.hpp to
    compute exact results. Used by stage_d. Includes power-of-2 constant
    detection and routing to scale_expansion_pow2.
  - *stage_a_error_bound.hpp*: computes Shewchuk's original stage A error
    bound for an expression tree at compile time. Produces a coefficient c
    such that the error is bounded by c·ε·M.
  - *forward_error_bound.hpp*: computes the improved error bound from the
    paper using robust_rules. More refined than stage A — produces tighter
    bounds by tracking error propagation more carefully. Defines
    forward_error_semi_static and forward_error_bound_expression. This is
    the paper's main contribution.
  - *semi_static_filter.hpp*: implements the semi-static filter. Takes an
    expression and an error bound expression, evaluates both at runtime,
    and returns +1, -1, 0, or sign_uncertain. Used by both stage A and the
    forward error bound filter.
  - *static_filter.hpp*: like the semi-static filter but the error bound is
    computed once at construction time from global input bounds, not per
    call. Faster but requires knowing input magnitude bounds in advance.
  - *almost_static_filter.hpp*: a stateful filter that starts with a static
    filter and updates it as new inputs arrive, tightening the bound over
    time.
  - *stage_a.hpp*: convenience aliases that combine stage_a_error_bound with
    semi_static_filter, static_filter, or almost_static_filter. Gives you
    stage_a_semi_static<expr>, stage_a_static<expr>,
    stage_a_almost_static<expr>.
  - *stage_b.hpp*: implements stage B: checks if all leaf differences are
    exactly representable, then runs simplified expansion arithmetic if so.
    Falls back to sign_uncertain if any leaf difference is not exact.
    Buffer size inspectable at compile time via `results_size`.
  - *stage_d.hpp*: implements stage D: full expansion arithmetic that always
    produces the exact sign. The guaranteed fallback. Never returns
    sign_uncertain. Buffer size inspectable at compile time via
    `results_size`. Uses thread-local heap storage for large predicates to
    avoid stack overflow.
  - *staged_predicate.hpp*: chains multiple stages together. Tries each in
    order, returns the first definitive result. The apply method you call
    on pred{}.
  - *simple_orient2d.hpp*: a hand-crafted semi-static and static filter
    specifically for the 2D orientation predicate, based on Ozaki et al.
    Slightly faster than the generic approach for this specific predicate.
    Also defines the phi and theta constants used in the error bound.
  - *fpg_error_bound.hpp*: implements an alternative error bound derivation
    inspired by the FPG (Floating Point Generator) approach by Meyer and
    Pion. Groups translation-invariant input variables together for a
    tighter bound. An alternative to forward_error_bound.
  - *interval_error_bound.hpp*: transforms a semi-static error bound
    expression into a static one using interval arithmetic. Takes upper and
    lower bounds on inputs and produces a single conservative error bound
    constant.
  - *predicate_approximation.hpp*: evaluates just the floating point
    approximation of an expression, returning its sign. The cheapest
    possible evaluation — used as the first step inside the filter before
    comparing against the error bound.
  - *signs_only_filter.hpp*: a filter that tries to determine the sign
    purely from the signs of subexpressions, without any error bound
    computation. Very cheap but only succeeds in simple cases.
  - *zero_pattern.hpp*: a filter that checks whether the result is provably
    zero based on the structure of the inputs. Cheap first check before
    doing any arithmetic.
  - *scaled_expansion_arithmetic.hpp*: an experimental extension of
    expansion arithmetic that handles very large magnitude differences
    between inputs. Work in progress.

