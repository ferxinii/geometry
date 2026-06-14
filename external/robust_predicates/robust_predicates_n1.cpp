#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_b.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;

template <int N>
struct int_const : grp::static_constant_interface<double> {
    static constexpr double value  = static_cast<double>(N);
    static constexpr bool non_negative = (N >= 0);
};


// ---------------------------------------------------------------------------
// powertest_n1_k2_D
// | xa-xc   (xa-xc)^2-(wa-wc) |
// | xb-xc   (xb-xc)^2-(wb-wc) |
// ---------------------------------------------------------------------------
namespace powertest_n1_k2_D_impl {
    constexpr auto xa = grp::_1;
    constexpr auto wa = grp::_2;
    constexpr auto xb = grp::_3;
    constexpr auto wb = grp::_4;
    constexpr auto xc = grp::_5;
    constexpr auto wc = grp::_6;

    constexpr auto da  = xa - xc;
    constexpr auto db  = xb - xc;
    constexpr auto dwa = wa - wc;
    constexpr auto dwb = wb - wc;

    constexpr auto la = da*da - dwa;
    constexpr auto lb = db*db - dwb;

    constexpr auto expr = da*lb - db*la;

    using filter   = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact    = grp::stage_d<expr, double>;
    // using stage_b  = grp::stage_b<expr, double>;
    using pred     = grp::staged_predicate<filter, exact>;
}


extern "C" int powertest_n1_k2(double xa, double wa,
                               double xb, double wb,
                               double xc, double wc)
{
    int orient_sign = xa>xb? 1 : xa<xb? -1 : 0;
    if (orient_sign == 0) return 0;

    int D_sign = powertest_n1_k2_D_impl::pred{}.apply(xa, wa,
                                                      xb, wb,
                                                      xc, wc);

    int out = - orient_sign * D_sign;
    if (out > 0) return 1;
    else if (out < 0) return -1;
    else return 0;
}


// ---------------------------------------------------------------------------
// powertest_n1_k1 
// sign(||xb - xa||^2 + wa - wb)  (orthogonal hypersphere to a has -wa)
// ---------------------------------------------------------------------------
namespace powertest_n1_k1_impl {
    constexpr auto xa = grp::_1;
    constexpr auto wa = grp::_2;
    constexpr auto xb = grp::_3;
    constexpr auto wb = grp::_4;

    constexpr auto dx = xb - xa;
    constexpr auto dw = wa - wb;
    constexpr auto expr = dx*dx + dw;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    // using stage_b  = grp::stage_b<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int powertest_n1_k1(double xa, double wa,
                                double xb, double wb)
{
    return powertest_n1_k1_impl::pred{}.apply(xa, wa, xb, wb);
}


// ---------------------------------------------------------------------------
// ---------------------------- RADIUS VARIANTS ------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// orthow_n1_k1
// ---------------------------------------------------------------------------
namespace orthow_n1_k1_impl {
    constexpr auto wa    = grp::_1;
    constexpr auto alpha = grp::_2;

    constexpr auto expr = wa - alpha;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    // using stage_b  = grp::stage_b<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int orthow_n1_k1(double xa, double wa,
                            double alpha)
{
    (void)xa;
    return orthow_n1_k1_impl::pred{}.apply(wa, alpha);
}

// ---------------------------------------------------------------------------
// orthow_n1_k2
// ---------------------------------------------------------------------------
namespace orthow_n1_k2_impl {
    constexpr auto xa    = grp::_1;
    constexpr auto wa    = grp::_2;
    constexpr auto xb    = grp::_3;
    constexpr auto wb    = grp::_4;
    constexpr auto alpha = grp::_5;

    constexpr auto dx  = xb - xa;
    constexpr auto dx2 = dx * dx;
    constexpr auto dw  = wb - wa;
    constexpr auto b   = dx2 - dw;

    constexpr auto a1  = wa + alpha;
    constexpr auto a1_4 = a1 * int_const<4>{};

    constexpr auto expr = a1_4 * dx2 - b * b;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    // using stage_b= grp::stage_b<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int orthow_n1_k2(double xa, double wa,
                            double xb, double wb,
                            double alpha)
{
    int orient_sign = xa>xb? 1 : xa<xb? -1 : 0;
    if (orient_sign == 0) return 0;

    return - orthow_n1_k2_impl::pred{}.apply(xa, wa,
                                             xb, wb,
                                             alpha);
}



#ifdef ROBUST_PREDICATES_PRINT_SIZE
template <std::size_t N>
struct [[deprecated("results_size — see template argument")]] show_stage_d_size {};
using _size_powertest_n1_k1 = show_stage_d_size<powertest_n1_k1_impl::exact::results_size>*;
using _size_powertest_n1_k2 = show_stage_d_size<powertest_n1_k2_D_impl::exact::results_size>*;
using _size_orthow_n1_k1 = show_stage_d_size<orthow_n1_k1_impl::exact::results_size>*;
using _size_orthow_n1_k2 = show_stage_d_size<orthow_n1_k2_impl::exact::results_size>*;

template <std::size_t N>
struct [[deprecated("results_size — see template argument")]] show_stage_b_size {};
using _size_powertest_n1_k1_b = show_stage_b_size<powertest_n1_k1_impl::stage_b::results_size>*;
using _size_powertest_n1_k2_b = show_stage_b_size<powertest_n1_k2_D_impl::stage_b::results_size>*;
using _size_orthow_n1_k1_b = show_stage_b_size<orthow_n1_k1_impl::stage_b::results_size>*;
using _size_orthow_n1_k2_b = show_stage_b_size<orthow_n1_k2_impl::stage_b::results_size>*;
#endif
