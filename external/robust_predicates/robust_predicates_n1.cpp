#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;


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
    using pred     = grp::staged_predicate<filter, exact>;
}

namespace powertest_n1_k2_orient_impl {
    constexpr auto xa = grp::_1;
    constexpr auto xb = grp::_2;

    constexpr auto expr  = xa - xb;

    using filter   = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact    = grp::stage_d<expr, double>;
    using pred     = grp::staged_predicate<filter, exact>;
} 


extern "C" int powertest_n1_k2(double xa, double wa,
                               double xb, double wb,
                               double xc, double wc)
{
    int orient_sign = powertest_n1_k2_orient_impl::pred{}.apply(xa, xb);
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
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int powertest_n1_k1(double xa, double wa,
                                double xb, double wb)
{
    return powertest_n1_k1_impl::pred{}.apply(xa, wa, xb, wb);
}


// ---------------------------------------------------------------------------
// ---------------------------- ALPHA VARIANTS -------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// powertest_n1_k1_alpha
// ---------------------------------------------------------------------------
namespace powertest_n1_k1_alpha_impl {
    constexpr auto xa    = grp::_1;
    constexpr auto wa    = grp::_2;
    constexpr auto xb    = grp::_3;
    constexpr auto wb    = grp::_4;
    constexpr auto alph1 = grp::_5;
    constexpr auto zero  = grp::_6;  // I need to make alpha a subtraction, 
                                     // otherwise if it is a leaf node, the
                                     // library crashes 

    constexpr auto dx   = xb - xa;
    constexpr auto dw   = wa - wb;  // orthogonal hypersphere to a has -wa
    constexpr auto alpha = alph1 - zero;
    constexpr auto expr = dx*dx + dw - alpha;  

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int powertest_n1_k1_alpha(double xa, double wa,
                                     double xb, double wb,
                                     double alpha)
{
    return powertest_n1_k1_alpha_impl::pred{}.apply(xa, wa, xb, wb, alpha, 0.0);
}

// ---------------------------------------------------------------------------
// powertest_n1_k2_alpha
// ---------------------------------------------------------------------------
namespace powertest_n1_k2_D_alpha_impl {
    constexpr auto xa    = grp::_1;
    constexpr auto wa    = grp::_2;
    constexpr auto xb    = grp::_3;
    constexpr auto wb    = grp::_4;
    constexpr auto xc    = grp::_5;
    constexpr auto wc    = grp::_6;
    constexpr auto alpha = grp::_7;

    constexpr auto da  = xa - xc;
    constexpr auto db  = xb - xc;
    constexpr auto dwa = wa - wc;
    constexpr auto dwb = wb - wc;

    constexpr auto la = da*da - dwa + alpha;
    constexpr auto lb = db*db - dwb + alpha;

    constexpr auto expr = da*lb - db*la;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int powertest_n1_k2_alpha(double xa, double wa,
                                     double xb, double wb,
                                     double xc, double wc,
                                     double alpha)
{
    int orient_sign = powertest_n1_k2_orient_impl::pred{}.apply(xa, xb);
    if (orient_sign == 0) return 0;

    int D_sign = powertest_n1_k2_D_alpha_impl::pred{}.apply(xa, wa,
                                                            xb, wb,
                                                            xc, wc,
                                                            alpha);

    int out = -orient_sign * D_sign;
    if (out > 0) return 1;
    else if (out < 0) return -1;
    else return 0;
}

