#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;

namespace orient2d_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto bx = grp::_3;
    constexpr auto by = grp::_4;
    constexpr auto cx = grp::_5;
    constexpr auto cy = grp::_6;

    constexpr auto dax  = ax - cx;
    constexpr auto day  = ay - cy;
    constexpr auto dbx  = bx - cx;
    constexpr auto dby  = by - cy;

    constexpr auto expr = dax*dby - dbx*day;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}


// ---------------------------------------------------------------------------
// incircle
// | ax-dx  ay-dy  (ax-dx)^2+(ay-dy)^2 |
// | bx-dx  by-dy  (bx-dx)^2+(by-dy)^2 |
// | cx-dx  cy-dy  (cx-dx)^2+(cy-dy)^2 |
// ---------------------------------------------------------------------------
namespace powertest_n2_k3_unweighted_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto bx = grp::_3;
    constexpr auto by = grp::_4;
    constexpr auto cx = grp::_5;
    constexpr auto cy = grp::_6;
    constexpr auto dx = grp::_7;
    constexpr auto dy = grp::_8;

    constexpr auto dax  = ax - dx;
    constexpr auto day  = ay - dy;
    constexpr auto dbx  = bx - dx;
    constexpr auto dby  = by - dy;
    constexpr auto dcx  = cx - dx;
    constexpr auto dcy  = cy - dy;

    constexpr auto la = dax*dax + day*day;
    constexpr auto lb = dbx*dbx + dby*dby;
    constexpr auto lc = dcx*dcx + dcy*dcy;
    
    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(la),
        decltype(dbx), decltype(dby), decltype(lb),
        decltype(dcx), decltype(dcy), decltype(lc)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n2_k3_unweighted(double ax, double ay,
                                          double bx, double by,
                                          double cx, double cy,
                                          double dx, double dy)
{
    int orient_sign = orient2d_impl::pred{}.apply(ax, ay,
                                                  bx, by,
                                                  cx, cy);
    if (orient_sign == 0) return 0;

    int D_sign = powertest_n2_k3_unweighted_impl::pred{}.apply(ax, ay,
                                                               bx, by,
                                                               cx, cy, 
                                                               dx, dy);
    int out = -orient_sign * D_sign;
    if (out > 0) return 1;
    else if (out < 0) return -1;
    else return 0;
}



// ---------------------------------------------------------------------------
// powertest_n2_k3
// | ax-dx  ay-dy  (ax-dx)^2+(ay-dy)^2-(wa-wd) |
// | bx-dx  by-dy  (bx-dx)^2+(by-dy)^2-(wb-wd) |
// | cx-dx  cy-dy  (cx-dx)^2+(cy-dy)^2-(wc-wd) |
// ---------------------------------------------------------------------------
namespace powertest_n2_k3_D_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto wa = grp::_3;
    constexpr auto bx = grp::_4;
    constexpr auto by = grp::_5;
    constexpr auto wb = grp::_6;
    constexpr auto cx = grp::_7;
    constexpr auto cy = grp::_8;
    constexpr auto wc = grp::_9;
    constexpr auto dx = grp::_10;
    constexpr auto dy = grp::_11;
    constexpr auto wd = grp::_12;

    constexpr auto dax = ax - dx;
    constexpr auto day = ay - dy;
    constexpr auto dwa = wa - wd;
    constexpr auto dbx = bx - dx;
    constexpr auto dby = by - dy;
    constexpr auto dwb = wb - wd;
    constexpr auto dcx = cx - dx;
    constexpr auto dcy = cy - dy;
    constexpr auto dwc = wc - wd;

    constexpr auto la = dax*dax + day*day - dwa;
    constexpr auto lb = dbx*dbx + dby*dby - dwb;
    constexpr auto lc = dcx*dcx + dcy*dcy - dwc;

    // Use built-in det<> to avoid writing expansion manually
    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(la),
        decltype(dbx), decltype(dby), decltype(lb),
        decltype(dcx), decltype(dcy), decltype(lc)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n2_k3(double ax, double ay, double wa,
                               double bx, double by, double wb,
                               double cx, double cy, double wc,
                               double dx, double dy, double wd)
{
    int orient_sign = orient2d_impl::pred{}.apply(ax, ay, bx, by, cx, cy);
    if (orient_sign == 0) return 0;

    int D_sign      = powertest_n2_k3_D_impl::pred{}.apply(ax, ay, wa,
                                                           bx, by, wb,
                                                           cx, cy, wc,
                                                           dx, dy, wd);
    int out = - orient_sign * D_sign;
    if (out > 0) return 1;
    else if (out < 0) return -1;
    else return 0;
}

// ---------------------------------------------------------------------------
// powertest_n2_k1
// sign((ax-bx)^2 + (ay-by)^2 + wa - wb)
// ---------------------------------------------------------------------------
namespace powertest_n2_k1_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto wa = grp::_3;
    constexpr auto bx = grp::_4;
    constexpr auto by = grp::_5;
    constexpr auto wb = grp::_6;

    constexpr auto dx = bx - ax;
    constexpr auto dy = by - ay;
    constexpr auto dw = wa - wb;
    constexpr auto expr = dx*dx + dy*dy + dw;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int powertest_n2_k1(double ax, double ay, double wa,
                                double bx, double by, double wb)
{
    return powertest_n2_k1_impl::pred{}.apply(ax, ay, wa, bx, by, wb);
}



// ---------------------------------------------------------------------------
// powertest_n2_k2
// Orthogonal circle power predicate for 2 weighted points in R^2.
// n-k = 0 (even), (-1)^0 = +1, so sign(pi) = sign(D).
//
// | ax-cx   ay-cy   (ax-cx)^2+(ay-cy)^2-(wa-wc)         |
// | bx-cx   by-cy   (bx-cx)^2+(by-cy)^2-(wb-wc)         |
// | by-ay   ax-bx   2*((ax-cx)(by-cy) - (ay-cy)(bx-cx)) |
// ---------------------------------------------------------------------------
namespace powertest_n2_k2_D_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto wa = grp::_3;
    constexpr auto bx = grp::_4;
    constexpr auto by = grp::_5;
    constexpr auto wb = grp::_6;
    constexpr auto cx = grp::_7;
    constexpr auto cy = grp::_8;
    constexpr auto wc = grp::_9;

    constexpr auto dax = ax - cx;
    constexpr auto day = ay - cy;
    constexpr auto dwa = wa - wc;
    constexpr auto dbx = bx - cx;
    constexpr auto dby = by - cy;
    constexpr auto dwb = wb - wc;

    // Lifted coordinates
    constexpr auto la = dax*dax + day*day - dwa;
    constexpr auto lb = dbx*dbx + dby*dby - dwb;

    // v1 = cross(b - a) = (by-ay, ax-bx) 
    constexpr auto v1x = by - ay;
    constexpr auto v1y = ax - bx;

    // mu1 = 2 * ((ax-cx)(by-cy) - (ay-cy)(bx-cx))
    // Use (expr)+(expr) instead of int_const<2>*expr to avoid zero-detection issues
    constexpr auto mu1_half = dax*dby - day*dbx;
    constexpr auto mu1 = mu1_half + mu1_half;

    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(la),
        decltype(dbx), decltype(dby), decltype(lb),
        decltype(v1x), decltype(v1y), decltype(mu1)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n2_k2(double ax, double ay, double wa,
                                double bx, double by, double wb,
                                double cx, double cy, double wc)
{
    return powertest_n2_k2_D_impl::pred{}.apply(ax, ay, wa,
                                                bx, by, wb,
                                                cx, cy, wc);
}






// ---------------------------------------------------------------------------
// powertest_n2_k1_alpha
// ---------------------------------------------------------------------------
namespace powertest_n2_k1_alpha_impl {
    constexpr auto ax    = grp::_1;
    constexpr auto ay    = grp::_2;
    constexpr auto wa    = grp::_3;
    constexpr auto bx    = grp::_4;
    constexpr auto by    = grp::_5;
    constexpr auto wb    = grp::_6;
    constexpr auto alph1 = grp::_7;
    constexpr auto zero  = grp::_8;  // I need to make alpha a subtraction, 
                                     // otherwise if it is a leaf node, the
                                     // library crashes 

    constexpr auto dx   = bx - ax;
    constexpr auto dy   = by - ay;
    constexpr auto dw = wa - wb;  // orthogonal has -wa
    constexpr auto alpha = alph1 - zero;
    constexpr auto expr = dx*dx + dy*dy + dw - alpha;  

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int powertest_n2_k1_alpha(double ax, double ay, double wa,
                                     double bx, double by, double wb,
                                     double alpha)
{
    return powertest_n2_k1_alpha_impl::pred{}.apply(ax, ay, wa, bx, by, wb, alpha, 0.0);
}


// ---------------------------------------------------------------------------
// powertest_n2_k2_alpha
// ---------------------------------------------------------------------------
namespace powertest_n2_k2_D_alpha_impl {
    constexpr auto ax    = grp::_1;
    constexpr auto ay    = grp::_2;
    constexpr auto wa    = grp::_3;
    constexpr auto bx    = grp::_4;
    constexpr auto by    = grp::_5;
    constexpr auto wb    = grp::_6;
    constexpr auto cx    = grp::_7;
    constexpr auto cy    = grp::_8;
    constexpr auto wc    = grp::_9;
    constexpr auto alpha = grp::_10;

    constexpr auto dax = ax - cx;
    constexpr auto day = ay - cy;
    constexpr auto dwa = wa - wc;
    constexpr auto dbx = bx - cx;
    constexpr auto dby = by - cy;
    constexpr auto dwb = wb - wc;

    constexpr auto la = dax*dax + day*day - dwa + alpha;
    constexpr auto lb = dbx*dbx + dby*dby - dwb + alpha;

    constexpr auto v1x = by - ay;
    constexpr auto v1y = ax - bx;

    constexpr auto mu1_h = dax*dby - day*dbx;
    constexpr auto mu1 = mu1_h + mu1_h;

    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(la),
        decltype(dbx), decltype(dby), decltype(lb),
        decltype(v1x), decltype(v1y), decltype(mu1)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n2_k2_alpha(double ax, double ay, double wa,
                                     double bx, double by, double wb,
                                     double cx, double cy, double wc,
                                     double alpha)
{
    return powertest_n2_k2_D_alpha_impl::pred{}.apply(ax, ay, wa,
                                                      bx, by, wb,
                                                      cx, cy, wc,
                                                      alpha);
}


// ---------------------------------------------------------------------------
// powertest_n2_k3_alpha
// ---------------------------------------------------------------------------
namespace powertest_n2_k3_D_alpha_impl {
    constexpr auto ax    = grp::_1;
    constexpr auto ay    = grp::_2;
    constexpr auto wa    = grp::_3;
    constexpr auto bx    = grp::_4;
    constexpr auto by    = grp::_5;
    constexpr auto wb    = grp::_6;
    constexpr auto cx    = grp::_7;
    constexpr auto cy    = grp::_8;
    constexpr auto wc    = grp::_9;
    constexpr auto dx    = grp::_10;
    constexpr auto dy    = grp::_11;
    constexpr auto wd    = grp::_12;
    constexpr auto alpha = grp::_13;

    constexpr auto dax = ax - dx;
    constexpr auto day = ay - dy;
    constexpr auto dwa = wa - wd;
    constexpr auto dbx = bx - dx;
    constexpr auto dby = by - dy;
    constexpr auto dwb = wb - wd;
    constexpr auto dcx = cx - dx;
    constexpr auto dcy = cy - dy;
    constexpr auto dwc = wc - wd;

    constexpr auto la = dax*dax + day*day - dwa + alpha;
    constexpr auto lb = dbx*dbx + dby*dby - dwb + alpha;
    constexpr auto lc = dcx*dcx + dcy*dcy - dwc + alpha;

    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(la),
        decltype(dbx), decltype(dby), decltype(lb),
        decltype(dcx), decltype(dcy), decltype(lc)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n2_k3_alpha(double ax, double ay, double wa,
                                     double bx, double by, double wb,
                                     double cx, double cy, double wc,
                                     double dx, double dy, double wd,
                                     double alpha)
{
    int orient_sign = orient2d_impl::pred{}.apply(ax, ay, bx, by, cx, cy);
    if (orient_sign == 0) return 0;

    int D_sign = powertest_n2_k3_D_alpha_impl::pred{}.apply(ax, ay, wa,
                                                            bx, by, wb,
                                                            cx, cy, wc,
                                                            dx, dy, wd,
                                                            alpha);
    int out = -orient_sign * D_sign;
    if (out > 0) return 1;
    else if (out < 0) return -1;
    else return 0;
}


