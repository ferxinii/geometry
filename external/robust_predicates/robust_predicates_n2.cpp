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
// orient2d
// | ax-cx  ay-cy |
// | bx-cx  by-cy |
// ---------------------------------------------------------------------------
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
    // using stage_b     = grp::stage_b<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int orient2d(double ax, double ay,
                        double bx, double by,
                        double cx, double cy)
{
    return orient2d_impl::pred{}.apply(ax, ay, bx, by, cx, cy);
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
    // using stage_b     = grp::stage_b<expr, double>;
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
    // using stage_b     = grp::stage_b<expr, double>;
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
    constexpr auto mu1 = mu1_half * int_const<2>{};

    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(la),
        decltype(dbx), decltype(dby), decltype(lb),
        decltype(v1x), decltype(v1y), decltype(mu1)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
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
    // using stage_b     = grp::stage_b<expr, double>;
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
// ---------------------------- RADIUS VARIANTS ------------------------------
// ---------------------------------------------------------------------------

namespace orthow_n2_k1_impl {
    constexpr auto wa    = grp::_1;
    constexpr auto alpha = grp::_2;

    constexpr auto expr = wa - alpha;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int orthow_n2_k1(double ax, double ay, double wa,
                            double alpha)
{
    (void)ax;  (void)ay; 
    return orthow_n2_k1_impl::pred{}.apply(wa, alpha);
}


// ---------------------------------------------------------------------------
// orthow_n2_k2
// ---------------------------------------------------------------------------
namespace orthow_n2_k2_impl {
    constexpr auto ax    = grp::_1;
    constexpr auto ay    = grp::_2;
    constexpr auto wa    = grp::_3;
    constexpr auto bx    = grp::_4;
    constexpr auto by    = grp::_5;
    constexpr auto wb    = grp::_6;
    constexpr auto alpha = grp::_7;

    constexpr auto dx  = bx - ax;
    constexpr auto dx2 = dx * dx;
    constexpr auto dy  = by - ay;
    constexpr auto dy2 = dy * dy;
    constexpr auto dw  = wb - wa;

    constexpr auto G = dx2 + dy2;
    constexpr auto b = dx2 + dy2 - dw;
    constexpr auto a1  = wa + alpha;
    constexpr auto a1_4 = a1 * int_const<4>{};

    constexpr auto expr = a1_4 * G - b * b;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int orthow_n2_k2(double ax, double ay, double wa,
                            double bx, double by, double wb,
                            double alpha)
{
    int orient_sign = ((ax == bx) && (ay == by)) ? 0 : 1;
    if (orient_sign == 0) return 0;

    return - orthow_n2_k2_impl::pred{}.apply(ax, ay, wa,
                                             bx, by, wb,
                                             alpha);
}


// ---------------------------------------------------------------------------
// orthow_n2_k3
// ---------------------------------------------------------------------------
namespace orthow_n2_k3_impl {
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

    constexpr auto dbx  = bx - ax;
    constexpr auto dbx2 = dbx * dbx;
    constexpr auto dby  = by - ay;
    constexpr auto dby2 = dby * dby;
    constexpr auto dwb  = wb - wa;

    constexpr auto dcx  = cx - ax;
    constexpr auto dcx2 = dcx * dcx;
    constexpr auto dcy  = cy - ay;
    constexpr auto dcy2 = dcy * dcy;
    constexpr auto dwc  = wc - wa;

    constexpr auto a1  = wa + alpha;
    constexpr auto a1_4 = a1 * int_const<4>{};

    constexpr auto b1 = dbx2 + dby2 - dwb;
    constexpr auto b2 = dcx2 + dcy2 - dwc;
    constexpr auto G11 = dbx2 + dby2;
    constexpr auto G12 = dbx * dcx + dby * dcy;
    constexpr auto G22 = dcx2 + dcy2;

    using expr_t = grp::det <
        decltype(a1_4), decltype(b1),  decltype(b2),
        decltype(b1),   decltype(G11), decltype(G12),
        decltype(b2),   decltype(G12), decltype(G22)
    >;
    constexpr auto expr = expr_t{};

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int orthow_n2_k3(double ax, double ay, double wa,
                            double bx, double by, double wb,
                            double cx, double cy, double wc,
                            double alpha)
{
    int orient_sign = orient2d_impl::pred{}.apply(ax, ay, bx, by, cx, cy);
    if (orient_sign == 0) return 0;

    return - orthow_n2_k3_impl::pred{}.apply(ax, ay, wa,
                                             bx, by, wb,
                                             cx, cy, wc,
                                             alpha);
}


template <std::size_t N>
struct [[deprecated("results_size — see template argument")]] show_stage_d_size {};
using _size_powertest_n2_k1 = show_stage_d_size<powertest_n2_k1_impl::exact::results_size>*;
using _size_powertest_n2_k2 = show_stage_d_size<powertest_n2_k2_D_impl::exact::results_size>*;
using _size_powertest_n2_k3 = show_stage_d_size<powertest_n2_k3_D_impl::exact::results_size>*;
using _size_orthow_n2_k1 = show_stage_d_size<orthow_n2_k1_impl::exact::results_size>*;
using _size_orthow_n2_k2 = show_stage_d_size<orthow_n2_k2_impl::exact::results_size>*;
using _size_orthow_n2_k3 = show_stage_d_size<orthow_n2_k3_impl::exact::results_size>*;

// template <std::size_t N>
// struct [[deprecated("results_size — see template argument")]] show_stage_b_size {};
// using _size_powertest_n2_k1_b = show_stage_b_size<powertest_n2_k1_impl::stage_b::results_size>*;
// using _size_powertest_n2_k2_b = show_stage_b_size<powertest_n2_k2_D_impl::stage_b::results_size>*;
// using _size_powertest_n2_k3_b = show_stage_b_size<powertest_n2_k3_D_impl::stage_b::results_size>*;
// using _size_orthow_n2_k1_b = show_stage_b_size<orthow_n2_k1_impl::stage_b::results_size>*;
// using _size_orthow_n2_k2_b = show_stage_b_size<orthow_n2_k2_impl::stage_b::results_size>*;
// using _size_orthow_n2_k3_b = show_stage_b_size<orthow_n2_k3_impl::stage_b::results_size>*;
