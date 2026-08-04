#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_b.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;


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



