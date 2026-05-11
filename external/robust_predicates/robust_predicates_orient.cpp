#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
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
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int orient2d(double ax, double ay,
                        double bx, double by,
                        double cx, double cy)
{
    return orient2d_impl::pred{}.apply(ax, ay, bx, by, cx, cy);
}

// ---------------------------------------------------------------------------
// orient3d
// | ax-dx  ay-dy  az-dz |
// | bx-dx  by-dy  bz-dz |
// | cx-dx  cy-dy  cz-dz |
// ---------------------------------------------------------------------------
namespace orient3d_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto az = grp::_3;
    constexpr auto bx = grp::_4;
    constexpr auto by = grp::_5;
    constexpr auto bz = grp::_6;
    constexpr auto cx = grp::_7;
    constexpr auto cy = grp::_8;
    constexpr auto cz = grp::_9;
    constexpr auto dx = grp::_10;
    constexpr auto dy = grp::_11;
    constexpr auto dz = grp::_12;

    constexpr auto dax = ax - dx;
    constexpr auto day = ay - dy;
    constexpr auto daz = az - dz;
    constexpr auto dbx = bx - dx;
    constexpr auto dby = by - dy;
    constexpr auto dbz = bz - dz;
    constexpr auto dcx = cx - dx;
    constexpr auto dcy = cy - dy;
    constexpr auto dcz = cz - dz;

    // Use built-in det<> to avoid writing expansion manually
    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(daz),
        decltype(dbx), decltype(dby), decltype(dbz),
        decltype(dcx), decltype(dcy), decltype(dcz)
    >;
    constexpr auto expr = expr_t{};


    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int orient3d(double ax, double ay, double az,
                        double bx, double by, double bz,
                        double cx, double cy, double cz,
                        double dx, double dy, double dz)
{
    return orient3d_impl::pred{}.apply(ax, ay, az, bx, by, bz,
                                       cx, cy, cz, dx, dy, dz);
}
