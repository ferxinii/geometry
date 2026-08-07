#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_b.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;



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
    // using stage_b     = grp::stage_b<expr, double>;
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



// ---------------------------------------------------------------------------
// insphere
// | ax-ex  ay-ey  az-ez  (ax-ex)^2+(ay-ey)^2+(az-ez)^2 |
// | bx-ex  by-ey  bz-ez  (bx-ex)^2+(by-ey)^2+(bz-ez)^2 |
// | cx-ex  cy-ey  cz-ez  (cx-ex)^2+(cy-ey)^2+(cz-ez)^2 |
// | dx-ex  dy-ey  dz-ez  (dx-ex)^2+(dy-ey)^2+(dz-ez)^2 |
// ---------------------------------------------------------------------------
namespace powertest_n3_k4_unweighted_impl {
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
    constexpr auto ex = grp::_13;
    constexpr auto ey = grp::_14;
    constexpr auto ez = grp::_15;

    constexpr auto dax = ax - ex;
    constexpr auto day = ay - ey;
    constexpr auto daz = az - ez;
    constexpr auto dbx = bx - ex;
    constexpr auto dby = by - ey;
    constexpr auto dbz = bz - ez;
    constexpr auto dcx = cx - ex;
    constexpr auto dcy = cy - ey;
    constexpr auto dcz = cz - ez;
    constexpr auto ddx = dx - ex;
    constexpr auto ddy = dy - ey;
    constexpr auto ddz = dz - ez;
    
    constexpr auto la = dax*dax + day*day + daz*daz;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz;
    constexpr auto lc = dcx*dcx + dcy*dcy + dcz*dcz;
    constexpr auto ld = ddx*ddx + ddy*ddy + ddz*ddz;

    // Use built-in det<> to avoid writing expansion manually
    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(daz), decltype(la),
        decltype(dbx), decltype(dby), decltype(dbz), decltype(lb),
        decltype(dcx), decltype(dcy), decltype(dcz), decltype(lc),
        decltype(ddx), decltype(ddy), decltype(ddz), decltype(ld)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n3_k4_unweighted(double ax, double ay, double az,
                                          double bx, double by, double bz,
                                          double cx, double cy, double cz,
                                          double dx, double dy, double dz,
                                          double ex, double ey, double ez)
{
    int orient_sign = orient3d_impl::pred{}.apply(ax, ay, az,
                                                  bx, by, bz,
                                                  cx, cy, cz,
                                                  dx, dy, dz);
    if (orient_sign == 0) return 0;

    int D_sign = powertest_n3_k4_unweighted_impl::pred{}.apply(ax, ay, az,
                                                               bx, by, bz,
                                                               cx, cy, cz,
                                                               dx, dy, dz,
                                                               ex, ey, ez);
    int out = -orient_sign * D_sign;
    if (out > 0) return 1;
    else if (out < 0) return -1;
    else return 0;
}


#ifdef ROBUST_PREDICATES_PRINT_SIZE
template <std::size_t N>
struct [[deprecated("results_size — see template argument")]] show_stage_d_size {};
using _size_orient3d = show_stage_d_size<orient3d_impl::exact::results_size>*;
using _size_powertest_n3_k4_unweighted = show_stage_d_size<powertest_n3_k4_unweighted_impl::exact::results_size>*;
#endif

