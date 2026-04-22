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

// ---------------------------------------------------------------------------
// incircle
// | ax-dx  ay-dy  (ax-dx)^2+(ay-dy)^2 |
// | bx-dx  by-dy  (bx-dx)^2+(by-dy)^2 |
// | cx-dx  cy-dy  (cx-dx)^2+(cy-dy)^2 |
// ---------------------------------------------------------------------------
namespace incircle_impl {
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

extern "C" int incircle(double ax, double ay,
                        double bx, double by,
                        double cx, double cy,
                        double dx, double dy)
{
    return incircle_impl::pred{}.apply(ax, ay, bx, by, cx, cy, dx, dy);
}

// ---------------------------------------------------------------------------
// insphere
// | ax-ex  ay-ey  az-ez  (ax-ex)^2+(ay-ey)^2+(az-ez)^2 |
// | bx-ex  by-ey  bz-ez  (bx-ex)^2+(by-ey)^2+(bz-ez)^2 |
// | cx-ex  cy-ey  cz-ez  (cx-ex)^2+(cy-ey)^2+(cz-ez)^2 |
// | dx-ex  dy-ey  dz-ez  (dx-ex)^2+(dy-ey)^2+(dz-ez)^2 |
// ---------------------------------------------------------------------------
namespace insphere_impl {
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
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int insphere(double ax, double ay, double az,
                        double bx, double by, double bz,
                        double cx, double cy, double cz,
                        double dx, double dy, double dz,
                        double ex, double ey, double ez)
{
    return insphere_impl::pred{}.apply(ax, ay, az, bx, by, bz,
                                       cx, cy, cz, dx, dy, dz,
                                       ex, ey, ez);
}





// ---------------------------------------------------------------------------
// powertest1d
// | xa-xc   (xa-xc)^2-(wa-wc) |
// | xb-xc   (xb-xc)^2-(wb-wc) |
// ---------------------------------------------------------------------------
namespace powertest1d_impl {
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


extern "C" int powertest1d(double xa, double wa,
                           double xb, double wb,
                           double xc, double wc)
{
    return powertest1d_impl::pred{}.apply(xa, wb, xb, wb, xc, wc);
}



// ---------------------------------------------------------------------------
// powertest2d
// | ax-dx  ay-dy  (ax-dx)^2+(ay-dy)^2-(wa-wd) |
// | bx-dx  by-dy  (bx-dx)^2+(by-dy)^2-(wb-wd) |
// | cx-dx  cy-dy  (cx-dx)^2+(cy-dy)^2-(wc-wd) |
// ---------------------------------------------------------------------------
namespace powertest2d_impl {
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

extern "C" int powertest2d(double ax, double ay, double wa,
                           double bx, double by, double wb,
                           double cx, double cy, double wc,
                           double dx, double dy, double wd)
{
    return powertest2d_impl::pred{}.apply(ax, ay, wa,
                                          bx, by, wb,
                                          cx, cy, wc,
                                          dx, dy, wd);
}


// ---------------------------------------------------------------------------
// powertest3d
// | ax-ex  ay-ey  az-ez  (ax-ex)^2+(ay-ey)^2+(az-ez)^2-(wa-we) |
// | bx-ex  by-ey  bz-ez  (bx-ex)^2+(by-ey)^2+(bz-ez)^2-(wb-we) |
// | cx-ex  cy-ey  cz-ez  (cx-ex)^2+(cy-ey)^2+(cz-ez)^2-(wc-we) |
// | dx-ex  dy-ey  dz-ez  (dx-ex)^2+(dy-ey)^2+(dz-ez)^2-(wd-we) |
// ---------------------------------------------------------------------------
namespace powertest3d_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto az = grp::_3;
    constexpr auto wa = grp::_4;
    constexpr auto bx = grp::_5;
    constexpr auto by = grp::_6;
    constexpr auto bz = grp::_7;
    constexpr auto wb = grp::_8;
    constexpr auto cx = grp::_9;
    constexpr auto cy = grp::_10;
    constexpr auto cz = grp::_11;
    constexpr auto wc = grp::_12;
    constexpr auto dx = grp::_13;
    constexpr auto dy = grp::_14;
    constexpr auto dz = grp::_15;
    constexpr auto wd = grp::_16;
    constexpr auto ex = grp::_17;
    constexpr auto ey = grp::_18;
    constexpr auto ez = grp::_19;
    constexpr auto we = grp::_20;

    constexpr auto dax = ax - ex;
    constexpr auto day = ay - ey;
    constexpr auto daz = az - ez;
    constexpr auto dwa = wa - we;
    constexpr auto dbx = bx - ex;
    constexpr auto dby = by - ey;
    constexpr auto dbz = bz - ez;
    constexpr auto dwb = wb - we;
    constexpr auto dcx = cx - ex;
    constexpr auto dcy = cy - ey;
    constexpr auto dcz = cz - ez;
    constexpr auto dwc = wc - we;
    constexpr auto ddx = dx - ex;
    constexpr auto ddy = dy - ey;
    constexpr auto ddz = dz - ez;
    constexpr auto dwd = wd - we;

    constexpr auto la = dax*dax + day*day + daz*daz - dwa;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb;
    constexpr auto lc = dcx*dcx + dcy*dcy + dcz*dcz - dwc;
    constexpr auto ld = ddx*ddx + ddy*ddy + ddz*ddz - dwd;

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
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest3d(double ax, double ay, double az, double wa,
                           double bx, double by, double bz, double wb,
                           double cx, double cy, double cz, double wc,
                           double dx, double dy, double dz, double wd,
                           double ex, double ey, double ez, double we)
{
    return powertest3d_impl::pred{}.apply(ax, ay, az, wa,
                                          bx, by, bz, wb,
                                          cx, cy, cz, wc,
                                          dx, dy, dz, wd,
                                          ex, ey, ez, we);
}

