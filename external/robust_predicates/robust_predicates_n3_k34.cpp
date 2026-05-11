#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;

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



// ---------------------------------------------------------------------------
// powertest_n3_k4
// | ax-ex  ay-ey  az-ez  (ax-ex)^2+(ay-ey)^2+(az-ez)^2-(wa-we) |
// | bx-ex  by-ey  bz-ez  (bx-ex)^2+(by-ey)^2+(bz-ez)^2-(wb-we) |
// | cx-ex  cy-ey  cz-ez  (cx-ex)^2+(cy-ey)^2+(cz-ez)^2-(wc-we) |
// | dx-ex  dy-ey  dz-ez  (dx-ex)^2+(dy-ey)^2+(dz-ez)^2-(wd-we) |
// ---------------------------------------------------------------------------
namespace powertest_n3_k4_D_impl {
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

extern "C" int powertest_n3_k4(double ax, double ay, double az, double wa,
                               double bx, double by, double bz, double wb,
                               double cx, double cy, double cz, double wc,
                               double dx, double dy, double dz, double wd,
                               double ex, double ey, double ez, double we)
{
    int orient_sign = orient3d_impl::pred{}.apply(ax, ay, az,
                                                  bx, by, bz,
                                                  cx, cy, cz,
                                                  dx, dy, dz);
    if (orient_sign == 0) return 0;

    int D_sign      = powertest_n3_k4_D_impl::pred{}.apply(ax, ay, az, wa,
                                                           bx, by, bz, wb,
                                                           cx, cy, cz, wc,
                                                           dx, dy, dz, wd,
                                                           ex, ey, ez, we);
    int out = - orient_sign * D_sign;
    if (out > 0) return 1;
    else if (out < 0) return -1;
    else return 0;
}


// ---------------------------------------------------------------------------
// powertest_n3_k3
// Orthogonal sphere power predicate for 3 weighted points in R^3.
// n-k = 0 (even), (-1)^0 = +1, so sign(pi) = sign(D).
// v1 = cross(b-a, c-a)
// mu1 = 2*det(a-d, b-d, c-d)
//
// | ax-dx   ay-dy   az-dz   |a-d|^2-(wa-wd) |
// | bx-dx   by-dy   bz-dz   |b-d|^2-(wb-wd) |
// | cx-dx   cy-dy   cz-dz   |c-d|^2-(wc-wd) |
// | v1x     v1y     v1z     mu1             |
// ---------------------------------------------------------------------------
namespace powertest_n3_k3_impl {
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

    constexpr auto dax = ax - dx;
    constexpr auto day = ay - dy;
    constexpr auto daz = az - dz;
    constexpr auto dwa = wa - wd;
    constexpr auto dbx = bx - dx;
    constexpr auto dby = by - dy;
    constexpr auto dbz = bz - dz;
    constexpr auto dwb = wb - wd;
    constexpr auto dcx = cx - dx;
    constexpr auto dcy = cy - dy;
    constexpr auto dcz = cz - dz;
    constexpr auto dwc = wc - wd;
    
    // Lifted coordinates
    constexpr auto la = dax*dax + day*day + daz*daz - dwa;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb;
    constexpr auto lc = dcx*dcx + dcy*dcy + dcz*dcz - dwc;

    // v1 = cross(b - a, c-a) 
    constexpr auto v1x = (by-ay)*(cz-az) - (bz-az)*(cy-ay);
    constexpr auto v1y = (bz-az)*(cx-ax) - (bx-ax)*(cz-az);
    constexpr auto v1z = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);

    // mu1 = 2 * det(a-d, b-d, c-d): expand manually and double by addition
    // (avoids nested grp::det inside grp::det which can cause zero-detection issues)
    constexpr auto mu1_half = dax*(dby*dcz - dbz*dcy)
                            - day*(dbx*dcz - dbz*dcx)
                            + daz*(dbx*dcy - dby*dcx);
    constexpr auto mu1 = mu1_half + mu1_half;

    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(daz), decltype(la),
        decltype(dbx), decltype(dby), decltype(dbz), decltype(lb),
        decltype(dcx), decltype(dcy), decltype(dcz), decltype(lc),
        decltype(v1x), decltype(v1y), decltype(v1z), decltype(mu1)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n3_k3(double ax, double ay, double az, double wa,
                               double bx, double by, double bz, double wb,
                               double cx, double cy, double cz, double wc,
                               double dx, double dy, double dz, double wd)
{
    return powertest_n3_k3_impl::pred{}.apply(ax, ay, az, wa,
                                              bx, by, bz, wb,
                                              cx, cy, cz, wc,
                                              dx, dy, dz, wd);
}




// ---------------------------------------------------------------------------
// powertest_n3_k3_alpha
// ---------------------------------------------------------------------------
namespace powertest_n3_k3_alpha_impl {
    constexpr auto ax    = grp::_1;
    constexpr auto ay    = grp::_2;
    constexpr auto az    = grp::_3;
    constexpr auto wa    = grp::_4;
    constexpr auto bx    = grp::_5;
    constexpr auto by    = grp::_6;
    constexpr auto bz    = grp::_7;
    constexpr auto wb    = grp::_8;
    constexpr auto cx    = grp::_9;
    constexpr auto cy    = grp::_10;
    constexpr auto cz    = grp::_11;
    constexpr auto wc    = grp::_12;
    constexpr auto dx    = grp::_13;
    constexpr auto dy    = grp::_14;
    constexpr auto dz    = grp::_15;
    constexpr auto wd    = grp::_16;
    constexpr auto alpha = grp::_17;

    constexpr auto dax = ax - dx;
    constexpr auto day = ay - dy;
    constexpr auto daz = az - dz;
    constexpr auto dwa = wa - wd;
    constexpr auto dbx = bx - dx;
    constexpr auto dby = by - dy;
    constexpr auto dbz = bz - dz;
    constexpr auto dwb = wb - wd;
    constexpr auto dcx = cx - dx;
    constexpr auto dcy = cy - dy;
    constexpr auto dcz = cz - dz;
    constexpr auto dwc = wc - wd;

    constexpr auto la = dax*dax + day*day + daz*daz - dwa + alpha;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb + alpha;
    constexpr auto lc = dcx*dcx + dcy*dcy + dcz*dcz - dwc + alpha;

    constexpr auto v1x = (by-ay)*(cz-az) - (bz-az)*(cy-ay);
    constexpr auto v1y = (bz-az)*(cx-ax) - (bx-ax)*(cz-az);
    constexpr auto v1z = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);

    using mu1_aux_t = grp::det <
        decltype(dax), decltype(day), decltype(daz),
        decltype(dbx), decltype(dby), decltype(dbz),
        decltype(dcx), decltype(dcy), decltype(dcz)
    >;
    // expand and double by addition to avoid zero-detection issues
    constexpr auto mu1_half = dax*(dby*dcz - dbz*dcy)
                            - day*(dbx*dcz - dbz*dcx)
                            + daz*(dbx*dcy - dby*dcx);
    constexpr auto mu1 = mu1_half + mu1_half;

    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(daz), decltype(la),
        decltype(dbx), decltype(dby), decltype(dbz), decltype(lb),
        decltype(dcx), decltype(dcy), decltype(dcz), decltype(lc),
        decltype(v1x), decltype(v1y), decltype(v1z), decltype(mu1)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n3_k3_alpha(double ax, double ay, double az, double wa,
                                     double bx, double by, double bz, double wb,
                                     double cx, double cy, double cz, double wc,
                                     double dx, double dy, double dz, double wd,
                                     double alpha)
{
    return powertest_n3_k3_alpha_impl::pred{}.apply(ax, ay, az, wa,
                                                    bx, by, bz, wb,
                                                    cx, cy, cz, wc,
                                                    dx, dy, dz, wd,
                                                    alpha);
}


// ---------------------------------------------------------------------------
// powertest_n3_k4_alpha
// ---------------------------------------------------------------------------
namespace powertest_n3_k4_D_alpha_impl {
    constexpr auto ax    = grp::_1;
    constexpr auto ay    = grp::_2;
    constexpr auto az    = grp::_3;
    constexpr auto wa    = grp::_4;
    constexpr auto bx    = grp::_5;
    constexpr auto by    = grp::_6;
    constexpr auto bz    = grp::_7;
    constexpr auto wb    = grp::_8;
    constexpr auto cx    = grp::_9;
    constexpr auto cy    = grp::_10;
    constexpr auto cz    = grp::_11;
    constexpr auto wc    = grp::_12;
    constexpr auto dx    = grp::_13;
    constexpr auto dy    = grp::_14;
    constexpr auto dz    = grp::_15;
    constexpr auto wd    = grp::_16;
    constexpr auto ex    = grp::_17;
    constexpr auto ey    = grp::_18;
    constexpr auto ez    = grp::_19;
    constexpr auto we    = grp::_20;
    constexpr auto alpha = grp::_21;

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

    constexpr auto la = dax*dax + day*day + daz*daz - dwa + alpha;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb + alpha;
    constexpr auto lc = dcx*dcx + dcy*dcy + dcz*dcz - dwc + alpha;
    constexpr auto ld = ddx*ddx + ddy*ddy + ddz*ddz - dwd + alpha;

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

extern "C" int powertest_n3_k4_alpha(double ax, double ay, double az, double wa,
                                     double bx, double by, double bz, double wb,
                                     double cx, double cy, double cz, double wc,
                                     double dx, double dy, double dz, double wd,
                                     double ex, double ey, double ez, double we,
                                     double alpha)
{
    int orient_sign = orient3d_impl::pred{}.apply(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz);
    if (orient_sign == 0) return 0;

    int D_sign = powertest_n3_k4_D_alpha_impl::pred{}.apply(ax, ay, az, wa,
                                                            bx, by, bz, wb,
                                                            cx, cy, cz, wc,
                                                            dx, dy, dz, wd,
                                                            ex, ey, ez, we,
                                                            alpha);
    int out = -orient_sign * D_sign;
    if (out > 0) return 1;
    else if (out < 0) return -1;
    else return 0;
}

