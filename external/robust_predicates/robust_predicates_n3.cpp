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


// ---------------------------------------------------------------------------
// powertest_n3_k1
// sign((ax-bx)^2 + (ay-by)^2 + (az-bz)^2 + wa - wb)
// ---------------------------------------------------------------------------
namespace powertest_n3_k1_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto az = grp::_3;
    constexpr auto wa = grp::_4;
    constexpr auto bx = grp::_5;
    constexpr auto by = grp::_6;
    constexpr auto bz = grp::_7;
    constexpr auto wb = grp::_8;

    constexpr auto dx = bx - ax;
    constexpr auto dy = by - ay;
    constexpr auto dz = bz - az;
    constexpr auto dw = wa - wb;
    constexpr auto expr = dx*dx + dy*dy + dz*dz + dw;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int powertest_n3_k1(double ax, double ay, double az, double wa,
                                double bx, double by, double bz, double wb)
{
    return powertest_n3_k1_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb);
}


// ---------------------------------------------------------------------------
// powertest_n3_k2
// Orthogonal sphere power predicate for 2 weighted points in R^3.
// n-k = 1 (odd), (-1)^0 = +1, so sign(pi) = sign(orient(a,b,ej1,ej2)) * sign(D)
//
// Basis pair chosen by smallest absolute component of b-a (least parallel).
//
// | ax-cx   ay-cy   az-cz   |a-c|^2-(wa-wc) |
// | bx-cx   by-cy   bz-cz   |b-c|^2-(wb-wc) |
// | v1x     v1y     v1z     mu1             |
// | v2x     v2y     v2z     mu2             |
// ---------------------------------------------------------------------------
// --- ej1=ex, ej2=ey ---
namespace powertest_n3_k2_D_xy_impl {
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

    constexpr auto dax = ax - cx;
    constexpr auto day = ay - cy;
    constexpr auto daz = az - cz;
    constexpr auto dwa = wa - wc;
    constexpr auto dbx = bx - cx;
    constexpr auto dby = by - cy;
    constexpr auto dbz = bz - cz;
    constexpr auto dwb = wb - wc;

    constexpr auto la = dax*dax + day*day + daz*daz - dwa;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb;

    // v1 = cross(b-a, ej2) = (az-bz, 0, bx-ax)
    constexpr auto v1x = az - bz;
    // v1y = 0  (not needed)
    constexpr auto v1z = bx - ax;
    // v2 = cross(b-a, ej1) = (0, bz-az, ay-by)
    // v2x = 0  (not needed)
    constexpr auto v2y = bz - az;
    constexpr auto v2z = ay - by;

    // mu1 = 2*det(da, db, ey) = 2*(daz*dbx - dax*dbz)
    constexpr auto mu1_h = daz*dbx - dax*dbz;
    constexpr auto mu1 = mu1_h * int_const<2>{};
    // mu2 = 2*det(da, db, ex) = 2*(day*dbz - daz*dby)
    constexpr auto mu2_h = day*dbz - daz*dby;
    constexpr auto mu2 = mu2_h * int_const<2>{};

    constexpr auto det2_ab_xy  = dax*dby - dbx*day;
    constexpr auto det2_v1z_v2z = v1z*mu2 - v2z*mu1;
    constexpr auto det2_dbz_v1z = dbz*mu1 - v1z*lb;
    constexpr auto det2_daz_v1z = daz*mu1 - v1z*la;

    using det3_t = grp::det <
        decltype(day), decltype(daz), decltype(la),
        decltype(dby), decltype(dbz), decltype(lb),
        decltype(v2y), decltype(v2z), decltype(mu2)
    >;
    constexpr auto det3 = det3_t{};

    constexpr auto expr = det2_ab_xy * det2_v1z_v2z
                        + v2y * (dax*det2_dbz_v1z - dbx*det2_daz_v1z)
                        + v1x * det3;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace powertest_n3_k2_orient_xy_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto az = grp::_3;
    constexpr auto wa = grp::_4;  // unused but keeps interface uniform
    constexpr auto bx = grp::_5;
    constexpr auto by = grp::_6;
    constexpr auto bz = grp::_7;
    constexpr auto wb = grp::_8;  // unused but keeps interface uniform
    constexpr auto cx = grp::_9;  // unused but keeps interface uniform
    constexpr auto cy = grp::_10; // unused but keeps interface uniform
    constexpr auto cz = grp::_11; // unused but keeps interface uniform
    constexpr auto wc = grp::_12; // unused but keeps interface uniform

    // orient3(a, b, ex, ey): (ax+ay-1)*bz - (bx+by-1)*az
    constexpr auto expr = (ax + ay)*bz - (bx + by)*az + (az - bz);

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

// --- ej1=ey, ej2=ez ---
namespace powertest_n3_k2_D_yz_impl {
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

    constexpr auto dax = ax - cx;
    constexpr auto day = ay - cy;
    constexpr auto daz = az - cz;
    constexpr auto dwa = wa - wc;
    constexpr auto dbx = bx - cx;
    constexpr auto dby = by - cy;
    constexpr auto dbz = bz - cz;
    constexpr auto dwb = wb - wc;

    constexpr auto la = dax*dax + day*day + daz*daz - dwa;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb;

    // v1 = cross(b-a, ej2) = cross(b-a, ez) = (by-ay, ax-bx, 0)
    constexpr auto v1x = by - ay;
    constexpr auto v1y = ax - bx;
    // v1z = 0  (not needed)
    // v2 = cross(b-a, ej1) = cross(b-a, ey) = (az-bz, 0, bx-ax)
    constexpr auto v2x = az - bz;
    // v2y = 0  (not needed)
    constexpr auto v2z = bx - ax;

    // mu1 = 2*det(da, db, ez) = 2*(dax*dby - day*dbx)
    constexpr auto mu1_h = dax*dby - day*dbx;
    constexpr auto mu1 = mu1_h * int_const<2>{};
    // mu2 = 2*det(da, db, ey) = 2*(daz*dbx - dax*dbz)
    constexpr auto mu2_h = daz*dbx - dax*dbz;
    constexpr auto mu2 = mu2_h * int_const<2>{};

    constexpr auto det2_ab_yz  = day*dbz - dby*daz;
    constexpr auto det2_v1x_v2x = v1x*mu2 - v2x*mu1;
    constexpr auto det2_dbx_v1x = dbx*mu1 - v1x*lb;
    constexpr auto det2_dax_v1x = dax*mu1 - v1x*la;

    using det3_t = grp::det <
        decltype(daz), decltype(dax), decltype(la),
        decltype(dbz), decltype(dbx), decltype(lb),
        decltype(v2z), decltype(v2x), decltype(mu2)
    >;
    constexpr auto det3 = det3_t{};

    constexpr auto expr = det2_ab_yz * det2_v1x_v2x
                        + v2z * (day*det2_dbx_v1x - dby*det2_dax_v1x)
                        + v1y * det3;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace powertest_n3_k2_orient_yz_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto az = grp::_3;
    constexpr auto wa = grp::_4;  // unused
    constexpr auto bx = grp::_5;
    constexpr auto by = grp::_6;
    constexpr auto bz = grp::_7;
    constexpr auto wb = grp::_8;  // unused
    constexpr auto cx = grp::_9;  // unused
    constexpr auto cy = grp::_10; // unused
    constexpr auto cz = grp::_11; // unused
    constexpr auto wc = grp::_12; // unused

    // orient3(a, b, ey, ez): (ay+az-1)*bx - (by+bz-1)*ax
    constexpr auto expr = (ay + az)*bx - (by + bz)*ax + (ax - bx);

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

// --- ej1=ez, ej2=ex ---
namespace powertest_n3_k2_D_zx_impl {
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

    constexpr auto dax = ax - cx;
    constexpr auto day = ay - cy;
    constexpr auto daz = az - cz;
    constexpr auto dwa = wa - wc;
    constexpr auto dbx = bx - cx;
    constexpr auto dby = by - cy;
    constexpr auto dbz = bz - cz;
    constexpr auto dwb = wb - wc;

    constexpr auto la = dax*dax + day*day + daz*daz - dwa;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb;

    // v1 = cross(b-a, ej2) = cross(b-a, ex) = (0, bz-az, ay-by)
    // v1x = 0  (not needed)
    constexpr auto v1y = bz - az;
    constexpr auto v1z = ay - by;
    // v2 = cross(b-a, ej1) = cross(b-a, ez) = (by-ay, ax-bx, 0)
    constexpr auto v2x = by - ay;
    constexpr auto v2y = ax - bx;
    // v2z = 0  (not needed)

    // mu1 = 2*det(da, db, ex) = 2*(day*dbz - daz*dby)
    constexpr auto mu1_h = day*dbz - daz*dby;
    constexpr auto mu1 = mu1_h * int_const<2>{};
    // mu2 = 2*det(da, db, ez) = 2*(dax*dby - day*dbx)
    constexpr auto mu2_h = dax*dby - day*dbx;
    constexpr auto mu2 = mu2_h * int_const<2>{};

    constexpr auto det2_ab_zx  = daz*dbx - dbz*dax;
    constexpr auto det2_v1y_v2y = v1y*mu2 - v2y*mu1;
    constexpr auto det2_dby_v1y = dby*mu1 - v1y*lb;
    constexpr auto det2_day_v1y = day*mu1 - v1y*la;

    using det3_t = grp::det <
        decltype(dax), decltype(day), decltype(la),
        decltype(dbx), decltype(dby), decltype(lb),
        decltype(v2x), decltype(v2y), decltype(mu2)
    >;
    constexpr auto det3 = det3_t{};

    constexpr auto expr = det2_ab_zx * det2_v1y_v2y
                        + v2x * (daz*det2_dby_v1y - dbz*det2_day_v1y)
                        + v1z * det3;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace powertest_n3_k2_orient_zx_impl {
    constexpr auto ax = grp::_1;
    constexpr auto ay = grp::_2;
    constexpr auto az = grp::_3;
    constexpr auto wa = grp::_4;  // unused
    constexpr auto bx = grp::_5;
    constexpr auto by = grp::_6;
    constexpr auto bz = grp::_7;
    constexpr auto wb = grp::_8;  // unused
    constexpr auto cx = grp::_9;  // unused
    constexpr auto cy = grp::_10; // unused
    constexpr auto cz = grp::_11; // unused
    constexpr auto wc = grp::_12; // unused

    // orient3(a, b, ez, ex) = (ax+az-1)*by - (bx+bz-1)*ay
    constexpr auto expr = (ax + az)*by - (bx + bz)*ay + (ay - by);

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n3_k2(double ax, double ay, double az, double wa,
                               double bx, double by, double bz, double wb,
                               double cx, double cy, double cz, double wc)
{
    double abx = std::abs(bx - ax);
    double aby = std::abs(by - ay);
    double abz = std::abs(bz - az);

    int D_sign, orient_sign;

    if (abz >= abx && abz >= aby) {
        // z is LARGEST -> ej1=ex, ej2=ey (exclude ez which would be parallel)
        orient_sign = powertest_n3_k2_orient_xy_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc);
        if (orient_sign == 0) return 0;
        D_sign      = powertest_n3_k2_D_xy_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc);
    } else if (abx >= aby && abx >= abz) {
        // x is LARGEST -> ej1=ey, ej2=ez (exclude ex which would be parallel)
        orient_sign = powertest_n3_k2_orient_yz_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc);
        if (orient_sign == 0) return 0;
        D_sign      = powertest_n3_k2_D_yz_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc);
    } else {
        // y is smallest -> ej1=ez, ej2=ex
        orient_sign = powertest_n3_k2_orient_zx_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc);
        if (orient_sign == 0) return 0;
        D_sign      = powertest_n3_k2_D_zx_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc);
    }

    int out = orient_sign * D_sign;
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
    constexpr auto mu1 = mu1_half * int_const<2>{};

    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(daz), decltype(la),
        decltype(dbx), decltype(dby), decltype(dbz), decltype(lb),
        decltype(dcx), decltype(dcy), decltype(dcz), decltype(lc),
        decltype(v1x), decltype(v1y), decltype(v1z), decltype(mu1)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
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
    // using stage_b     = grp::stage_b<expr, double>;
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
// ---------------------------- RADIUS VARIANTS ------------------------------
// ---------------------------------------------------------------------------

namespace orthow_n3_k1_impl {
    constexpr auto wa    = grp::_1;
    constexpr auto alpha = grp::_2;

    constexpr auto expr = wa - alpha;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int orthow_n3_k1(double ax, double ay, double az, double wa,
                            double alpha)
{
    (void)ax;  (void)ay; 
    return orthow_n3_k1_impl::pred{}.apply(wa, alpha);
}



// ---------------------------------------------------------------------------
// orthow_n2_k2
// ---------------------------------------------------------------------------
namespace orthow_n3_k2_impl {
    constexpr auto ax    = grp::_1;
    constexpr auto ay    = grp::_2;
    constexpr auto az    = grp::_3;
    constexpr auto wa    = grp::_4;
    constexpr auto bx    = grp::_5;
    constexpr auto by    = grp::_6;
    constexpr auto bz    = grp::_7;
    constexpr auto wb    = grp::_8;
    constexpr auto alpha = grp::_9;

    constexpr auto dx  = bx - ax;
    constexpr auto dx2 = dx * dx;
    constexpr auto dy  = by - ay;
    constexpr auto dy2 = dy * dy;
    constexpr auto dz  = bz - az;
    constexpr auto dz2 = dz * dz;
    constexpr auto dw  = wb - wa;

    constexpr auto G = dx2 + dy2 + dz2;
    constexpr auto b = dx2 + dy2 + dz2 - dw;
    constexpr auto a1  = wa + alpha;
    constexpr auto a1_4 = a1 * int_const<4>{};

    constexpr auto expr = a1_4 * G - b * b;

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    // using stage_b     = grp::stage_b<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int orthow_n3_k2(double ax, double ay, double az, double wa,
                            double bx, double by, double bz, double wb,
                            double alpha)
{
    int orient_sign = ((ax == bx) && (ay == by) && (az == bz)) ? 0 : 1;
    if (orient_sign == 0) return 0;

    return - orthow_n3_k2_impl::pred{}.apply(ax, ay, az, wa,
                                             bx, by, bz, wb,
                                             alpha);
}



// ---------------------------------------------------------------------------
// orthow_n2_k3
// ---------------------------------------------------------------------------
namespace orthow_n3_k3_impl {
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
    constexpr auto alpha = grp::_13;

    constexpr auto dbx  = bx - ax;
    constexpr auto dbx2 = dbx * dbx;
    constexpr auto dby  = by - ay;
    constexpr auto dby2 = dby * dby;
    constexpr auto dbz  = bz - az;
    constexpr auto dbz2 = dbz * dbz;
    constexpr auto dwb  = wb - wa;

    constexpr auto dcx  = cx - ax;
    constexpr auto dcx2 = dcx * dcx;
    constexpr auto dcy  = cy - ay;
    constexpr auto dcy2 = dcy * dcy;
    constexpr auto dcz  = cz - az;
    constexpr auto dcz2 = dcz * dcz;
    constexpr auto dwc  = wc - wa;

    constexpr auto a1  = wa + alpha;
    constexpr auto a1_4 = a1 * int_const<4>{};

    constexpr auto b1 = dbx2 + dby2 + dbz2 - dwb;
    constexpr auto b2 = dcx2 + dcy2 + dcz2 - dwc;
    constexpr auto G11 = dbx2 + dby2 + dbz2;
    constexpr auto G12 = dbx * dcx + dby * dcy + dbz * dcz;
    constexpr auto G22 = dcx2 + dcy2 + dcz2;

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



extern "C" int orthow_n3_k3(double ax, double ay, double az, double wa,
                            double bx, double by, double bz, double wb,
                            double cx, double cy, double cz, double wc,
                            double alpha)
{
    double nx = (by-ay)*(cz-az) - (bz-az)*(cy-ay);
    double ny = (bz-az)*(cx-ax) - (bx-ax)*(cz-az);
    double nz = (bx-ax)*(cy-ay) - (by-ay)*(cx-ax);
    double a_nx = fabs(nx), a_ny = fabs(ny), a_nz = fabs(nz);
    int orient_sign;
    if (a_nx >= a_ny && a_nx >= a_nz)
        orient_sign = orient3d_impl::pred{}.apply(ax,ay,az, bx,by,bz, cx,cy,cz, ax+1,ay,az);
    else if (a_ny >= a_nz)
        orient_sign = orient3d_impl::pred{}.apply(ax,ay,az, bx,by,bz, cx,cy,cz, ax,ay+1,az);
    else
        orient_sign = orient3d_impl::pred{}.apply(ax,ay,az, bx,by,bz, cx,cy,cz, ax,ay,az+1);

    if (orient_sign == 0) return 0;

    return - orthow_n3_k3_impl::pred{}.apply(ax, ay, az, wa,
                                             bx, by, bz, wb,
                                             cx, cy, cz, wc,
                                             alpha);
}


// ---------------------------------------------------------------------------
// orthow_n3_k4
// ---------------------------------------------------------------------------
namespace orthow_n3_k4_impl {
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

    constexpr auto dbx  = bx - ax;
    constexpr auto dbx2 = dbx * dbx;
    constexpr auto dby  = by - ay;
    constexpr auto dby2 = dby * dby;
    constexpr auto dbz  = bz - az;
    constexpr auto dbz2 = dbz * dbz;
    constexpr auto dwb  = wb - wa;

    constexpr auto dcx  = cx - ax;
    constexpr auto dcx2 = dcx * dcx;
    constexpr auto dcy  = cy - ay;
    constexpr auto dcy2 = dcy * dcy;
    constexpr auto dcz  = cz - az;
    constexpr auto dcz2 = dcz * dcz;
    constexpr auto dwc  = wc - wa;

    constexpr auto ddx  = dx - ax;
    constexpr auto ddx2 = ddx * ddx;
    constexpr auto ddy  = dy - ay;
    constexpr auto ddy2 = ddy * ddy;
    constexpr auto ddz  = dz - az;
    constexpr auto ddz2 = ddz * ddz;
    constexpr auto dwd  = wd - wa;

    constexpr auto a1  = wa + alpha;
    constexpr auto a1_4 = a1 * int_const<4>{};

    constexpr auto b1 = dbx2 + dby2 + dbz2 - dwb;
    constexpr auto b2 = dcx2 + dcy2 + dcz2 - dwc;
    constexpr auto b3 = ddx2 + ddy2 + ddz2 - dwd;
    constexpr auto G11 = dbx2 + dby2 + dbz2;
    constexpr auto G12 = dbx * dcx + dby * dcy + dbz * dcz;
    constexpr auto G13 = dbx * ddx + dby * ddy + dbz * ddz;
    constexpr auto G22 = dcx2 + dcy2 + dcz2;
    constexpr auto G23 = dcx * ddx + dcy * ddy + dcz * ddz;
    constexpr auto G33 = ddx2 + ddy2 + ddz2;

    using expr_t = grp::det <
        decltype(a1_4), decltype(b1),  decltype(b2),  decltype(b3),
        decltype(b1),   decltype(G11), decltype(G12), decltype(G13),
        decltype(b2),   decltype(G12), decltype(G22), decltype(G23),
        decltype(b3),   decltype(G13), decltype(G23), decltype(G33)
    >;
    constexpr auto expr = expr_t{};

    using filter  = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact   = grp::stage_d<expr, double>;
    using stage_b = grp::stage_b<expr, double>;
    using pred    = grp::staged_predicate<filter, stage_b, exact>;
}



extern "C" int orthow_n3_k4(double ax, double ay, double az, double wa,
                            double bx, double by, double bz, double wb,
                            double cx, double cy, double cz, double wc,
                            double dx, double dy, double dz, double wd,
                            double alpha)
{
    int orient_sign = orient3d_impl::pred{}.apply(ax,ay,az, bx,by,bz, cx,cy,cz, dx,dy,dz);
    if (orient_sign == 0) return 0;

    return - orthow_n3_k4_impl::pred{}.apply(ax, ay, az, wa,
                                             bx, by, bz, wb,
                                             cx, cy, cz, wc,
                                             dx, dy, dz, wd,
                                             alpha);
}


template <std::size_t N>
struct [[deprecated("results_size — see template argument")]] show_stage_d_size {};
using _size_powertest_n3_k1 = show_stage_d_size<powertest_n3_k1_impl::exact::results_size>*;
using _size_powertest_n3_k2 = show_stage_d_size<powertest_n3_k2_D_xy_impl::exact::results_size>*;
using _size_powertest_n3_k3 = show_stage_d_size<powertest_n3_k3_impl::exact::results_size>*;
using _size_powertest_n3_k4 = show_stage_d_size<powertest_n3_k4_D_impl::exact::results_size>*;
using _size_orthow_n3_k1 = show_stage_d_size<orthow_n3_k1_impl::exact::results_size>*;
using _size_orthow_n3_k2 = show_stage_d_size<orthow_n3_k2_impl::exact::results_size>*;
using _size_orthow_n3_k3 = show_stage_d_size<orthow_n3_k3_impl::exact::results_size>*;
using _size_orthow_n3_k4 = show_stage_d_size<orthow_n3_k4_impl::exact::results_size>*;


// template <std::size_t N>
// struct [[deprecated("results_size — see template argument")]] show_stage_b_size {};
// using _size_powertest_n3_k1_b = show_stage_b_size<powertest_n3_k1_impl::stage_b::results_size>*;
// using _size_powertest_n3_k2_b = show_stage_b_size<powertest_n3_k2_D_xy_impl::stage_b::results_size>*;
// using _size_powertest_n3_k3_b = show_stage_b_size<powertest_n3_k3_impl::stage_b::results_size>*;
// using _size_powertest_n3_k4_b = show_stage_b_size<powertest_n3_k4_D_impl::stage_b::results_size>*;
// using _size_orthow_n3_k1_b = show_stage_b_size<orthow_n3_k1_impl::stage_b::results_size>*;
// using _size_orthow_n3_k2_b = show_stage_b_size<orthow_n3_k2_impl::stage_b::results_size>*;
// using _size_orthow_n3_k3_b = show_stage_b_size<orthow_n3_k3_impl::stage_b::results_size>*;
// using _size_orthow_n3_k4_b = show_stage_b_size<orthow_n3_k4_impl::stage_b::results_size>*;

