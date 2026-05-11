#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;

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
    constexpr auto mu1 = mu1_h + mu1_h;
    // mu2 = 2*det(da, db, ex) = 2*(day*dbz - daz*dby)
    constexpr auto mu2_h = day*dbz - daz*dby;
    constexpr auto mu2 = mu2_h + mu2_h;

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
    constexpr auto mu1 = mu1_h + mu1_h;
    // mu2 = 2*det(da, db, ey) = 2*(daz*dbx - dax*dbz)
    constexpr auto mu2_h = daz*dbx - dax*dbz;
    constexpr auto mu2 = mu2_h + mu2_h;

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
    constexpr auto mu1 = mu1_h + mu1_h;
    // mu2 = 2*det(da, db, ez) = 2*(dax*dby - day*dbx)
    constexpr auto mu2_h = dax*dby - day*dbx;
    constexpr auto mu2 = mu2_h + mu2_h;

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
// powertest_n3_k1_alpha
// ---------------------------------------------------------------------------
namespace powertest_n3_k1_alpha_impl {
    constexpr auto ax    = grp::_1;
    constexpr auto ay    = grp::_2;
    constexpr auto az    = grp::_3;
    constexpr auto wa    = grp::_4;
    constexpr auto bx    = grp::_5;
    constexpr auto by    = grp::_6;
    constexpr auto bz    = grp::_7;
    constexpr auto wb    = grp::_8;
    constexpr auto alph1 = grp::_9;
    constexpr auto zero  = grp::_10;  // I need to make alpha a subtraction, 
                                      // otherwise if it is a leaf node, the
                                      // library crashes 

    constexpr auto dx   = bx - ax;
    constexpr auto dy   = by - ay;
    constexpr auto dz   = bz - az;
    constexpr auto dw = wa - wb;   // orthosphere has wa negative
    constexpr auto alpha = alph1 - zero;
    constexpr auto expr = dx*dx + dy*dy + dz*dz + dw - alpha;  

    using filter = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact  = grp::stage_d<expr, double>;
    using pred   = grp::staged_predicate<filter, exact>;
}

extern "C" int powertest_n3_k1_alpha(double ax, double ay, double az, double wa,
                                     double bx, double by, double bz, double wb,
                                     double alpha)
{
    return powertest_n3_k1_alpha_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, alpha, 0.0);
}


// ---------------------------------------------------------------------------
// powertest_n3_k2_alpha
// ---------------------------------------------------------------------------
namespace powertest_n3_k2_D_xy_alpha_impl {
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

    constexpr auto dax = ax - cx;
    constexpr auto day = ay - cy;
    constexpr auto daz = az - cz;
    constexpr auto dwa = wa - wc;
    constexpr auto dbx = bx - cx;
    constexpr auto dby = by - cy;
    constexpr auto dbz = bz - cz;
    constexpr auto dwb = wb - wc;

    constexpr auto la = dax*dax + day*day + daz*daz - dwa + alpha;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb + alpha;

    constexpr auto v1x = az - bz;
    constexpr auto v1z = bx - ax;
    constexpr auto v2y = bz - az;
    constexpr auto v2z = ay - by;

    constexpr auto mu1_h = daz*dbx - dax*dbz;
    constexpr auto mu1 = mu1_h + mu1_h;
    constexpr auto mu2_h = day*dbz - daz*dby;
    constexpr auto mu2 = mu2_h + mu2_h;

    constexpr auto det2_ab_xy   = dax*dby - dbx*day;
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
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace powertest_n3_k2_D_yz_alpha_impl {
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

    constexpr auto dax = ax - cx;
    constexpr auto day = ay - cy;
    constexpr auto daz = az - cz;
    constexpr auto dwa = wa - wc;
    constexpr auto dbx = bx - cx;
    constexpr auto dby = by - cy;
    constexpr auto dbz = bz - cz;
    constexpr auto dwb = wb - wc;

    constexpr auto la = dax*dax + day*day + daz*daz - dwa + alpha;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb + alpha;

    constexpr auto v1x = by - ay;
    constexpr auto v1y = ax - bx;
    constexpr auto v2x = az - bz;
    constexpr auto v2z = bx - ax;

    constexpr auto mu1_h = dax*dby - day*dbx;
    constexpr auto mu1 = mu1_h + mu1_h;
    constexpr auto mu2_h = daz*dbx - dax*dbz;
    constexpr auto mu2 = mu2_h + mu2_h;

    constexpr auto det2_ab_yz   = day*dbz - dby*daz;
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
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace powertest_n3_k2_D_zx_alpha_impl {
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

    constexpr auto dax = ax - cx;
    constexpr auto day = ay - cy;
    constexpr auto daz = az - cz;
    constexpr auto dwa = wa - wc;
    constexpr auto dbx = bx - cx;
    constexpr auto dby = by - cy;
    constexpr auto dbz = bz - cz;
    constexpr auto dwb = wb - wc;

    constexpr auto la = dax*dax + day*day + daz*daz - dwa + alpha;
    constexpr auto lb = dbx*dbx + dby*dby + dbz*dbz - dwb + alpha;

    constexpr auto v1y = bz - az;
    constexpr auto v1z = ay - by;
    constexpr auto v2x = by - ay;
    constexpr auto v2y = ax - bx;

    constexpr auto mu1_h = day*dbz - daz*dby;
    constexpr auto mu1 = mu1_h + mu1_h;
    constexpr auto mu2_h = dax*dby - day*dbx;
    constexpr auto mu2 = mu2_h + mu2_h;

    constexpr auto det2_ab_zx   = daz*dbx - dbz*dax;
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
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int powertest_n3_k2_alpha(double ax, double ay, double az, double wa,
                                     double bx, double by, double bz, double wb,
                                     double cx, double cy, double cz, double wc,
                                     double alpha)
{
    double abx = std::abs(bx - ax);
    double aby = std::abs(by - ay);
    double abz = std::abs(bz - az);

    int D_sign, orient_sign;

    if (abz >= abx && abz >= aby) {
        orient_sign = powertest_n3_k2_orient_xy_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc);
        if (orient_sign == 0) return 0;
        D_sign = powertest_n3_k2_D_xy_alpha_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc, alpha);
    } else if (abx >= aby && abx >= abz) {
        orient_sign = powertest_n3_k2_orient_yz_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc);
        if (orient_sign == 0) return 0;
        D_sign = powertest_n3_k2_D_yz_alpha_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc, alpha);
    } else {
        orient_sign = powertest_n3_k2_orient_zx_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc);
        if (orient_sign == 0) return 0;
        D_sign = powertest_n3_k2_D_zx_alpha_impl::pred{}.apply(ax, ay, az, wa, bx, by, bz, wb, cx, cy, cz, wc, alpha);
    }

    int out = orient_sign * D_sign;
    if (out > 0) return 1;
    else if (out < 0) return -1;
    else return 0;
}


