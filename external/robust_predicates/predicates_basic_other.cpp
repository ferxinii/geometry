#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_b.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;


// ---------------------------------------------------------------------------
// incircle3d -- is d on the circumcircle of triangle (a,b,c), for COPLANAR
// input.  = coplanar-reduced 3D insphere determinant of (a,b,c, w, d) with
// w = a + n, n = (b-a)x(c-a).  Translating by d (A=a-d,B=b-d,C=c-d, W=A+n) the
// insphere 4x4 is det[[A,SA],[B,SB],[C,SC],[W,SW]]; subtract row A from row W to
// get row (n, g), g=2A.n+|n|^2; split the offset column -> the g part factors as
// g*det[A,B,C] = g*(+-orient3d(a,b,c,d)) = 0 under the coplanarity precondition.
// The remaining part, expanded along the (now-zero-offset) squared-norm column:
//   expr = -SA*det[B,C,n] + SB*det[A,C,n] - SC*det[A,B,n]   (degree 6)
// zero iff d cospherical with {a,b,c,w}; with a,b,c,d coplanar and w off-plane
// that sphere meets the plane in circumcircle(a,b,c), so zero iff concyclic.
// Inputs: a(_1-_3), b(_4-_6), c(_7-_9), d(_10-_12)
// ---------------------------------------------------------------------------
namespace incircle3d_impl {
    constexpr auto ax = grp::_1;   constexpr auto ay = grp::_2;   constexpr auto az = grp::_3;
    constexpr auto bx = grp::_4;   constexpr auto by = grp::_5;   constexpr auto bz = grp::_6;
    constexpr auto cx = grp::_7;   constexpr auto cy = grp::_8;   constexpr auto cz = grp::_9;
    constexpr auto dx = grp::_10;  constexpr auto dy = grp::_11;  constexpr auto dz = grp::_12;

    // translate by d
    constexpr auto Ax = ax - dx;  constexpr auto Ay = ay - dy;  constexpr auto Az = az - dz;
    constexpr auto Bx = bx - dx;  constexpr auto By = by - dy;  constexpr auto Bz = bz - dz;
    constexpr auto Cx = cx - dx;  constexpr auto Cy = cy - dy;  constexpr auto Cz = cz - dz;
    // n = (b-a) x (c-a)
    constexpr auto ux = bx - ax;  constexpr auto uy = by - ay;  constexpr auto uz = bz - az;
    constexpr auto vx = cx - ax;  constexpr auto vy = cy - ay;  constexpr auto vz = cz - az;
    constexpr auto nx = uy*vz - uz*vy;
    constexpr auto ny = uz*vx - ux*vz;
    constexpr auto nz = ux*vy - uy*vx;
    // squared norms of the translated points
    constexpr auto SA = Ax*Ax + Ay*Ay + Az*Az;
    constexpr auto SB = Bx*Bx + By*By + Bz*Bz;
    constexpr auto SC = Cx*Cx + Cy*Cy + Cz*Cz;
    // 3x3 minors [., ., n]
    constexpr auto Mbcn = grp::det<decltype(Bx),decltype(By),decltype(Bz),
                                   decltype(Cx),decltype(Cy),decltype(Cz),
                                   decltype(nx),decltype(ny),decltype(nz)>{};
    constexpr auto Macn = grp::det<decltype(Ax),decltype(Ay),decltype(Az),
                                   decltype(Cx),decltype(Cy),decltype(Cz),
                                   decltype(nx),decltype(ny),decltype(nz)>{};
    constexpr auto Mabn = grp::det<decltype(Ax),decltype(Ay),decltype(Az),
                                   decltype(Bx),decltype(By),decltype(Bz),
                                   decltype(nx),decltype(ny),decltype(nz)>{};
    constexpr auto expr = SB*Macn - SA*Mbcn - SC*Mabn;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}


extern "C" int orient2d(double ax, double ay, double bx, double by, double cx, double cy);

extern "C" int incircle3d(double ax, double ay, double az,
                          double bx, double by, double bz,
                          double cx, double cy, double cz,
                          double dx, double dy, double dz)
{
    // a,b,c collinear (n = 0) => no circumcircle: the determinant would be a
    // spurious 0, so report the sentinel instead.  Collinear in 3D iff collinear
    // in all three axis projections.
    if (orient2d(ax,ay, bx,by, cx,cy) == 0 &&
        orient2d(ay,az, by,bz, cy,cz) == 0 &&
        orient2d(az,ax, bz,bx, cz,cx) == 0)
        return 2;
    return incircle3d_impl::pred{}.apply(ax,ay,az, bx,by,bz, cx,cy,cz, dx,dy,dz);
}


// ---------------------------------------------------------------------------
// orient3d_dd -- sign of det3(a-b, c-d, e-f): orientation of three explicit
// DIFFERENCE vectors.  Generalizes orient3d (orient3d(a,b,c,d) ==
// orient3d_dd(a,d, b,d, c,d)).  Needed by the weighted-DT sentinel reductions
// (vor3d PLAN_DE_PREDICATES.md): rows mix finite-point differences (F-o) with
// sentinel coordinate differences (s-t), which are not expressible as a plain
// orient3d of four points.  Degree 3, same filter behavior as orient3d.
// Inputs: a(_1-_3), b(_4-_6), c(_7-_9), d(_10-_12), e(_13-_15), f(_16-_18)
// ---------------------------------------------------------------------------
namespace orient3d_dd_impl {
    constexpr auto ax = grp::_1;   constexpr auto ay = grp::_2;   constexpr auto az = grp::_3;
    constexpr auto bx = grp::_4;   constexpr auto by = grp::_5;   constexpr auto bz = grp::_6;
    constexpr auto cx = grp::_7;   constexpr auto cy = grp::_8;   constexpr auto cz = grp::_9;
    constexpr auto dx = grp::_10;  constexpr auto dy = grp::_11;  constexpr auto dz = grp::_12;
    constexpr auto ex = grp::_13;  constexpr auto ey = grp::_14;  constexpr auto ez = grp::_15;
    constexpr auto fx = grp::_16;  constexpr auto fy = grp::_17;  constexpr auto fz = grp::_18;

    constexpr auto r0x = ax - bx;
    constexpr auto r0y = ay - by;
    constexpr auto r0z = az - bz;
    constexpr auto r1x = cx - dx;
    constexpr auto r1y = cy - dy;
    constexpr auto r1z = cz - dz;
    constexpr auto r2x = ex - fx;
    constexpr auto r2y = ey - fy;
    constexpr auto r2z = ez - fz;

    using expr_t = grp::det<
        decltype(r0x), decltype(r0y), decltype(r0z),
        decltype(r1x), decltype(r1y), decltype(r1z),
        decltype(r2x), decltype(r2y), decltype(r2z)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int orient3d_dd(double ax, double ay, double az,
                           double bx, double by, double bz,
                           double cx, double cy, double cz,
                           double dx, double dy, double dz,
                           double ex, double ey, double ez,
                           double fx, double fy, double fz)
{
    return orient3d_dd_impl::pred{}.apply(ax,ay,az, bx,by,bz, cx,cy,cz,
                                          dx,dy,dz, ex,ey,ez, fx,fy,fz);
}


#ifdef ROBUST_PREDICATES_PRINT_SIZE
template <std::size_t N>
struct [[deprecated("results_size — see template argument")]] show_stage_d_size {};
using _size_incircle3d = show_stage_d_size<incircle3d_impl::exact::results_size>*;
using _size_orient3d_dd = show_stage_d_size<orient3d_dd_impl::exact::results_size>*;
#endif


