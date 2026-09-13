#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;

extern "C" int orient2d(double ax, double ay, double bx, double by, double cx, double cy);

// ---------------------------------------------------------------------------
// power_incircle3d -- weighted analogue of incircle3d (predicates_basic_other.cpp):
// is d co-orthocircular with the weighted triple (a,b,c), for COPLANAR input
// (a,b,c,d)?  Zero iff yes.  Same coplanar-reduction derivation as
// incircle3d (translate by d, build an off-plane witness W = A+n so the
// standard 4-point insphere determinant's witness-row contribution
// annihilates identically under the coplanarity precondition, det[A,B,C]=0)
// with ONE substitution: each defining point's squared-distance-from-d term
// SA=|A|^2 (etc.) becomes SA - wa + wd -- the same leaf substitution
// powertest_n3_k4 already uses for its own defining-point/query weight
// difference.  The witness's own weight never appears in the final formula
// (it multiplies the same det[A,B,C] term that vanishes under coplanarity,
// regardless of its value), so no witness weight is needed at all -- this
// generalizes incircle3d exactly, with no 2D projection anywhere.
// Returns 2 (sentinel, matching incircle3d) if a,b,c are collinear (no
// circumcircle to test against); this check is unweighted (collinearity is
// purely geometric).  Does NOT check a,b,c,d coplanarity itself -- same
// contract as incircle3d, the caller (power_points_concyclic) does that.
//
// LIBRARY BUG WORKAROUND (see predicates_lp_power.cpp): `<expansion> -
// <scalar>` doesn't compile for this framework (missing expansion_plus
// overload). Every weight substitution is written as a scalar-minus-scalar
// difference computed first (dwa = wd - wa), then added once (SA + dwa) --
// never SA - wa + wd.
//
// Inputs: a(_1-_3,wa=_4), b(_5-_7,wb=_8), c(_9-_11,wc=_12), d(_13-_15,wd=_16)
// ---------------------------------------------------------------------------
namespace power_incircle3d_impl {
    constexpr auto ax = grp::_1;   constexpr auto ay = grp::_2;   constexpr auto az = grp::_3;
    constexpr auto wa = grp::_4;
    constexpr auto bx = grp::_5;   constexpr auto by = grp::_6;   constexpr auto bz = grp::_7;
    constexpr auto wb = grp::_8;
    constexpr auto cx = grp::_9;   constexpr auto cy = grp::_10;  constexpr auto cz = grp::_11;
    constexpr auto wc = grp::_12;
    constexpr auto dx = grp::_13;  constexpr auto dy = grp::_14;  constexpr auto dz = grp::_15;
    constexpr auto wd = grp::_16;

    // translate by d (unweighted, geometric)
    constexpr auto Ax = ax - dx;  constexpr auto Ay = ay - dy;  constexpr auto Az = az - dz;
    constexpr auto Bx = bx - dx;  constexpr auto By = by - dy;  constexpr auto Bz = bz - dz;
    constexpr auto Cx = cx - dx;  constexpr auto Cy = cy - dy;  constexpr auto Cz = cz - dz;
    // n = (b-a) x (c-a) (unweighted, geometric)
    constexpr auto ux = bx - ax;  constexpr auto uy = by - ay;  constexpr auto uz = bz - az;
    constexpr auto vx = cx - ax;  constexpr auto vy = cy - ay;  constexpr auto vz = cz - az;
    constexpr auto nx = uy*vz - uz*vy;
    constexpr auto ny = uz*vx - ux*vz;
    constexpr auto nz = ux*vy - uy*vx;
    // raw squared norms of the translated points
    constexpr auto SA = Ax*Ax + Ay*Ay + Az*Az;
    constexpr auto SB = Bx*Bx + By*By + Bz*Bz;
    constexpr auto SC = Cx*Cx + Cy*Cy + Cz*Cz;
    // weighted substitution: SXw = SX - wx + wd, via a scalar diff computed
    // first (library-bug workaround).
    constexpr auto dwa = wd - wa;
    constexpr auto dwb = wd - wb;
    constexpr auto dwc = wd - wc;
    constexpr auto SAw = SA + dwa;
    constexpr auto SBw = SB + dwb;
    constexpr auto SCw = SC + dwc;
    // 3x3 minors [., ., n] -- pure orientation, unweighted, untouched.
    constexpr auto Mbcn = grp::det<decltype(Bx),decltype(By),decltype(Bz),
                                   decltype(Cx),decltype(Cy),decltype(Cz),
                                   decltype(nx),decltype(ny),decltype(nz)>{};
    constexpr auto Macn = grp::det<decltype(Ax),decltype(Ay),decltype(Az),
                                   decltype(Cx),decltype(Cy),decltype(Cz),
                                   decltype(nx),decltype(ny),decltype(nz)>{};
    constexpr auto Mabn = grp::det<decltype(Ax),decltype(Ay),decltype(Az),
                                   decltype(Bx),decltype(By),decltype(Bz),
                                   decltype(nx),decltype(ny),decltype(nz)>{};
    constexpr auto expr = SBw*Macn - SAw*Mbcn - SCw*Mabn;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_incircle3d(double ax, double ay, double az, double wa,
                                double bx, double by, double bz, double wb,
                                double cx, double cy, double cz, double wc,
                                double dx, double dy, double dz, double wd)
{
    // a,b,c collinear (n = 0): same sentinel and same (unweighted,
    // geometric) check as incircle3d.
    if (orient2d(ax,ay, bx,by, cx,cy) == 0 &&
        orient2d(ay,az, by,bz, cy,cz) == 0 &&
        orient2d(az,ax, bz,bx, cz,cx) == 0)
        return 2;
    return power_incircle3d_impl::pred{}.apply(ax,ay,az,wa, bx,by,bz,wb,
                                               cx,cy,cz,wc, dx,dy,dz,wd);
}
