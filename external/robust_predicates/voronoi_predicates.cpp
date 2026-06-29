#include <mpfr.h>
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;

// ---------------------------------------------------------------------------
// lp_feasible_T0_T1_S  -- orient_sign(TRI0, TRI1, seed_i)
// sign(|ti-A|^2 - |s-A|^2)
// Inputs: A(_1-_3), s(_4-_6), ti(_7-_9)
// ---------------------------------------------------------------------------
namespace lp_feasible_T0_T1_S_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto sx  = grp::_4;
    constexpr auto sy  = grp::_5;
    constexpr auto sz  = grp::_6;
    constexpr auto tix = grp::_7;
    constexpr auto tiy = grp::_8;
    constexpr auto tiz = grp::_9;

    constexpr auto dsx  = sx  - Ax;
    constexpr auto dsy  = sy  - Ay;
    constexpr auto dsz  = sz  - Az;
    constexpr auto dtix = tix - Ax;
    constexpr auto dtiy = tiy - Ay;
    constexpr auto dtiz = tiz - Az;

    constexpr auto expr = dtix*dtix + dtiy*dtiy + dtiz*dtiz
                        - dsx*dsx - dsy*dsy - dsz*dsz;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_feasible_T0_T1_S(double Ax, double Ay, double Az,
                                 double sx, double sy, double sz,
                                 double tix, double tiy, double tiz)
{
    return lp_feasible_T0_T1_S_impl::pred{}.apply(Ax, Ay, Az, sx, sy, sz, tix, tiy, tiz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T1_T2_S  -- orient_sign(TRI1, TRI2, seed_i)
// sign(|ti-B|^2 - |s-B|^2)
// Inputs: B(_1-_3), s(_4-_6), ti(_7-_9)
// ---------------------------------------------------------------------------
namespace lp_feasible_T1_T2_S_impl {
    constexpr auto Bx  = grp::_1;
    constexpr auto By  = grp::_2;
    constexpr auto Bz  = grp::_3;
    constexpr auto sx  = grp::_4;
    constexpr auto sy  = grp::_5;
    constexpr auto sz  = grp::_6;
    constexpr auto tix = grp::_7;
    constexpr auto tiy = grp::_8;
    constexpr auto tiz = grp::_9;

    constexpr auto dsx  = sx  - Bx;
    constexpr auto dsy  = sy  - By;
    constexpr auto dsz  = sz  - Bz;
    constexpr auto dtix = tix - Bx;
    constexpr auto dtiy = tiy - By;
    constexpr auto dtiz = tiz - Bz;

    constexpr auto expr = dtix*dtix + dtiy*dtiy + dtiz*dtiz
                        - dsx*dsx - dsy*dsy - dsz*dsz;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_feasible_T1_T2_S(double Bx, double By, double Bz,
                                 double sx, double sy, double sz,
                                 double tix, double tiy, double tiz)
{
    return lp_feasible_T1_T2_S_impl::pred{}.apply(Bx, By, Bz, sx, sy, sz, tix, tiy, tiz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T0_T2_S  -- orient_sign(TRI0, TRI2, seed_i)
// sign(|ti-C|^2 - |s-C|^2)
// Inputs: C(_1-_3), s(_4-_6), ti(_7-_9)
// ---------------------------------------------------------------------------
namespace lp_feasible_T0_T2_S_impl {
    constexpr auto Cx  = grp::_1;
    constexpr auto Cy  = grp::_2;
    constexpr auto Cz  = grp::_3;
    constexpr auto sx  = grp::_4;
    constexpr auto sy  = grp::_5;
    constexpr auto sz  = grp::_6;
    constexpr auto tix = grp::_7;
    constexpr auto tiy = grp::_8;
    constexpr auto tiz = grp::_9;

    constexpr auto dsx  = sx  - Cx;
    constexpr auto dsy  = sy  - Cy;
    constexpr auto dsz  = sz  - Cz;
    constexpr auto dtix = tix - Cx;
    constexpr auto dtiy = tiy - Cy;
    constexpr auto dtiz = tiz - Cz;

    constexpr auto expr = dtix*dtix + dtiy*dtiy + dtiz*dtiz
                        - dsx*dsx - dsy*dsy - dsz*dsz;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_feasible_T0_T2_S(double Cx, double Cy, double Cz,
                                 double sx, double sy, double sz,
                                 double tix, double tiy, double tiz)
{
    return lp_feasible_T0_T2_S_impl::pred{}.apply(Cx, Cy, Cz, sx, sy, sz, tix, tiy, tiz);
}


// ---------------------------------------------------------------------------
// lp_D_T0_S
// sign(-(ti-s)*(C-A))  =  sign(D(T0,i))  =  sign(-b_i)
// Inputs: A(_1-_3), C(_4-_6), s(_7-_9), ti(_10-_12)
// ---------------------------------------------------------------------------
namespace lp_D_T0_S_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Cx  = grp::_4;
    constexpr auto Cy  = grp::_5;
    constexpr auto Cz  = grp::_6;
    constexpr auto sx  = grp::_7;
    constexpr auto sy  = grp::_8;
    constexpr auto sz  = grp::_9;
    constexpr auto tix = grp::_10;
    constexpr auto tiy = grp::_11;
    constexpr auto tiz = grp::_12;

    constexpr auto dtx = sx - tix;
    constexpr auto dty = sy - tiy;
    constexpr auto dtz = sz - tiz;
    constexpr auto dCx = Cx - Ax;
    constexpr auto dCy = Cy - Ay;
    constexpr auto dCz = Cz - Az;

    constexpr auto expr = dtx*dCx + dty*dCy + dtz*dCz;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_D_T0_S(double Ax, double Ay, double Az,
                          double Cx, double Cy, double Cz,
                          double sx, double sy, double sz,
                          double tix, double tiy, double tiz)
{
    return lp_D_T0_S_impl::pred{}.apply(Ax, Ay, Az, Cx, Cy, Cz, sx, sy, sz, tix, tiy, tiz);
}


// ---------------------------------------------------------------------------
// lp_D_T1_S
// sign((ti-s)*(B-A))  =  sign(D(T1,i))  =  sign(a_i)
// Inputs: A(_1-_3), B(_4-_6), s(_7-_9), ti(_10-_12)
// ---------------------------------------------------------------------------
namespace lp_D_T1_S_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    constexpr auto sx  = grp::_7;
    constexpr auto sy  = grp::_8;
    constexpr auto sz  = grp::_9;
    constexpr auto tix = grp::_10;
    constexpr auto tiy = grp::_11;
    constexpr auto tiz = grp::_12;

    constexpr auto dtx = tix - sx;
    constexpr auto dty = tiy - sy;
    constexpr auto dtz = tiz - sz;
    constexpr auto dBx = Bx - Ax;
    constexpr auto dBy = By - Ay;
    constexpr auto dBz = Bz - Az;

    constexpr auto expr = dtx*dBx + dty*dBy + dtz*dBz;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_D_T1_S(double Ax, double Ay, double Az,
                          double Bx, double By, double Bz,
                          double sx, double sy, double sz,
                          double tix, double tiy, double tiz)
{
    return lp_D_T1_S_impl::pred{}.apply(Ax, Ay, Az, Bx, By, Bz, sx, sy, sz, tix, tiy, tiz);
}


// ---------------------------------------------------------------------------
// lp_D_T2_S
// sign((ti-s)*(C-B))  =  sign(D(T2,i))  =  sign(b_i - a_i)
// Inputs: B(_1-_3), C(_4-_6), s(_7-_9), ti(_10-_12)
// ---------------------------------------------------------------------------
namespace lp_D_T2_S_impl {
    constexpr auto Bx  = grp::_1;
    constexpr auto By  = grp::_2;
    constexpr auto Bz  = grp::_3;
    constexpr auto Cx  = grp::_4;
    constexpr auto Cy  = grp::_5;
    constexpr auto Cz  = grp::_6;
    constexpr auto sx  = grp::_7;
    constexpr auto sy  = grp::_8;
    constexpr auto sz  = grp::_9;
    constexpr auto tix = grp::_10;
    constexpr auto tiy = grp::_11;
    constexpr auto tiz = grp::_12;

    constexpr auto dtx = tix - sx;
    constexpr auto dty = tiy - sy;
    constexpr auto dtz = tiz - sz;
    constexpr auto dCBx = Cx - Bx;
    constexpr auto dCBy = Cy - By;
    constexpr auto dCBz = Cz - Bz;

    constexpr auto expr = dtx*dCBx + dty*dCBy + dtz*dCBz;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_D_T2_S(double Bx, double By, double Bz,
                          double Cx, double Cy, double Cz,
                          double sx, double sy, double sz,
                          double tix, double tiy, double tiz)
{
    return lp_D_T2_S_impl::pred{}.apply(Bx, By, Bz, Cx, Cy, Cz, sx, sy, sz, tix, tiy, tiz);
}


// ---------------------------------------------------------------------------
// lp_det2
// sign(a_j*b_k - a_k*b_j)
//   a_m = (tm-s).(B-A),  b_m = (tm-s).(C-A)
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9), s(_10-_12), tj(_13-_15), tk(_16-_18)
// ---------------------------------------------------------------------------
namespace lp_det2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    constexpr auto Cx  = grp::_7;
    constexpr auto Cy  = grp::_8;
    constexpr auto Cz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    constexpr auto tkx = grp::_16;
    constexpr auto tky = grp::_17;
    constexpr auto tkz = grp::_18;

    constexpr auto dBAx = Bx - Ax;
    constexpr auto dBAy = By - Ay;
    constexpr auto dBAz = Bz - Az;
    constexpr auto dCAx = Cx - Ax;
    constexpr auto dCAy = Cy - Ay;
    constexpr auto dCAz = Cz - Az;

    constexpr auto djx = tjx - sx;
    constexpr auto djy = tjy - sy;
    constexpr auto djz = tjz - sz;
    constexpr auto dkx = tkx - sx;
    constexpr auto dky = tky - sy;
    constexpr auto dkz = tkz - sz;

    constexpr auto a_j = djx*dBAx + djy*dBAy + djz*dBAz;
    constexpr auto b_j = djx*dCAx + djy*dCAy + djz*dCAz;
    constexpr auto a_k = dkx*dBAx + dky*dBAy + dkz*dBAz;
    constexpr auto b_k = dkx*dCAx + dky*dCAy + dkz*dCAz;

    using expr_t = grp::det<decltype(a_j), decltype(b_j),
                            decltype(a_k), decltype(b_k)>;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_det2(double Ax, double Ay, double Az,
                        double Bx, double By, double Bz,
                        double Cx, double Cy, double Cz,
                        double sx, double sy, double sz,
                        double tjx, double tjy, double tjz,
                        double tkx, double tky, double tkz)
{
    return lp_det2_impl::pred{}.apply(Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz,
                                      sx, sy, sz, tjx, tjy, tjz, tkx, tky, tkz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T0_S_S  -- orient_sign(TRI0, seed_j, seed_l)
// sign(b_j) * sign(b_j*c_l2 - c_j2*b_l)
//   b_m = (tm-s)*(C-A),  c_m2 = |tm-A|^2 - |s-A|^2
// Inputs: A(_1-_3), B(_4-_6, unused), C(_7-_9), s(_10-_12), tj(_13-_15), tl(_16-_18)
// ---------------------------------------------------------------------------
namespace lp_feasible_T0_S_S_factor1_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    // _4 _5 _6 = B, unused
    constexpr auto Cx  = grp::_7;
    constexpr auto Cy  = grp::_8;
    constexpr auto Cz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    // _16 _17 _18 = tl, unused

    constexpr auto dCAx = Cx - Ax;
    constexpr auto dCAy = Cy - Ay;
    constexpr auto dCAz = Cz - Az;
    constexpr auto djx  = tjx - sx;
    constexpr auto djy  = tjy - sy;
    constexpr auto djz  = tjz - sz;

    constexpr auto expr = djx*dCAx + djy*dCAy + djz*dCAz;  // b_j

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace lp_feasible_T0_S_S_factor2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    // _4 _5 _6 = B, unused
    constexpr auto Cx  = grp::_7;
    constexpr auto Cy  = grp::_8;
    constexpr auto Cz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    constexpr auto tlx = grp::_16;
    constexpr auto tly = grp::_17;
    constexpr auto tlz = grp::_18;

    constexpr auto dCAx  = Cx  - Ax;
    constexpr auto dCAy  = Cy  - Ay;
    constexpr auto dCAz  = Cz  - Az;
    constexpr auto djx   = tjx - sx;
    constexpr auto djy   = tjy - sy;
    constexpr auto djz   = tjz - sz;
    constexpr auto dlx   = tlx - sx;
    constexpr auto dly   = tly - sy;
    constexpr auto dlz   = tlz - sz;
    constexpr auto dsAx  = sx  - Ax;
    constexpr auto dsAy  = sy  - Ay;
    constexpr auto dsAz  = sz  - Az;
    constexpr auto dtjAx = tjx - Ax;
    constexpr auto dtjAy = tjy - Ay;
    constexpr auto dtjAz = tjz - Az;
    constexpr auto dtlAx = tlx - Ax;
    constexpr auto dtlAy = tly - Ay;
    constexpr auto dtlAz = tlz - Az;

    constexpr auto b_j = djx*dCAx + djy*dCAy + djz*dCAz;
    constexpr auto c_j = dtjAx*dtjAx + dtjAy*dtjAy + dtjAz*dtjAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;
    constexpr auto b_l = dlx*dCAx + dly*dCAy + dlz*dCAz;
    constexpr auto c_l = dtlAx*dtlAx + dtlAy*dtlAy + dtlAz*dtlAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;

    // b_j*c_l2 - c_j2*b_l
    using det2_t = grp::det<decltype(b_j), decltype(c_j),
                            decltype(b_l), decltype(c_l)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_feasible_T0_S_S(double Ax, double Ay, double Az,
                                double Bx, double By, double Bz,
                                double Cx, double Cy, double Cz,
                                double sx, double sy, double sz,
                                double tjx, double tjy, double tjz,
                                double tlx, double tly, double tlz)
{
    return lp_feasible_T0_S_S_factor1_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                                       sx,sy,sz, tjx,tjy,tjz, tlx,tly,tlz)
         * lp_feasible_T0_S_S_factor2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                                       sx,sy,sz, tjx,tjy,tjz, tlx,tly,tlz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T1_S_S  -- orient_sign(TRI1, seed_j, seed_l)
// sign(a_j) * sign(a_j*c_l2 - c_j2*a_l)
//   a_m = (tm-s)*(B-A),  c_m2 = |tm-A|^2 - |s-A|^2
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9, unused), s(_10-_12), tj(_13-_15), tl(_16-_18)
// ---------------------------------------------------------------------------
namespace lp_feasible_T1_S_S_factor1_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    // _7 _8 _9 = C, unused
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    // _16 _17 _18 = tl, unused

    constexpr auto dBAx = Bx - Ax;
    constexpr auto dBAy = By - Ay;
    constexpr auto dBAz = Bz - Az;
    constexpr auto djx  = tjx - sx;
    constexpr auto djy  = tjy - sy;
    constexpr auto djz  = tjz - sz;

    constexpr auto expr = djx*dBAx + djy*dBAy + djz*dBAz;  // a_j

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace lp_feasible_T1_S_S_factor2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    // _7 _8 _9 = C, unused
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    constexpr auto tlx = grp::_16;
    constexpr auto tly = grp::_17;
    constexpr auto tlz = grp::_18;

    constexpr auto dBAx  = Bx  - Ax;
    constexpr auto dBAy  = By  - Ay;
    constexpr auto dBAz  = Bz  - Az;
    constexpr auto djx   = tjx - sx;
    constexpr auto djy   = tjy - sy;
    constexpr auto djz   = tjz - sz;
    constexpr auto dlx   = tlx - sx;
    constexpr auto dly   = tly - sy;
    constexpr auto dlz   = tlz - sz;
    constexpr auto dsAx  = sx  - Ax;
    constexpr auto dsAy  = sy  - Ay;
    constexpr auto dsAz  = sz  - Az;
    constexpr auto dtjAx = tjx - Ax;
    constexpr auto dtjAy = tjy - Ay;
    constexpr auto dtjAz = tjz - Az;
    constexpr auto dtlAx = tlx - Ax;
    constexpr auto dtlAy = tly - Ay;
    constexpr auto dtlAz = tlz - Az;

    constexpr auto a_j = djx*dBAx + djy*dBAy + djz*dBAz;
    constexpr auto c_j = dtjAx*dtjAx + dtjAy*dtjAy + dtjAz*dtjAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;
    constexpr auto a_l = dlx*dBAx + dly*dBAy + dlz*dBAz;
    constexpr auto c_l = dtlAx*dtlAx + dtlAy*dtlAy + dtlAz*dtlAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;

    // a_j*c_l2 - c_j2*a_l
    using det2_t = grp::det<decltype(a_j), decltype(c_j),
                            decltype(a_l), decltype(c_l)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_feasible_T1_S_S(double Ax, double Ay, double Az,
                                double Bx, double By, double Bz,
                                double Cx, double Cy, double Cz,
                                double sx, double sy, double sz,
                                double tjx, double tjy, double tjz,
                                double tlx, double tly, double tlz)
{
    return lp_feasible_T1_S_S_factor1_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                                       sx,sy,sz, tjx,tjy,tjz, tlx,tly,tlz)
         * lp_feasible_T1_S_S_factor2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                                       sx,sy,sz, tjx,tjy,tjz, tlx,tly,tlz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T2_S_S  -- orient_sign(TRI2, seed_j, seed_l)
// sign(a_j-b_j) * sign(a_l*c_j2 - 2*a_l*b_j - b_l*c_j2 + 2*b_l*a_j + c_l2*b_j - c_l2*a_j)
//   a_m=(tm-s)*(B-A), b_m=(tm-s)*(C-A), c_m2=|tm-A|^2-|s-A|^2
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9), s(_10-_12), tj(_13-_15), tl(_16-_18)
// ---------------------------------------------------------------------------
namespace lp_feasible_T2_S_S_factor1_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    constexpr auto Cx  = grp::_7;
    constexpr auto Cy  = grp::_8;
    constexpr auto Cz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    // _16 _17 _18 = tl, unused

    constexpr auto dBAx = Bx - Ax;
    constexpr auto dBAy = By - Ay;
    constexpr auto dBAz = Bz - Az;
    constexpr auto dCAx = Cx - Ax;
    constexpr auto dCAy = Cy - Ay;
    constexpr auto dCAz = Cz - Az;
    constexpr auto djx  = tjx - sx;
    constexpr auto djy  = tjy - sy;
    constexpr auto djz  = tjz - sz;

    constexpr auto a_j = djx*dBAx + djy*dBAy + djz*dBAz;
    constexpr auto b_j = djx*dCAx + djy*dCAy + djz*dCAz;
    constexpr auto expr = b_j - a_j;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace lp_feasible_T2_S_S_factor2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    constexpr auto Cx  = grp::_7;
    constexpr auto Cy  = grp::_8;
    constexpr auto Cz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    constexpr auto tlx = grp::_16;
    constexpr auto tly = grp::_17;
    constexpr auto tlz = grp::_18;

    constexpr auto dBAx  = Bx  - Ax;
    constexpr auto dBAy  = By  - Ay;
    constexpr auto dBAz  = Bz  - Az;
    constexpr auto dCAx  = Cx  - Ax;
    constexpr auto dCAy  = Cy  - Ay;
    constexpr auto dCAz  = Cz  - Az;
    constexpr auto djx   = tjx - sx;
    constexpr auto djy   = tjy - sy;
    constexpr auto djz   = tjz - sz;
    constexpr auto dlx   = tlx - sx;
    constexpr auto dly   = tly - sy;
    constexpr auto dlz   = tlz - sz;
    constexpr auto dsAx  = sx  - Ax;
    constexpr auto dsAy  = sy  - Ay;
    constexpr auto dsAz  = sz  - Az;
    constexpr auto dtjAx = tjx - Ax;
    constexpr auto dtjAy = tjy - Ay;
    constexpr auto dtjAz = tjz - Az;
    constexpr auto dtlAx = tlx - Ax;
    constexpr auto dtlAy = tly - Ay;
    constexpr auto dtlAz = tlz - Az;

    constexpr auto a_j = djx*dBAx + djy*dBAy + djz*dBAz;
    constexpr auto b_j = djx*dCAx + djy*dCAy + djz*dCAz;
    constexpr auto c_j = dtjAx*dtjAx + dtjAy*dtjAy + dtjAz*dtjAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;
    constexpr auto a_l = dlx*dBAx + dly*dBAy + dlz*dBAz;
    constexpr auto b_l = dlx*dCAx + dly*dCAy + dlz*dCAz;
    constexpr auto c_l = dtlAx*dtlAx + dtlAy*dtlAy + dtlAz*dtlAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;

    // a_l*c_j2 - 2*a_l*b_j - b_l*c_j2 + 2*b_l*a_j + c_l2*b_j - c_l2*a_j
    constexpr auto expr = a_l*c_j - (a_l*b_j + a_l*b_j)
                        - b_l*c_j + (b_l*a_j + b_l*a_j)
                        + c_l*b_j - c_l*a_j;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_feasible_T2_S_S(double Ax, double Ay, double Az,
                                double Bx, double By, double Bz,
                                double Cx, double Cy, double Cz,
                                double sx, double sy, double sz,
                                double tjx, double tjy, double tjz,
                                double tlx, double tly, double tlz)
{
    return lp_feasible_T2_S_S_factor1_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                                       sx,sy,sz, tjx,tjy,tjz, tlx,tly,tlz)
         * lp_feasible_T2_S_S_factor2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                                       sx,sy,sz, tjx,tjy,tjz, tlx,tly,tlz);
}


// ---------------------------------------------------------------------------
// lp_det3
// sign(det | a_j b_j c_j |)
//          | a_k b_k c_k |
//          | a_l b_l c_l |
//   a_m=(tm-s).(B-A), b_m=(tm-s).(C-A), c_m=|tm-A|^2-|s-A|^2
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9), s(_10-_12),
//         tj(_13-_15), tk(_16-_18), tl(_19-_21)
// ---------------------------------------------------------------------------
namespace lp_det3_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    constexpr auto Cx  = grp::_7;
    constexpr auto Cy  = grp::_8;
    constexpr auto Cz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    constexpr auto tkx = grp::_16;
    constexpr auto tky = grp::_17;
    constexpr auto tkz = grp::_18;
    constexpr auto tlx = grp::_19;
    constexpr auto tly = grp::_20;
    constexpr auto tlz = grp::_21;

    constexpr auto dBAx = Bx - Ax;
    constexpr auto dBAy = By - Ay;
    constexpr auto dBAz = Bz - Az;
    constexpr auto dCAx = Cx - Ax;
    constexpr auto dCAy = Cy - Ay;
    constexpr auto dCAz = Cz - Az;

    constexpr auto dsAx  = sx  - Ax;
    constexpr auto dsAy  = sy  - Ay;
    constexpr auto dsAz  = sz  - Az;

    constexpr auto djx = tjx - sx;
    constexpr auto djy = tjy - sy;
    constexpr auto djz = tjz - sz;
    constexpr auto dkx = tkx - sx;
    constexpr auto dky = tky - sy;
    constexpr auto dkz = tkz - sz;
    constexpr auto dlx = tlx - sx;
    constexpr auto dly = tly - sy;
    constexpr auto dlz = tlz - sz;

    constexpr auto dtjAx = tjx - Ax;
    constexpr auto dtjAy = tjy - Ay;
    constexpr auto dtjAz = tjz - Az;
    constexpr auto dtkAx = tkx - Ax;
    constexpr auto dtkAy = tky - Ay;
    constexpr auto dtkAz = tkz - Az;
    constexpr auto dtlAx = tlx - Ax;
    constexpr auto dtlAy = tly - Ay;
    constexpr auto dtlAz = tlz - Az;

    constexpr auto a_j = djx*dBAx + djy*dBAy + djz*dBAz;
    constexpr auto b_j = djx*dCAx + djy*dCAy + djz*dCAz;
    constexpr auto c_j = dtjAx*dtjAx + dtjAy*dtjAy + dtjAz*dtjAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;
    constexpr auto a_k = dkx*dBAx + dky*dBAy + dkz*dBAz;
    constexpr auto b_k = dkx*dCAx + dky*dCAy + dkz*dCAz;
    constexpr auto c_k = dtkAx*dtkAx + dtkAy*dtkAy + dtkAz*dtkAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;
    constexpr auto a_l = dlx*dBAx + dly*dBAy + dlz*dBAz;
    constexpr auto b_l = dlx*dCAx + dly*dCAy + dlz*dCAz;
    constexpr auto c_l = dtlAx*dtlAx + dtlAy*dtlAy + dtlAz*dtlAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;

    using expr_t = grp::det<
        decltype(a_j), decltype(b_j), decltype(c_j),
        decltype(a_k), decltype(b_k), decltype(c_k),
        decltype(a_l), decltype(b_l), decltype(c_l)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_det3(double Ax, double Ay, double Az,
                        double Bx, double By, double Bz,
                        double Cx, double Cy, double Cz,
                        double sx, double sy, double sz,
                        double tjx, double tjy, double tjz,
                        double tkx, double tky, double tkz,
                        double tlx, double tly, double tlz)
{
    return lp_det3_impl::pred{}.apply(Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz,
                                      sx, sy, sz, tjx, tjy, tjz, tkx, tky, tkz,
                                      tlx, tly, tlz);
}



// ---------------------------------------------------------------------------
// lp_feasible_T0_S_T2  -- orient_sign(TRI0, seed_j, TRI2)
// sign(-b_j) * sign(|tj-C|^2 - |s-C|^2)
//   = lp_D_T0_S * feasible_T0_T2_S
// Inputs: A(_1-_3), C(_4-_6), s(_7-_9), tj(_10-_12)
// ---------------------------------------------------------------------------
extern "C" int lp_feasible_T0_S_T2(double Ax, double Ay, double Az,
                                 double Cx, double Cy, double Cz,
                                 double sx, double sy, double sz,
                                 double tjx, double tjy, double tjz)
{
    return lp_D_T0_S_impl::pred{}.apply(Ax,Ay,Az, Cx,Cy,Cz, sx,sy,sz, tjx,tjy,tjz)
         * lp_feasible_T0_T2_S_impl::pred{}.apply(Cx,Cy,Cz, sx,sy,sz, tjx,tjy,tjz);
}

// ---------------------------------------------------------------------------
// lp_feasible_S_S_T0  -- orient_sign(seed_j, seed_k, TRI0)
// sign(a_j*b_k-a_k*b_j) * sign(b_j*c_k2-c_j2*b_k)
//   a_m=(tm-s)*(B-A), b_m=(tm-s)*(C-A), c_m2=|tm-A|^2-|s-A|^2
// factor1 reuses lp_det2_impl
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9), s(_10-_12), tj(_13-_15), tk(_16-_18)
// ---------------------------------------------------------------------------
namespace lp_feasible_S_S_T0_factor2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    // _4 _5 _6 = B, unused
    constexpr auto Cx  = grp::_7;
    constexpr auto Cy  = grp::_8;
    constexpr auto Cz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    constexpr auto tkx = grp::_16;
    constexpr auto tky = grp::_17;
    constexpr auto tkz = grp::_18;

    constexpr auto dCAx = Cx - Ax;
    constexpr auto dCAy = Cy - Ay;
    constexpr auto dCAz = Cz - Az;

    constexpr auto djx = tjx - sx;
    constexpr auto djy = tjy - sy;
    constexpr auto djz = tjz - sz;
    constexpr auto dkx = tkx - sx;
    constexpr auto dky = tky - sy;
    constexpr auto dkz = tkz - sz;

    constexpr auto dsAx  = sx  - Ax;
    constexpr auto dsAy  = sy  - Ay;
    constexpr auto dsAz  = sz  - Az;
    constexpr auto dtjAx = tjx - Ax;
    constexpr auto dtjAy = tjy - Ay;
    constexpr auto dtjAz = tjz - Az;
    constexpr auto dtkAx = tkx - Ax;
    constexpr auto dtkAy = tky - Ay;
    constexpr auto dtkAz = tkz - Az;

    constexpr auto b_j = djx*dCAx + djy*dCAy + djz*dCAz;
    constexpr auto c_j = dtjAx*dtjAx + dtjAy*dtjAy + dtjAz*dtjAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;
    constexpr auto b_k = dkx*dCAx + dky*dCAy + dkz*dCAz;
    constexpr auto c_k = dtkAx*dtkAx + dtkAy*dtkAy + dtkAz*dtkAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;

    // -(b_j*c_k2 - c_j2*b_k)
    using det2_t = grp::det<decltype(b_k), decltype(c_k),
                            decltype(b_j), decltype(c_j)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_feasible_S_S_T0(double Ax, double Ay, double Az,
                                double Bx, double By, double Bz,
                                double Cx, double Cy, double Cz,
                                double sx, double sy, double sz,
                                double tjx, double tjy, double tjz,
                                double tkx, double tky, double tkz)
{
    return lp_det2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                      sx,sy,sz, tjx,tjy,tjz, tkx,tky,tkz)
         * lp_feasible_S_S_T0_factor2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                                       sx,sy,sz, tjx,tjy,tjz, tkx,tky,tkz);
}


// ---------------------------------------------------------------------------
// lp_feasible_S_S_T2  -- orient_sign(seed_j, seed_k, TRI2)
// sign(a_j*b_k-a_k*b_j) * sign(c_j2*b_k-b_j*c_k2+a_j*c_k2-c_j2*a_k-2*a_j*b_k+2*b_j*a_k)
//   a_m=(tm-s)*(B-A), b_m=(tm-s)*(C-A), c_m2=|tm-A|^2-|s-A|^2
// factor1 reuses lp_det2_impl
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9), s(_10-_12), tj(_13-_15), tk(_16-_18)
// ---------------------------------------------------------------------------
namespace lp_feasible_S_S_T2_factor2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    constexpr auto Cx  = grp::_7;
    constexpr auto Cy  = grp::_8;
    constexpr auto Cz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    constexpr auto tkx = grp::_16;
    constexpr auto tky = grp::_17;
    constexpr auto tkz = grp::_18;

    constexpr auto dBAx = Bx - Ax;
    constexpr auto dBAy = By - Ay;
    constexpr auto dBAz = Bz - Az;
    constexpr auto dCAx = Cx - Ax;
    constexpr auto dCAy = Cy - Ay;
    constexpr auto dCAz = Cz - Az;

    constexpr auto djx = tjx - sx;
    constexpr auto djy = tjy - sy;
    constexpr auto djz = tjz - sz;
    constexpr auto dkx = tkx - sx;
    constexpr auto dky = tky - sy;
    constexpr auto dkz = tkz - sz;

    constexpr auto dsAx  = sx  - Ax;
    constexpr auto dsAy  = sy  - Ay;
    constexpr auto dsAz  = sz  - Az;
    constexpr auto dtjAx = tjx - Ax;
    constexpr auto dtjAy = tjy - Ay;
    constexpr auto dtjAz = tjz - Az;
    constexpr auto dtkAx = tkx - Ax;
    constexpr auto dtkAy = tky - Ay;
    constexpr auto dtkAz = tkz - Az;

    constexpr auto a_j = djx*dBAx + djy*dBAy + djz*dBAz;
    constexpr auto b_j = djx*dCAx + djy*dCAy + djz*dCAz;
    constexpr auto c_j = dtjAx*dtjAx + dtjAy*dtjAy + dtjAz*dtjAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;
    constexpr auto a_k = dkx*dBAx + dky*dBAy + dkz*dBAz;
    constexpr auto b_k = dkx*dCAx + dky*dCAy + dkz*dCAz;
    constexpr auto c_k = dtkAx*dtkAx + dtkAy*dtkAy + dtkAz*dtkAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;

    // -(c_j2*b_k - b_j*c_k2 + a_j*c_k2 - c_j2*a_k - 2*a_j*b_k + 2*b_j*a_k)
    constexpr auto expr = b_j*c_k - c_j*b_k - a_j*c_k + c_j*a_k
                        + (a_j*b_k + a_j*b_k) - (b_j*a_k + b_j*a_k);

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_feasible_S_S_T2(double Ax, double Ay, double Az,
                                double Bx, double By, double Bz,
                                double Cx, double Cy, double Cz,
                                double sx, double sy, double sz,
                                double tjx, double tjy, double tjz,
                                double tkx, double tky, double tkz)
{
    return lp_det2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                      sx,sy,sz, tjx,tjy,tjz, tkx,tky,tkz)
         * lp_feasible_S_S_T2_factor2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                                       sx,sy,sz, tjx,tjy,tjz, tkx,tky,tkz);
}


// ---------------------------------------------------------------------------
// lp_feasible_S_S_T1  -- orient_sign(seed_j, seed_k, TRI1)
// sign(a_j*b_k-a_k*b_j) * sign(c_j2*a_k-a_j*c_k2)
//   a_m=(tm-s)*(B-A), c_m2=|tm-A|^2-|s-A|^2
// factor1 reuses lp_det2_impl
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9, unused), s(_10-_12), tj(_13-_15), tk(_16-_18)
// ---------------------------------------------------------------------------
namespace lp_feasible_S_S_T1_factor2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    // _7 _8 _9 = C, unused
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto tjx = grp::_13;
    constexpr auto tjy = grp::_14;
    constexpr auto tjz = grp::_15;
    constexpr auto tkx = grp::_16;
    constexpr auto tky = grp::_17;
    constexpr auto tkz = grp::_18;

    constexpr auto dBAx = Bx - Ax;
    constexpr auto dBAy = By - Ay;
    constexpr auto dBAz = Bz - Az;

    constexpr auto djx = tjx - sx;
    constexpr auto djy = tjy - sy;
    constexpr auto djz = tjz - sz;
    constexpr auto dkx = tkx - sx;
    constexpr auto dky = tky - sy;
    constexpr auto dkz = tkz - sz;

    constexpr auto dsAx  = sx  - Ax;
    constexpr auto dsAy  = sy  - Ay;
    constexpr auto dsAz  = sz  - Az;
    constexpr auto dtjAx = tjx - Ax;
    constexpr auto dtjAy = tjy - Ay;
    constexpr auto dtjAz = tjz - Az;
    constexpr auto dtkAx = tkx - Ax;
    constexpr auto dtkAy = tky - Ay;
    constexpr auto dtkAz = tkz - Az;

    constexpr auto a_j = djx*dBAx + djy*dBAy + djz*dBAz;
    constexpr auto c_j = dtjAx*dtjAx + dtjAy*dtjAy + dtjAz*dtjAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;
    constexpr auto a_k = dkx*dBAx + dky*dBAy + dkz*dBAz;
    constexpr auto c_k = dtkAx*dtkAx + dtkAy*dtkAy + dtkAz*dtkAz
                       - dsAx*dsAx   - dsAy*dsAy   - dsAz*dsAz;

    // -(c_j2*a_k - a_j*c_k2)
    using det2_t = grp::det<decltype(c_k), decltype(a_k),
                            decltype(c_j), decltype(a_j)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp_feasible_S_S_T1(double Ax, double Ay, double Az,
                                double Bx, double By, double Bz,
                                double Cx, double Cy, double Cz,
                                double sx, double sy, double sz,
                                double tjx, double tjy, double tjz,
                                double tkx, double tky, double tkz)
{
    return lp_det2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                      sx,sy,sz, tjx,tjy,tjz, tkx,tky,tkz)
         * lp_feasible_S_S_T1_factor2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                                                       sx,sy,sz, tjx,tjy,tjz, tkx,tky,tkz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T1_S_T2  -- orient_sign(TRI1, seed_j, TRI2)
// -lp_D_T1_S * feasible_T1_T2_S
//   = -sign(a_j) * sign(|tj-B|^2 - |s-B|^2)
// Inputs: A(_1-_3), B(_4-_6), s(_7-_9), tj(_10-_12)
// ---------------------------------------------------------------------------
extern "C" int lp_feasible_T1_S_T2(double Ax, double Ay, double Az,
                                 double Bx, double By, double Bz,
                                 double sx, double sy, double sz,
                                 double tjx, double tjy, double tjz)
{
    return -lp_D_T1_S_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, sx,sy,sz, tjx,tjy,tjz)
         *  lp_feasible_T1_T2_S_impl::pred{}.apply(Bx,By,Bz, sx,sy,sz, tjx,tjy,tjz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T2_S_T1  -- orient_sign(TRI2, seed_j, TRI1)
// lp_D_T2_S * feasible_T1_T2_S
//   = sign(b_j - a_j) * sign(|tj-B|^2 - |s-B|^2)
// Inputs: B(_1-_3), C(_4-_6), s(_7-_9), tj(_10-_12)
// ---------------------------------------------------------------------------
extern "C" int lp_feasible_T2_S_T1(double Bx, double By, double Bz,
                                 double Cx, double Cy, double Cz,
                                 double sx, double sy, double sz,
                                 double tjx, double tjy, double tjz)
{
    return lp_D_T2_S_impl::pred{}.apply(Bx,By,Bz, Cx,Cy,Cz, sx,sy,sz, tjx,tjy,tjz)
         * lp_feasible_T1_T2_S_impl::pred{}.apply(Bx,By,Bz, sx,sy,sz, tjx,tjy,tjz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T2_S_T0  -- orient_sign(TRI2, seed_j, TRI0)
// -lp_D_T2_S * feasible_T0_T2_S
//   = -sign(b_j - a_j) * sign(|tj-C|^2 - |s-C|^2)
// Inputs: B(_1-_3), C(_4-_6), s(_7-_9), tj(_10-_12)
// ---------------------------------------------------------------------------
extern "C" int lp_feasible_T2_S_T0(double Bx, double By, double Bz,
                                 double Cx, double Cy, double Cz,
                                 double sx, double sy, double sz,
                                 double tjx, double tjy, double tjz)
{
    return -lp_D_T2_S_impl::pred{}.apply(Bx,By,Bz, Cx,Cy,Cz, sx,sy,sz, tjx,tjy,tjz)
         *  lp_feasible_T0_T2_S_impl::pred{}.apply(Cx,Cy,Cz, sx,sy,sz, tjx,tjy,tjz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T0_S_T1  -- orient_sign(TRI0, seed_j, TRI1)
// -lp_D_T0_S * feasible_T0_T1_S
//   = -sign(-b_j) * sign(|tj-A|^2 - |s-A|^2)
// Inputs: A(_1-_3), C(_4-_6), s(_7-_9), tj(_10-_12)
// ---------------------------------------------------------------------------
extern "C" int lp_feasible_T0_S_T1(double Ax, double Ay, double Az,
                                 double Cx, double Cy, double Cz,
                                 double sx, double sy, double sz,
                                 double tjx, double tjy, double tjz)
{
    return -lp_D_T0_S_impl::pred{}.apply(Ax,Ay,Az, Cx,Cy,Cz, sx,sy,sz, tjx,tjy,tjz)
         *  lp_feasible_T0_T1_S_impl::pred{}.apply(Ax,Ay,Az, sx,sy,sz, tjx,tjy,tjz);
}


// ---------------------------------------------------------------------------
// lp_feasible_T1_S_T0  -- orient_sign(TRI1, seed_j, TRI0)
// lp_D_T1_S * feasible_T0_T1_S
//   = sign(a_j) * sign(|tj-A|^2 - |s-A|^2)
// Inputs: A(_1-_3), B(_4-_6), s(_7-_9), tj(_10-_12)
// ---------------------------------------------------------------------------
extern "C" int lp_feasible_T1_S_T0(double Ax, double Ay, double Az,
                                 double Bx, double By, double Bz,
                                 double sx, double sy, double sz,
                                 double tjx, double tjy, double tjz)
{
    return lp_D_T1_S_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz, sx,sy,sz, tjx,tjy,tjz)
         * lp_feasible_T0_T1_S_impl::pred{}.apply(Ax,Ay,Az, sx,sy,sz, tjx,tjy,tjz);
}


// ---------------------------------------------------------------------------
// cdt_dprime
// sign of D' -- the 4x4 determinant:
//   | a-Rv  dot(Rv-a, n) |   rows for {a, b, c, L}
//   | b-Rv  dot(Rv-b, n) |   n = (pb-pa) x (pc-pa)  (unnormalized)
//   | c-Rv  dot(Rv-c, n) |
//   | L-Rv  dot(Rv-L, n) |
// Inputs: a(_1-_3), b(_4-_6), c(_7-_9), L(_10-_12),
//         Rv(_13-_15), pa(_16-_18), pb(_19-_21), pc(_22-_24)
// ---------------------------------------------------------------------------
namespace cdt_dprime_impl {
    constexpr auto ax  = grp::_1;
    constexpr auto ay  = grp::_2;
    constexpr auto az  = grp::_3;
    constexpr auto bx  = grp::_4;
    constexpr auto by  = grp::_5;
    constexpr auto bz  = grp::_6;
    constexpr auto cx  = grp::_7;
    constexpr auto cy  = grp::_8;
    constexpr auto cz  = grp::_9;
    constexpr auto Lx  = grp::_10;
    constexpr auto Ly  = grp::_11;
    constexpr auto Lz  = grp::_12;
    constexpr auto Rvx = grp::_13;
    constexpr auto Rvy = grp::_14;
    constexpr auto Rvz = grp::_15;
    constexpr auto pax = grp::_16;
    constexpr auto pay = grp::_17;
    constexpr auto paz = grp::_18;
    constexpr auto pbx = grp::_19;
    constexpr auto pby = grp::_20;
    constexpr auto pbz = grp::_21;
    constexpr auto pcx = grp::_22;
    constexpr auto pcy = grp::_23;
    constexpr auto pcz = grp::_24;

    // first 3 columns: pts[i] - Rv
    constexpr auto dax = ax - Rvx;
    constexpr auto day = ay - Rvy;
    constexpr auto daz = az - Rvz;
    constexpr auto dbx = bx - Rvx;
    constexpr auto dby = by - Rvy;
    constexpr auto dbz = bz - Rvz;
    constexpr auto dcx = cx - Rvx;
    constexpr auto dcy = cy - Rvy;
    constexpr auto dcz = cz - Rvz;
    constexpr auto dLx = Lx - Rvx;
    constexpr auto dLy = Ly - Rvy;
    constexpr auto dLz = Lz - Rvz;

    // unnormalized face normal n = (pb-pa) x (pc-pa)
    constexpr auto dpbpax = pbx - pax;
    constexpr auto dpbpay = pby - pay;
    constexpr auto dpbpaz = pbz - paz;
    constexpr auto dpcpax = pcx - pax;
    constexpr auto dpcpay = pcy - pay;
    constexpr auto dpcpaz = pcz - paz;

    // negative normal
    constexpr auto nnx = dpbpaz * dpcpay - dpbpay * dpcpaz;
    constexpr auto nny = dpbpax * dpcpaz - dpbpaz * dpcpax;
    constexpr auto nnz = dpbpay * dpcpax - dpbpax * dpcpay;

    // 4th column: dot(Rv - pts[i], n) = -dot(pts[i]-Rv, n) = dot(., -n)
    constexpr auto col4_a = dax * nnx + day * nny + daz * nnz;
    constexpr auto col4_b = dbx * nnx + dby * nny + dbz * nnz;
    constexpr auto col4_c = dcx * nnx + dcy * nny + dcz * nnz;
    constexpr auto col4_L = dLx * nnx + dLy * nny + dLz * nnz;

    using expr_t = grp::det<
        decltype(dax), decltype(day), decltype(daz), decltype(col4_a),
        decltype(dbx), decltype(dby), decltype(dbz), decltype(col4_b),
        decltype(dcx), decltype(dcy), decltype(dcz), decltype(col4_c),
        decltype(dLx), decltype(dLy), decltype(dLz), decltype(col4_L)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int cdt_dprime(double ax,  double ay,  double az,
                           double bx,  double by,  double bz,
                           double cx,  double cy,  double cz,
                           double Lx,  double Ly,  double Lz,
                           double Rvx, double Rvy, double Rvz,
                           double pax, double pay, double paz,
                           double pbx, double pby, double pbz,
                           double pcx, double pcy, double pcz)
{
    return cdt_dprime_impl::pred{}.apply(
        ax,  ay,  az,
        bx,  by,  bz,
        cx,  cy,  cz,
        Lx,  Ly,  Lz,
        Rvx, Rvy, Rvz,
        pax, pay, paz,
        pbx, pby, pbz,
        pcx, pcy, pcz);
}


// ---------------------------------------------------------------------------
// MPFR helpers for cdt_sign_cross_dprime exact stage.
// All computations at prec bits; for prec=640 every intermediate result from
// double inputs is represented exactly (degree-11 poly in 53-bit inputs needs
// at most 11*53 + ceil(log2(1152)) = 594 bits; 640 gives a 46-bit margin).
// ---------------------------------------------------------------------------
static void mpfr_det3(mpfr_t result,
    mpfr_t m00, mpfr_t m01, mpfr_t m02,
    mpfr_t m10, mpfr_t m11, mpfr_t m12,
    mpfr_t m20, mpfr_t m21, mpfr_t m22,
    mpfr_prec_t prec)
{
    mpfr_t t0, t1, tmp;
    mpfr_init2(t0,  prec);
    mpfr_init2(t1,  prec);
    mpfr_init2(tmp, prec);

    // result = m00*(m11*m22 - m12*m21)
    //        - m01*(m10*m22 - m12*m20)
    //        + m02*(m10*m21 - m11*m20)

    mpfr_mul(t0, m11, m22, MPFR_RNDN);
    mpfr_mul(t1, m12, m21, MPFR_RNDN);
    mpfr_sub(tmp, t0, t1, MPFR_RNDN);
    mpfr_mul(result, m00, tmp, MPFR_RNDN);

    mpfr_mul(t0, m10, m22, MPFR_RNDN);
    mpfr_mul(t1, m12, m20, MPFR_RNDN);
    mpfr_sub(tmp, t0, t1, MPFR_RNDN);
    mpfr_mul(tmp, m01, tmp, MPFR_RNDN);
    mpfr_sub(result, result, tmp, MPFR_RNDN);

    mpfr_mul(t0, m10, m21, MPFR_RNDN);
    mpfr_mul(t1, m11, m20, MPFR_RNDN);
    mpfr_sub(tmp, t0, t1, MPFR_RNDN);
    mpfr_mul(tmp, m02, tmp, MPFR_RNDN);
    mpfr_add(result, result, tmp, MPFR_RNDN);

    mpfr_clear(t0);
    mpfr_clear(t1);
    mpfr_clear(tmp);
}

// Cofactor expansion along row 0. m[i][j] must be initialised at prec bits.
static void mpfr_det4(mpfr_t result, mpfr_t m[4][4], mpfr_prec_t prec)
{
    mpfr_t t0, t1;
    mpfr_init2(t0, prec);
    mpfr_init2(t1, prec);

    // j=0: +m[0][0] * det3(cols 1,2,3 of rows 1,2,3)
    mpfr_det3(t0,
        m[1][1], m[1][2], m[1][3],
        m[2][1], m[2][2], m[2][3],
        m[3][1], m[3][2], m[3][3], prec);
    mpfr_mul(result, m[0][0], t0, MPFR_RNDN);

    // j=1: -m[0][1] * det3(cols 0,2,3 of rows 1,2,3)
    mpfr_det3(t0,
        m[1][0], m[1][2], m[1][3],
        m[2][0], m[2][2], m[2][3],
        m[3][0], m[3][2], m[3][3], prec);
    mpfr_mul(t1, m[0][1], t0, MPFR_RNDN);
    mpfr_sub(result, result, t1, MPFR_RNDN);

    // j=2: +m[0][2] * det3(cols 0,1,3 of rows 1,2,3)
    mpfr_det3(t0,
        m[1][0], m[1][1], m[1][3],
        m[2][0], m[2][1], m[2][3],
        m[3][0], m[3][1], m[3][3], prec);
    mpfr_mul(t1, m[0][2], t0, MPFR_RNDN);
    mpfr_add(result, result, t1, MPFR_RNDN);

    // j=3: -m[0][3] * det3(cols 0,1,2 of rows 1,2,3)
    mpfr_det3(t0,
        m[1][0], m[1][1], m[1][2],
        m[2][0], m[2][1], m[2][2],
        m[3][0], m[3][1], m[3][2], prec);
    mpfr_mul(t1, m[0][3], t0, MPFR_RNDN);
    mpfr_sub(result, result, t1, MPFR_RNDN);

    mpfr_clear(t0);
    mpfr_clear(t1);
}

struct cdt_sign_cross_dprime_mpfr {
    static constexpr bool stateful = false;
    static constexpr bool updates  = false;

    template <typename Real>
    static int apply(
        Real aAx,  Real aAy,  Real aAz,
        Real bAx,  Real bAy,  Real bAz,
        Real cAx,  Real cAy,  Real cAz,
        Real LAx,  Real LAy,  Real LAz,
        Real RvAx, Real RvAy, Real RvAz,
        Real aBx,  Real aBy,  Real aBz,
        Real bBx,  Real bBy,  Real bBz,
        Real cBx,  Real cBy,  Real cBz,
        Real LBx,  Real LBy,  Real LBz,
        Real RvBx, Real RvBy, Real RvBz,
        Real pax,  Real pay,  Real paz,
        Real pbx,  Real pby,  Real pbz,
        Real pcx,  Real pcy,  Real pcz)
    {
        static constexpr mpfr_prec_t prec = 640;

        mpfr_t nnx, nny, nnz;
        mpfr_t m[4][4];
        mpfr_t t0, t1;
        mpfr_t D0A, DpA, D0B, DpB;
        mpfr_t prod1, prod2, diff;

        // --- init ---
        mpfr_init2(nnx, prec); mpfr_init2(nny, prec); mpfr_init2(nnz, prec);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                mpfr_init2(m[i][j], prec);
        mpfr_init2(t0,    prec); mpfr_init2(t1,    prec);
        mpfr_init2(D0A,   prec); mpfr_init2(DpA,   prec);
        mpfr_init2(D0B,   prec); mpfr_init2(DpB,   prec);
        mpfr_init2(prod1, prec); mpfr_init2(prod2, prec); mpfr_init2(diff, prec);

        // --- face normal nn = -(pb-pa) x (pc-pa), same sign as cdt_dprime_impl ---
        // nnx = (pbz-paz)*(pcy-pay) - (pby-pay)*(pcz-paz)
        mpfr_set_d(t0, (double)pbz, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)paz, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcy, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)pay, MPFR_RNDN);
        mpfr_mul(nnx, t0, t1, MPFR_RNDN);
        mpfr_set_d(t0, (double)pby, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)pay, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcz, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)paz, MPFR_RNDN);
        mpfr_mul(t0, t0, t1, MPFR_RNDN);
        mpfr_sub(nnx, nnx, t0, MPFR_RNDN);

        // nny = (pbx-pax)*(pcz-paz) - (pbz-paz)*(pcx-pax)
        mpfr_set_d(t0, (double)pbx, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)pax, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcz, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)paz, MPFR_RNDN);
        mpfr_mul(nny, t0, t1, MPFR_RNDN);
        mpfr_set_d(t0, (double)pbz, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)paz, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcx, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)pax, MPFR_RNDN);
        mpfr_mul(t0, t0, t1, MPFR_RNDN);
        mpfr_sub(nny, nny, t0, MPFR_RNDN);

        // nnz = (pby-pay)*(pcx-pax) - (pbx-pax)*(pcy-pay)
        mpfr_set_d(t0, (double)pby, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)pay, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcx, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)pax, MPFR_RNDN);
        mpfr_mul(nnz, t0, t1, MPFR_RNDN);
        mpfr_set_d(t0, (double)pbx, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)pax, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcy, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)pay, MPFR_RNDN);
        mpfr_mul(t0, t0, t1, MPFR_RNDN);
        mpfr_sub(nnz, nnz, t0, MPFR_RNDN);

        // --- helper: set m[row][0..2] = pt - rv ---
        auto set_row = [&](int row,
                           double px, double py, double pz,
                           double rvx, double rvy, double rvz) {
            mpfr_set_d(m[row][0], px,  MPFR_RNDN); mpfr_sub_d(m[row][0], m[row][0], rvx, MPFR_RNDN);
            mpfr_set_d(m[row][1], py,  MPFR_RNDN); mpfr_sub_d(m[row][1], m[row][1], rvy, MPFR_RNDN);
            mpfr_set_d(m[row][2], pz,  MPFR_RNDN); mpfr_sub_d(m[row][2], m[row][2], rvz, MPFR_RNDN);
        };

        // --- helper: fill 4th column as |d|^2 and compute det into dst ---
        auto det_D0 = [&](mpfr_t dst) {
            for (int i = 0; i < 4; i++) {
                mpfr_mul(t0, m[i][0], m[i][0], MPFR_RNDN);
                mpfr_mul(t1, m[i][1], m[i][1], MPFR_RNDN);
                mpfr_add(t0, t0, t1, MPFR_RNDN);
                mpfr_mul(t1, m[i][2], m[i][2], MPFR_RNDN);
                mpfr_add(m[i][3], t0, t1, MPFR_RNDN);
            }
            mpfr_det4(dst, m, prec);
        };

        // --- helper: fill 4th column as dot(d, nn) and compute det into dst ---
        auto det_Dp = [&](mpfr_t dst) {
            for (int i = 0; i < 4; i++) {
                mpfr_mul(t0, m[i][0], nnx, MPFR_RNDN);
                mpfr_mul(t1, m[i][1], nny, MPFR_RNDN);
                mpfr_add(t0, t0, t1, MPFR_RNDN);
                mpfr_mul(t1, m[i][2], nnz, MPFR_RNDN);
                mpfr_add(m[i][3], t0, t1, MPFR_RNDN);
            }
            mpfr_det4(dst, m, prec);
        };

        // --- event A: rows aA, bA, cA, LA relative to RvA ---
        set_row(0, (double)aAx, (double)aAy, (double)aAz, (double)RvAx, (double)RvAy, (double)RvAz);
        set_row(1, (double)bAx, (double)bAy, (double)bAz, (double)RvAx, (double)RvAy, (double)RvAz);
        set_row(2, (double)cAx, (double)cAy, (double)cAz, (double)RvAx, (double)RvAy, (double)RvAz);
        set_row(3, (double)LAx, (double)LAy, (double)LAz, (double)RvAx, (double)RvAy, (double)RvAz);
        det_D0(D0A);
        det_Dp(DpA);

        // --- event B: rows aB, bB, cB, LB relative to RvB ---
        set_row(0, (double)aBx, (double)aBy, (double)aBz, (double)RvBx, (double)RvBy, (double)RvBz);
        set_row(1, (double)bBx, (double)bBy, (double)bBz, (double)RvBx, (double)RvBy, (double)RvBz);
        set_row(2, (double)cBx, (double)cBy, (double)cBz, (double)RvBx, (double)RvBy, (double)RvBz);
        set_row(3, (double)LBx, (double)LBy, (double)LBz, (double)RvBx, (double)RvBy, (double)RvBz);
        det_D0(D0B);
        det_Dp(DpB);

        // --- P = D0_B * D'_A - D0_A * D'_B ---
        mpfr_mul(prod1, D0B, DpA, MPFR_RNDN);
        mpfr_mul(prod2, D0A, DpB, MPFR_RNDN);
        mpfr_sub(diff, prod1, prod2, MPFR_RNDN);

        int s = mpfr_sgn(diff);

        // --- clear ---
        mpfr_clear(nnx); mpfr_clear(nny); mpfr_clear(nnz);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                mpfr_clear(m[i][j]);
        mpfr_clear(t0);    mpfr_clear(t1);
        mpfr_clear(D0A);   mpfr_clear(DpA);
        mpfr_clear(D0B);   mpfr_clear(DpB);
        mpfr_clear(prod1); mpfr_clear(prod2); mpfr_clear(diff);

        return s;
    }
};

// ---------------------------------------------------------------------------
// cdt_sign_cross_dprime
// sign of D0_B*D'_A - D0_A*D'_B, used to compare flip times τ_A < τ_B.
// D0_X = insphere det for event X; D'_X = dprime det for event X.
// Inputs: aA(_1-_3), bA(_4-_6), cA(_7-_9), LA(_10-_12), RvA(_13-_15),
//         aB(_16-_18), bB(_19-_21), cB(_22-_24), LB(_25-_27), RvB(_28-_30),
//         pa(_31-_33), pb(_34-_36), pc(_37-_39)
// ---------------------------------------------------------------------------
namespace cdt_sign_cross_dprime_impl {
    constexpr auto aAx  = grp::_1;
    constexpr auto aAy  = grp::_2;
    constexpr auto aAz  = grp::_3;
    constexpr auto bAx  = grp::_4;
    constexpr auto bAy  = grp::_5;
    constexpr auto bAz  = grp::_6;
    constexpr auto cAx  = grp::_7;
    constexpr auto cAy  = grp::_8;
    constexpr auto cAz  = grp::_9;
    constexpr auto LAx  = grp::_10;
    constexpr auto LAy  = grp::_11;
    constexpr auto LAz  = grp::_12;
    constexpr auto RvAx = grp::_13;
    constexpr auto RvAy = grp::_14;
    constexpr auto RvAz = grp::_15;

    constexpr auto aBx  = grp::_16;
    constexpr auto aBy  = grp::_17;
    constexpr auto aBz  = grp::_18;
    constexpr auto bBx  = grp::_19;
    constexpr auto bBy  = grp::_20;
    constexpr auto bBz  = grp::_21;
    constexpr auto cBx  = grp::_22;
    constexpr auto cBy  = grp::_23;
    constexpr auto cBz  = grp::_24;
    constexpr auto LBx  = grp::_25;
    constexpr auto LBy  = grp::_26;
    constexpr auto LBz  = grp::_27;
    constexpr auto RvBx = grp::_28;
    constexpr auto RvBy = grp::_29;
    constexpr auto RvBz = grp::_30;

    constexpr auto pax  = grp::_31;
    constexpr auto pay  = grp::_32;
    constexpr auto paz  = grp::_33;
    constexpr auto pbx  = grp::_34;
    constexpr auto pby  = grp::_35;
    constexpr auto pbz  = grp::_36;
    constexpr auto pcx  = grp::_37;
    constexpr auto pcy  = grp::_38;
    constexpr auto pcz  = grp::_39;

    // -- rows for event A: p - RvA --
    constexpr auto daAx = aAx - RvAx;
    constexpr auto daAy = aAy - RvAy;
    constexpr auto daAz = aAz - RvAz;
    constexpr auto dbAx = bAx - RvAx;
    constexpr auto dbAy = bAy - RvAy;
    constexpr auto dbAz = bAz - RvAz;
    constexpr auto dcAx = cAx - RvAx;
    constexpr auto dcAy = cAy - RvAy;
    constexpr auto dcAz = cAz - RvAz;
    constexpr auto dLAx = LAx - RvAx;
    constexpr auto dLAy = LAy - RvAy;
    constexpr auto dLAz = LAz - RvAz;

    // -- rows for event B: p - RvB --
    constexpr auto daBx = aBx - RvBx;
    constexpr auto daBy = aBy - RvBy;
    constexpr auto daBz = aBz - RvBz;
    constexpr auto dbBx = bBx - RvBx;
    constexpr auto dbBy = bBy - RvBy;
    constexpr auto dbBz = bBz - RvBz;
    constexpr auto dcBx = cBx - RvBx;
    constexpr auto dcBy = cBy - RvBy;
    constexpr auto dcBz = cBz - RvBz;
    constexpr auto dLBx = LBx - RvBx;
    constexpr auto dLBy = LBy - RvBy;
    constexpr auto dLBz = LBz - RvBz;

    // -- D0_A: 4x4 insphere det for event A (4th column = squared norm) --
    constexpr auto daAw = daAx*daAx + daAy*daAy + daAz*daAz;
    constexpr auto dbAw = dbAx*dbAx + dbAy*dbAy + dbAz*dbAz;
    constexpr auto dcAw = dcAx*dcAx + dcAy*dcAy + dcAz*dcAz;
    constexpr auto dLAw = dLAx*dLAx + dLAy*dLAy + dLAz*dLAz;

    using D0_A_t = grp::det<
        decltype(daAx), decltype(daAy), decltype(daAz), decltype(daAw),
        decltype(dbAx), decltype(dbAy), decltype(dbAz), decltype(dbAw),
        decltype(dcAx), decltype(dcAy), decltype(dcAz), decltype(dcAw),
        decltype(dLAx), decltype(dLAy), decltype(dLAz), decltype(dLAw)
    >;
    constexpr auto D0_A = D0_A_t{};

    // -- D0_B: 4x4 insphere det for event B --
    constexpr auto daBw = daBx*daBx + daBy*daBy + daBz*daBz;
    constexpr auto dbBw = dbBx*dbBx + dbBy*dbBy + dbBz*dbBz;
    constexpr auto dcBw = dcBx*dcBx + dcBy*dcBy + dcBz*dcBz;
    constexpr auto dLBw = dLBx*dLBx + dLBy*dLBy + dLBz*dLBz;

    using D0_B_t = grp::det<
        decltype(daBx), decltype(daBy), decltype(daBz), decltype(daBw),
        decltype(dbBx), decltype(dbBy), decltype(dbBz), decltype(dbBw),
        decltype(dcBx), decltype(dcBy), decltype(dcBz), decltype(dcBw),
        decltype(dLBx), decltype(dLBy), decltype(dLBz), decltype(dLBw)
    >;
    constexpr auto D0_B = D0_B_t{};

    // -- face normal n = (pb-pa) x (pc-pa), stored negated as in cdt_dprime --
    constexpr auto dpbpax = pbx - pax;
    constexpr auto dpbpay = pby - pay;
    constexpr auto dpbpaz = pbz - paz;
    constexpr auto dpcpax = pcx - pax;
    constexpr auto dpcpay = pcy - pay;
    constexpr auto dpcpaz = pcz - paz;

    constexpr auto nnx = dpbpaz * dpcpay - dpbpay * dpcpaz;
    constexpr auto nny = dpbpax * dpcpaz - dpbpaz * dpcpax;
    constexpr auto nnz = dpbpay * dpcpax - dpbpax * dpcpay;

    // -- D'_A: 4th column = dot(p - RvA, -n) = dot(RvA - p, n) --
    constexpr auto col4_aA = daAx * nnx + daAy * nny + daAz * nnz;
    constexpr auto col4_bA = dbAx * nnx + dbAy * nny + dbAz * nnz;
    constexpr auto col4_cA = dcAx * nnx + dcAy * nny + dcAz * nnz;
    constexpr auto col4_LA = dLAx * nnx + dLAy * nny + dLAz * nnz;

    using Dp_A_t = grp::det<
        decltype(daAx), decltype(daAy), decltype(daAz), decltype(col4_aA),
        decltype(dbAx), decltype(dbAy), decltype(dbAz), decltype(col4_bA),
        decltype(dcAx), decltype(dcAy), decltype(dcAz), decltype(col4_cA),
        decltype(dLAx), decltype(dLAy), decltype(dLAz), decltype(col4_LA)
    >;
    constexpr auto Dp_A = Dp_A_t{};

    // -- D'_B: 4th column = dot(p - RvB, -n) --
    constexpr auto col4_aB = daBx * nnx + daBy * nny + daBz * nnz;
    constexpr auto col4_bB = dbBx * nnx + dbBy * nny + dbBz * nnz;
    constexpr auto col4_cB = dcBx * nnx + dcBy * nny + dcBz * nnz;
    constexpr auto col4_LB = dLBx * nnx + dLBy * nny + dLBz * nnz;

    using Dp_B_t = grp::det<
        decltype(daBx), decltype(daBy), decltype(daBz), decltype(col4_aB),
        decltype(dbBx), decltype(dbBy), decltype(dbBz), decltype(col4_bB),
        decltype(dcBx), decltype(dcBy), decltype(dcBz), decltype(col4_cB),
        decltype(dLBx), decltype(dLBy), decltype(dLBz), decltype(col4_LB)
    >;
    constexpr auto Dp_B = Dp_B_t{};

    // -- P = D0_B * D'_A - D0_A * D'_B, as a 2x2 det --
    using expr_t = grp::det<decltype(D0_B), decltype(D0_A),
                            decltype(Dp_B), decltype(Dp_A)>;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using pred        = grp::staged_predicate<semi_static, cdt_sign_cross_dprime_mpfr>;
}

extern "C" int cdt_sign_cross_dprime(
    double aAx,  double aAy,  double aAz,
    double bAx,  double bAy,  double bAz,
    double cAx,  double cAy,  double cAz,
    double LAx,  double LAy,  double LAz,
    double RvAx, double RvAy, double RvAz,
    double aBx,  double aBy,  double aBz,
    double bBx,  double bBy,  double bBz,
    double cBx,  double cBy,  double cBz,
    double LBx,  double LBy,  double LBz,
    double RvBx, double RvBy, double RvBz,
    double pax,  double pay,  double paz,
    double pbx,  double pby,  double pbz,
    double pcx,  double pcy,  double pcz)
{
    return cdt_sign_cross_dprime_impl::pred{}.apply(
        aAx,  aAy,  aAz,
        bAx,  bAy,  bAz,
        cAx,  cAy,  cAz,
        LAx,  LAy,  LAz,
        RvAx, RvAy, RvAz,
        aBx,  aBy,  aBz,
        bBx,  bBy,  bBz,
        cBx,  cBy,  cBz,
        LBx,  LBy,  LBz,
        RvBx, RvBy, RvBz,
        pax,  pay,  paz,
        pbx,  pby,  pbz,
        pcx,  pcy,  pcz);
}

// ---------------------------------------------------------------------------
// cdt_sign_cross_dprime_weighted
// Same predicate as cdt_sign_cross_dprime but with additive weight offsets
// kA[5] / kB[5] on each vertex's insphere 4th column (used for SoS cascade).
// D' is unaffected by weights.
// ---------------------------------------------------------------------------

struct cdt_sign_cross_dprime_weighted_mpfr {
    static constexpr bool stateful = false;
    static constexpr bool updates  = false;

    template <typename Real>
    static int apply(
        Real aAx,  Real aAy,  Real aAz,
        Real bAx,  Real bAy,  Real bAz,
        Real cAx,  Real cAy,  Real cAz,
        Real LAx,  Real LAy,  Real LAz,
        Real RvAx, Real RvAy, Real RvAz,
        Real aBx,  Real aBy,  Real aBz,
        Real bBx,  Real bBy,  Real bBz,
        Real cBx,  Real cBy,  Real cBz,
        Real LBx,  Real LBy,  Real LBz,
        Real RvBx, Real RvBy, Real RvBz,
        Real pax,  Real pay,  Real paz,
        Real pbx,  Real pby,  Real pbz,
        Real pcx,  Real pcy,  Real pcz,
        Real k_aA, Real k_bA, Real k_cA, Real k_LA, Real k_RvA,
        Real k_aB, Real k_bB, Real k_cB, Real k_LB, Real k_RvB)
    {
        static constexpr mpfr_prec_t prec = 704;

        mpfr_t nnx, nny, nnz;
        mpfr_t m[4][4];
        mpfr_t t0, t1;
        mpfr_t D0A, DpA, D0B, DpB;
        mpfr_t prod1, prod2, diff;

        mpfr_init2(nnx, prec); mpfr_init2(nny, prec); mpfr_init2(nnz, prec);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                mpfr_init2(m[i][j], prec);
        mpfr_init2(t0,    prec); mpfr_init2(t1,    prec);
        mpfr_init2(D0A,   prec); mpfr_init2(DpA,   prec);
        mpfr_init2(D0B,   prec); mpfr_init2(DpB,   prec);
        mpfr_init2(prod1, prec); mpfr_init2(prod2, prec); mpfr_init2(diff, prec);

        // face normal nn = -(pb-pa) x (pc-pa)
        mpfr_set_d(t0, (double)pbz, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)paz, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcy, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)pay, MPFR_RNDN);
        mpfr_mul(nnx, t0, t1, MPFR_RNDN);
        mpfr_set_d(t0, (double)pby, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)pay, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcz, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)paz, MPFR_RNDN);
        mpfr_mul(t0, t0, t1, MPFR_RNDN);
        mpfr_sub(nnx, nnx, t0, MPFR_RNDN);

        mpfr_set_d(t0, (double)pbx, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)pax, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcz, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)paz, MPFR_RNDN);
        mpfr_mul(nny, t0, t1, MPFR_RNDN);
        mpfr_set_d(t0, (double)pbz, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)paz, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcx, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)pax, MPFR_RNDN);
        mpfr_mul(t0, t0, t1, MPFR_RNDN);
        mpfr_sub(nny, nny, t0, MPFR_RNDN);

        mpfr_set_d(t0, (double)pby, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)pay, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcx, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)pax, MPFR_RNDN);
        mpfr_mul(nnz, t0, t1, MPFR_RNDN);
        mpfr_set_d(t0, (double)pbx, MPFR_RNDN); mpfr_sub_d(t0, t0, (double)pax, MPFR_RNDN);
        mpfr_set_d(t1, (double)pcy, MPFR_RNDN); mpfr_sub_d(t1, t1, (double)pay, MPFR_RNDN);
        mpfr_mul(t0, t0, t1, MPFR_RNDN);
        mpfr_sub(nnz, nnz, t0, MPFR_RNDN);

        auto set_row = [&](int row,
                           double px, double py, double pz,
                           double rvx, double rvy, double rvz) {
            mpfr_set_d(m[row][0], px,  MPFR_RNDN); mpfr_sub_d(m[row][0], m[row][0], rvx, MPFR_RNDN);
            mpfr_set_d(m[row][1], py,  MPFR_RNDN); mpfr_sub_d(m[row][1], m[row][1], rvy, MPFR_RNDN);
            mpfr_set_d(m[row][2], pz,  MPFR_RNDN); mpfr_sub_d(m[row][2], m[row][2], rvz, MPFR_RNDN);
        };

        // D0 with weight offsets: 4th col = |d|^2 + (k_i - k_Rv)
        auto det_D0_w = [&](mpfr_t dst, double k[4]) {
            for (int i = 0; i < 4; i++) {
                mpfr_mul(t0, m[i][0], m[i][0], MPFR_RNDN);
                mpfr_mul(t1, m[i][1], m[i][1], MPFR_RNDN);
                mpfr_add(t0, t0, t1, MPFR_RNDN);
                mpfr_mul(t1, m[i][2], m[i][2], MPFR_RNDN);
                mpfr_add(m[i][3], t0, t1, MPFR_RNDN);
                if (k[i] != 0.0) mpfr_add_d(m[i][3], m[i][3], k[i], MPFR_RNDN);
            }
            mpfr_det4(dst, m, prec);
        };

        auto det_Dp = [&](mpfr_t dst) {
            for (int i = 0; i < 4; i++) {
                mpfr_mul(t0, m[i][0], nnx, MPFR_RNDN);
                mpfr_mul(t1, m[i][1], nny, MPFR_RNDN);
                mpfr_add(t0, t0, t1, MPFR_RNDN);
                mpfr_mul(t1, m[i][2], nnz, MPFR_RNDN);
                mpfr_add(m[i][3], t0, t1, MPFR_RNDN);
            }
            mpfr_det4(dst, m, prec);
        };

        // event A
        set_row(0, (double)aAx, (double)aAy, (double)aAz, (double)RvAx, (double)RvAy, (double)RvAz);
        set_row(1, (double)bAx, (double)bAy, (double)bAz, (double)RvAx, (double)RvAy, (double)RvAz);
        set_row(2, (double)cAx, (double)cAy, (double)cAz, (double)RvAx, (double)RvAy, (double)RvAz);
        set_row(3, (double)LAx, (double)LAy, (double)LAz, (double)RvAx, (double)RvAy, (double)RvAz);
        {
            double kA[4] = { (double)k_aA - (double)k_RvA,
                             (double)k_bA - (double)k_RvA,
                             (double)k_cA - (double)k_RvA,
                             (double)k_LA - (double)k_RvA };
            det_D0_w(D0A, kA);
        }
        det_Dp(DpA);

        // event B
        set_row(0, (double)aBx, (double)aBy, (double)aBz, (double)RvBx, (double)RvBy, (double)RvBz);
        set_row(1, (double)bBx, (double)bBy, (double)bBz, (double)RvBx, (double)RvBy, (double)RvBz);
        set_row(2, (double)cBx, (double)cBy, (double)cBz, (double)RvBx, (double)RvBy, (double)RvBz);
        set_row(3, (double)LBx, (double)LBy, (double)LBz, (double)RvBx, (double)RvBy, (double)RvBz);
        {
            double kB[4] = { (double)k_aB - (double)k_RvB,
                             (double)k_bB - (double)k_RvB,
                             (double)k_cB - (double)k_RvB,
                             (double)k_LB - (double)k_RvB };
            det_D0_w(D0B, kB);
        }
        det_Dp(DpB);

        mpfr_mul(prod1, D0B, DpA, MPFR_RNDN);
        mpfr_mul(prod2, D0A, DpB, MPFR_RNDN);
        mpfr_sub(diff, prod1, prod2, MPFR_RNDN);
        int s = mpfr_sgn(diff);

        mpfr_clear(nnx); mpfr_clear(nny); mpfr_clear(nnz);
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                mpfr_clear(m[i][j]);
        mpfr_clear(t0);    mpfr_clear(t1);
        mpfr_clear(D0A);   mpfr_clear(DpA);
        mpfr_clear(D0B);   mpfr_clear(DpB);
        mpfr_clear(prod1); mpfr_clear(prod2); mpfr_clear(diff);
        return s;
    }
};

namespace cdt_sign_cross_dprime_weighted_impl {
    // coordinate inputs (_1–_39) same mapping as cdt_sign_cross_dprime_impl
    constexpr auto aAx  = grp::_1;  constexpr auto aAy  = grp::_2;  constexpr auto aAz  = grp::_3;
    constexpr auto bAx  = grp::_4;  constexpr auto bAy  = grp::_5;  constexpr auto bAz  = grp::_6;
    constexpr auto cAx  = grp::_7;  constexpr auto cAy  = grp::_8;  constexpr auto cAz  = grp::_9;
    constexpr auto LAx  = grp::_10; constexpr auto LAy  = grp::_11; constexpr auto LAz  = grp::_12;
    constexpr auto RvAx = grp::_13; constexpr auto RvAy = grp::_14; constexpr auto RvAz = grp::_15;
    constexpr auto aBx  = grp::_16; constexpr auto aBy  = grp::_17; constexpr auto aBz  = grp::_18;
    constexpr auto bBx  = grp::_19; constexpr auto bBy  = grp::_20; constexpr auto bBz  = grp::_21;
    constexpr auto cBx  = grp::_22; constexpr auto cBy  = grp::_23; constexpr auto cBz  = grp::_24;
    constexpr auto LBx  = grp::_25; constexpr auto LBy  = grp::_26; constexpr auto LBz  = grp::_27;
    constexpr auto RvBx = grp::_28; constexpr auto RvBy = grp::_29; constexpr auto RvBz = grp::_30;
    constexpr auto pax  = grp::_31; constexpr auto pay  = grp::_32; constexpr auto paz  = grp::_33;
    constexpr auto pbx  = grp::_34; constexpr auto pby  = grp::_35; constexpr auto pbz  = grp::_36;
    constexpr auto pcx  = grp::_37; constexpr auto pcy  = grp::_38; constexpr auto pcz  = grp::_39;

    // weight inputs: kA[0..4] = _40.._44, kB[0..4] = _45.._49
    constexpr auto k_aA  = grp::_40; constexpr auto k_bA  = grp::_41;
    constexpr auto k_cA  = grp::_42; constexpr auto k_LA  = grp::_43;
    constexpr auto k_RvA = grp::_44;
    constexpr auto k_aB  = grp::_45; constexpr auto k_bB  = grp::_46;
    constexpr auto k_cB  = grp::_47; constexpr auto k_LB  = grp::_48;
    constexpr auto k_RvB = grp::_49;

    // difference vectors (same as unweighted)
    constexpr auto daAx = aAx - RvAx; constexpr auto daAy = aAy - RvAy; constexpr auto daAz = aAz - RvAz;
    constexpr auto dbAx = bAx - RvAx; constexpr auto dbAy = bAy - RvAy; constexpr auto dbAz = bAz - RvAz;
    constexpr auto dcAx = cAx - RvAx; constexpr auto dcAy = cAy - RvAy; constexpr auto dcAz = cAz - RvAz;
    constexpr auto dLAx = LAx - RvAx; constexpr auto dLAy = LAy - RvAy; constexpr auto dLAz = LAz - RvAz;
    constexpr auto daBx = aBx - RvBx; constexpr auto daBy = aBy - RvBy; constexpr auto daBz = aBz - RvBz;
    constexpr auto dbBx = bBx - RvBx; constexpr auto dbBy = bBy - RvBy; constexpr auto dbBz = bBz - RvBz;
    constexpr auto dcBx = cBx - RvBx; constexpr auto dcBy = cBy - RvBy; constexpr auto dcBz = cBz - RvBz;
    constexpr auto dLBx = LBx - RvBx; constexpr auto dLBy = LBy - RvBy; constexpr auto dLBz = LBz - RvBz;

    // D0_A: 4th column = |d|^2 + (k_i - k_RvA)
    constexpr auto daAw = daAx*daAx + daAy*daAy + daAz*daAz + k_aA - k_RvA;
    constexpr auto dbAw = dbAx*dbAx + dbAy*dbAy + dbAz*dbAz + k_bA - k_RvA;
    constexpr auto dcAw = dcAx*dcAx + dcAy*dcAy + dcAz*dcAz + k_cA - k_RvA;
    constexpr auto dLAw = dLAx*dLAx + dLAy*dLAy + dLAz*dLAz + k_LA - k_RvA;

    using D0_A_t = grp::det<
        decltype(daAx), decltype(daAy), decltype(daAz), decltype(daAw),
        decltype(dbAx), decltype(dbAy), decltype(dbAz), decltype(dbAw),
        decltype(dcAx), decltype(dcAy), decltype(dcAz), decltype(dcAw),
        decltype(dLAx), decltype(dLAy), decltype(dLAz), decltype(dLAw)
    >;
    constexpr auto D0_A = D0_A_t{};

    // D0_B: 4th column = |d|^2 + (k_i - k_RvB)
    constexpr auto daBw = daBx*daBx + daBy*daBy + daBz*daBz + k_aB - k_RvB;
    constexpr auto dbBw = dbBx*dbBx + dbBy*dbBy + dbBz*dbBz + k_bB - k_RvB;
    constexpr auto dcBw = dcBx*dcBx + dcBy*dcBy + dcBz*dcBz + k_cB - k_RvB;
    constexpr auto dLBw = dLBx*dLBx + dLBy*dLBy + dLBz*dLBz + k_LB - k_RvB;

    using D0_B_t = grp::det<
        decltype(daBx), decltype(daBy), decltype(daBz), decltype(daBw),
        decltype(dbBx), decltype(dbBy), decltype(dbBz), decltype(dbBw),
        decltype(dcBx), decltype(dcBy), decltype(dcBz), decltype(dcBw),
        decltype(dLBx), decltype(dLBy), decltype(dLBz), decltype(dLBw)
    >;
    constexpr auto D0_B = D0_B_t{};

    // face normal and D' (unchanged from unweighted)
    constexpr auto dpbpax = pbx - pax; constexpr auto dpbpay = pby - pay; constexpr auto dpbpaz = pbz - paz;
    constexpr auto dpcpax = pcx - pax; constexpr auto dpcpay = pcy - pay; constexpr auto dpcpaz = pcz - paz;
    constexpr auto nnx = dpbpaz * dpcpay - dpbpay * dpcpaz;
    constexpr auto nny = dpbpax * dpcpaz - dpbpaz * dpcpax;
    constexpr auto nnz = dpbpay * dpcpax - dpbpax * dpcpay;

    constexpr auto col4_aA = daAx * nnx + daAy * nny + daAz * nnz;
    constexpr auto col4_bA = dbAx * nnx + dbAy * nny + dbAz * nnz;
    constexpr auto col4_cA = dcAx * nnx + dcAy * nny + dcAz * nnz;
    constexpr auto col4_LA = dLAx * nnx + dLAy * nny + dLAz * nnz;
    using Dp_A_t = grp::det<
        decltype(daAx), decltype(daAy), decltype(daAz), decltype(col4_aA),
        decltype(dbAx), decltype(dbAy), decltype(dbAz), decltype(col4_bA),
        decltype(dcAx), decltype(dcAy), decltype(dcAz), decltype(col4_cA),
        decltype(dLAx), decltype(dLAy), decltype(dLAz), decltype(col4_LA)
    >;
    constexpr auto Dp_A = Dp_A_t{};

    constexpr auto col4_aB = daBx * nnx + daBy * nny + daBz * nnz;
    constexpr auto col4_bB = dbBx * nnx + dbBy * nny + dbBz * nnz;
    constexpr auto col4_cB = dcBx * nnx + dcBy * nny + dcBz * nnz;
    constexpr auto col4_LB = dLBx * nnx + dLBy * nny + dLBz * nnz;
    using Dp_B_t = grp::det<
        decltype(daBx), decltype(daBy), decltype(daBz), decltype(col4_aB),
        decltype(dbBx), decltype(dbBy), decltype(dbBz), decltype(col4_bB),
        decltype(dcBx), decltype(dcBy), decltype(dcBz), decltype(col4_cB),
        decltype(dLBx), decltype(dLBy), decltype(dLBz), decltype(col4_LB)
    >;
    constexpr auto Dp_B = Dp_B_t{};

    using expr_t = grp::det<decltype(D0_B), decltype(D0_A),
                            decltype(Dp_B), decltype(Dp_A)>;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using pred        = grp::staged_predicate<semi_static, cdt_sign_cross_dprime_weighted_mpfr>;
}

extern "C" int cdt_sign_cross_dprime_weighted(
    double aAx,  double aAy,  double aAz,
    double bAx,  double bAy,  double bAz,
    double cAx,  double cAy,  double cAz,
    double LAx,  double LAy,  double LAz,
    double RvAx, double RvAy, double RvAz,
    double aBx,  double aBy,  double aBz,
    double bBx,  double bBy,  double bBz,
    double cBx,  double cBy,  double cBz,
    double LBx,  double LBy,  double LBz,
    double RvBx, double RvBy, double RvBz,
    double pax,  double pay,  double paz,
    double pbx,  double pby,  double pbz,
    double pcx,  double pcy,  double pcz,
    const double kA[5], const double kB[5])
{
    return cdt_sign_cross_dprime_weighted_impl::pred{}.apply(
        aAx,  aAy,  aAz,
        bAx,  bAy,  bAz,
        cAx,  cAy,  cAz,
        LAx,  LAy,  LAz,
        RvAx, RvAy, RvAz,
        aBx,  aBy,  aBz,
        bBx,  bBy,  bBz,
        cBx,  cBy,  cBz,
        LBx,  LBy,  LBz,
        RvBx, RvBy, RvBz,
        pax,  pay,  paz,
        pbx,  pby,  pbz,
        pcx,  pcy,  pcz,
        kA[0], kA[1], kA[2], kA[3], kA[4],
        kB[0], kB[1], kB[2], kB[3], kB[4]);
}

#ifdef ROBUST_PREDICATES_PRINT_SIZE
template <std::size_t N>
struct [[deprecated("results_size -- see template argument")]] show_stage_d_size {};
using _size_lp_feasible_T0_T1_S = show_stage_d_size<lp_feasible_T0_T1_S_impl::exact::results_size>*;
using _size_lp_feasible_T1_T2_S = show_stage_d_size<lp_feasible_T1_T2_S_impl::exact::results_size>*;
using _size_lp_feasible_T0_T2_S = show_stage_d_size<lp_feasible_T0_T2_S_impl::exact::results_size>*;
using _size_lp_D_T0_S     = show_stage_d_size<lp_D_T0_S_impl::exact::results_size>*;
using _size_lp_D_T1_S     = show_stage_d_size<lp_D_T1_S_impl::exact::results_size>*;
using _size_lp_D_T2_S     = show_stage_d_size<lp_D_T2_S_impl::exact::results_size>*;
using _size_lp_det2       = show_stage_d_size<lp_det2_impl::exact::results_size>*;
using _size_lp_feasible_T0_S_S_f1 = show_stage_d_size<lp_feasible_T0_S_S_factor1_impl::exact::results_size>*;
using _size_lp_feasible_T0_S_S_f2 = show_stage_d_size<lp_feasible_T0_S_S_factor2_impl::exact::results_size>*;
using _size_lp_feasible_T1_S_S_f1 = show_stage_d_size<lp_feasible_T1_S_S_factor1_impl::exact::results_size>*;
using _size_lp_feasible_T1_S_S_f2 = show_stage_d_size<lp_feasible_T1_S_S_factor2_impl::exact::results_size>*;
using _size_lp_feasible_T2_S_S_f1 = show_stage_d_size<lp_feasible_T2_S_S_factor1_impl::exact::results_size>*;
using _size_lp_feasible_T2_S_S_f2 = show_stage_d_size<lp_feasible_T2_S_S_factor2_impl::exact::results_size>*;
using _size_lp_det3       = show_stage_d_size<lp_det3_impl::exact::results_size>*;
using _size_lp_feasible_S_S_T0_f2 = show_stage_d_size<lp_feasible_S_S_T0_factor2_impl::exact::results_size>*;
using _size_lp_feasible_S_S_T2_f2 = show_stage_d_size<lp_feasible_S_S_T2_factor2_impl::exact::results_size>*;
using _size_lp_feasible_S_S_T1_f2 = show_stage_d_size<lp_feasible_S_S_T1_factor2_impl::exact::results_size>*;
// lp_feasible_T1_S_T2 reuses lp_D_T1_S_impl and lp_feasible_T1_T2_S_impl -- no new size entry
// lp_feasible_T0_S_T2 reuses lp_D_T0_S_impl and lp_feasible_T0_T2_S_impl -- no new size entry
using _size_cdt_dprime = show_stage_d_size<cdt_dprime_impl::exact::results_size>*;
// cdt_sign_cross_dprime uses the MPFR exact stage -- no stage_d results_size
#endif
