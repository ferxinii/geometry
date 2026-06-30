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
#endif
