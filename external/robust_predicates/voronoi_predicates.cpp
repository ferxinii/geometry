#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"

namespace grp = boost::geometry::detail::generic_robust_predicates;

// Compile-time integer constant leaf, for scaling a sub-expression by an exact
// integer (e.g. the offset-column x2 in lp3_feasible_TSS_T_gen). Cheaper in
// stage_d than duplicating the sub-expression (x + x). `non_negative` lets the
// forward-error filter treat it as a known-sign constant.
template <int N>
struct int_const : grp::static_constant_interface<double> {
    static constexpr double value = static_cast<double>(N);
    static constexpr bool non_negative = (N >= 0);
};

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


// ===========================================================================
// Tier B: 3D LP-feasibility predicates (lp3_*).  See lp3_predicates.tex.
// Constraint types: T (tet face, normal n_T) and S (seed bisector s->t).
// A vertex is a triple of constraint boundaries; feasible = sign(Delta)*sign(Gamma)
// with Delta = det of the three normals and Gamma the 4x4 [normal | offset].
// When a T constraint helps define the vertex its Gamma row has a zero offset;
// we never build a literal 0 matrix entry (the filter cannot bound a zero leaf),
// but cofactor-expand Gamma along the offset column into a zero-free sum of
// c*(3x3 minor) terms.  See lp3_predicates.tex (Remark: zero-free Gamma).
// ===========================================================================

// ---------------------------------------------------------------------------
// lp3_D_TSS -- slope Delta of the triple (T, S_i, S_j)
// sign det[ n_T ; t1-s ; t2-s ] = sign( n_T . ((t1-s) x (t2-s)) )
//   n_T = (Q-A) x (R-A) for tet face (A,Q,R); the caller applies the fixed
//   per-face orientation sign so that the tet interior is the feasible side.
// Also the factor1 reused inside lp3_feasible_TSS_S and lp3_feasible_TSS_T.
// Inputs: A(_1-_3), Q(_4-_6), R(_7-_9), s(_10-_12), t1(_13-_15), t2(_16-_18)
// ---------------------------------------------------------------------------
namespace lp3_D_TSS_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Qx  = grp::_4;
    constexpr auto Qy  = grp::_5;
    constexpr auto Qz  = grp::_6;
    constexpr auto Rx  = grp::_7;
    constexpr auto Ry  = grp::_8;
    constexpr auto Rz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto t1x = grp::_13;
    constexpr auto t1y = grp::_14;
    constexpr auto t1z = grp::_15;
    constexpr auto t2x = grp::_16;
    constexpr auto t2y = grp::_17;
    constexpr auto t2z = grp::_18;

    // n_T = (Q-A) x (R-A)
    constexpr auto dQAx = Qx - Ax;
    constexpr auto dQAy = Qy - Ay;
    constexpr auto dQAz = Qz - Az;
    constexpr auto dRAx = Rx - Ax;
    constexpr auto dRAy = Ry - Ay;
    constexpr auto dRAz = Rz - Az;
    constexpr auto nTx = dQAy*dRAz - dQAz*dRAy;
    constexpr auto nTy = dQAz*dRAx - dQAx*dRAz;
    constexpr auto nTz = dQAx*dRAy - dQAy*dRAx;

    // t1-s, t2-s
    constexpr auto d1x = t1x - sx;
    constexpr auto d1y = t1y - sy;
    constexpr auto d1z = t1z - sz;
    constexpr auto d2x = t2x - sx;
    constexpr auto d2y = t2y - sy;
    constexpr auto d2z = t2z - sz;

    using expr_t = grp::det<
        decltype(nTx), decltype(nTy), decltype(nTz),
        decltype(d1x), decltype(d1y), decltype(d1z),
        decltype(d2x), decltype(d2y), decltype(d2z)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp3_D_TSS(double Ax, double Ay, double Az,
                         double Qx, double Qy, double Qz,
                         double Rx, double Ry, double Rz,
                         double sx, double sy, double sz,
                         double t1x, double t1y, double t1z,
                         double t2x, double t2y, double t2z)
{
    return lp3_D_TSS_impl::pred{}.apply(Ax,Ay,Az, Qx,Qy,Qz, Rx,Ry,Rz,
                                        sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z);
}


// ---------------------------------------------------------------------------
// lp3_feasible_TTS_S -- vertex (T,T,S) [= bisector s->t cutting tet edge AB]
// tested against a second bisector s->u.
// feasible = -sign(gA*hB - gB*hA) * sign(gA - gB)
//   gX = |X-s|^2 - |X-t|^2   (edge-defining bisector, at edge endpoints)
//   hX = |X-s|^2 - |X-u|^2   (query bisector)
// factor1 = sign(gA - gB) (the 1D-edge slope); factor2 = sign(gA*hB - gB*hA).
// Inputs: A(_1-_3), B(_4-_6), s(_7-_9), t(_10-_12), u(_13-_15)
// ---------------------------------------------------------------------------
namespace lp3_feasible_TTS_S_factor1_impl {
    constexpr auto Ax = grp::_1;
    constexpr auto Ay = grp::_2;
    constexpr auto Az = grp::_3;
    constexpr auto Bx = grp::_4;
    constexpr auto By = grp::_5;
    constexpr auto Bz = grp::_6;
    constexpr auto sx = grp::_7;
    constexpr auto sy = grp::_8;
    constexpr auto sz = grp::_9;
    constexpr auto tx = grp::_10;
    constexpr auto ty = grp::_11;
    constexpr auto tz = grp::_12;
    // _13 _14 _15 = u, unused

    constexpr auto dAsx = Ax - sx;  constexpr auto dAsy = Ay - sy;  constexpr auto dAsz = Az - sz;
    constexpr auto dAtx = Ax - tx;  constexpr auto dAty = Ay - ty;  constexpr auto dAtz = Az - tz;
    constexpr auto dBsx = Bx - sx;  constexpr auto dBsy = By - sy;  constexpr auto dBsz = Bz - sz;
    constexpr auto dBtx = Bx - tx;  constexpr auto dBty = By - ty;  constexpr auto dBtz = Bz - tz;

    constexpr auto gA = dAsx*dAsx + dAsy*dAsy + dAsz*dAsz
                      - dAtx*dAtx - dAty*dAty - dAtz*dAtz;
    constexpr auto gB = dBsx*dBsx + dBsy*dBsy + dBsz*dBsz
                      - dBtx*dBtx - dBty*dBty - dBtz*dBtz;
    constexpr auto expr = gA - gB;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace lp3_feasible_TTS_S_factor2_impl {
    constexpr auto Ax = grp::_1;
    constexpr auto Ay = grp::_2;
    constexpr auto Az = grp::_3;
    constexpr auto Bx = grp::_4;
    constexpr auto By = grp::_5;
    constexpr auto Bz = grp::_6;
    constexpr auto sx = grp::_7;
    constexpr auto sy = grp::_8;
    constexpr auto sz = grp::_9;
    constexpr auto tx = grp::_10;
    constexpr auto ty = grp::_11;
    constexpr auto tz = grp::_12;
    constexpr auto ux = grp::_13;
    constexpr auto uy = grp::_14;
    constexpr auto uz = grp::_15;

    constexpr auto dAsx = Ax - sx;  constexpr auto dAsy = Ay - sy;  constexpr auto dAsz = Az - sz;
    constexpr auto dAtx = Ax - tx;  constexpr auto dAty = Ay - ty;  constexpr auto dAtz = Az - tz;
    constexpr auto dAux = Ax - ux;  constexpr auto dAuy = Ay - uy;  constexpr auto dAuz = Az - uz;
    constexpr auto dBsx = Bx - sx;  constexpr auto dBsy = By - sy;  constexpr auto dBsz = Bz - sz;
    constexpr auto dBtx = Bx - tx;  constexpr auto dBty = By - ty;  constexpr auto dBtz = Bz - tz;
    constexpr auto dBux = Bx - ux;  constexpr auto dBuy = By - uy;  constexpr auto dBuz = Bz - uz;

    constexpr auto gA = dAsx*dAsx + dAsy*dAsy + dAsz*dAsz
                      - dAtx*dAtx - dAty*dAty - dAtz*dAtz;
    constexpr auto hA = dAsx*dAsx + dAsy*dAsy + dAsz*dAsz
                      - dAux*dAux - dAuy*dAuy - dAuz*dAuz;
    constexpr auto gB = dBsx*dBsx + dBsy*dBsy + dBsz*dBsz
                      - dBtx*dBtx - dBty*dBty - dBtz*dBtz;
    constexpr auto hB = dBsx*dBsx + dBsy*dBsy + dBsz*dBsz
                      - dBux*dBux - dBuy*dBuy - dBuz*dBuz;

    // det[ gA hA ; gB hB ] = gA*hB - hA*gB = gA*hB - gB*hA
    using det2_t = grp::det<decltype(gA), decltype(hA),
                            decltype(gB), decltype(hB)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp3_feasible_TTS_S(double Ax, double Ay, double Az,
                                  double Bx, double By, double Bz,
                                  double sx, double sy, double sz,
                                  double tx, double ty, double tz,
                                  double ux, double uy, double uz)
{
    return - lp3_feasible_TTS_S_factor2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz,
                                                           sx,sy,sz, tx,ty,tz, ux,uy,uz)
           * lp3_feasible_TTS_S_factor1_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz,
                                                           sx,sy,sz, tx,ty,tz, ux,uy,uz);
}

// ---------------------------------------------------------------------------
// lp3_slope_TTS -- sign(gA - gB), gX = |X-s|^2 - |X-t|^2, on tet edge AB.
// The 1D-edge slope of bisector s->t: >0 => bisector value decreases A->B
// (a LOwer bound on the edge parameter), <0 => an upper bound. (= factor1 of
// lp3_feasible_TTS_S.)  Used by the 1D-LP edge clip (category d).
// ---------------------------------------------------------------------------
extern "C" int lp3_slope_TTS(double Ax, double Ay, double Az,
                             double Bx, double By, double Bz,
                             double sx, double sy, double sz,
                             double tx, double ty, double tz)
{
    return lp3_feasible_TTS_S_factor1_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz,
                                                         sx,sy,sz, tx,ty,tz, tx,ty,tz);
}

// ---------------------------------------------------------------------------
// lp3_det_TTS -- sign(g1A*g2B - g1B*g2A), giX = |X-s|^2 - |X-ti|^2, edge AB.
// The 2x2 determinant that orders the two bisector crossings of s->t1 and
// s->t2 on the edge: sign(lam1 - lam2) = -lp3_det_TTS * slope1 * slope2, with
// lam_i = g_iA/(g_iA - g_iB).  (= factor2 of lp3_feasible_TTS_S.)
// ---------------------------------------------------------------------------
extern "C" int lp3_det_TTS(double Ax, double Ay, double Az,
                           double Bx, double By, double Bz,
                           double sx, double sy, double sz,
                           double t1x, double t1y, double t1z,
                           double t2x, double t2y, double t2z)
{
    return lp3_feasible_TTS_S_factor2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz,
                                                         sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z);
}


// ---------------------------------------------------------------------------
// lp3_feasible_TTS_T -- vertex (T,T,S) [bisector s->t cutting tet edge AB]
// tested against a query tet face Tc=(Pc,Qc,Rc).  Surface edge-ordering
// predicate (see lp3_predicates.tex, Prop. "TTS vs T, signed"): the clip omitted
// this because TTS-vs-T is automatic within the vertex's own tet, but ordering
// two vertices on a shared edge needs the genuine sign against another face.
// feasible = -sign(gA*wB - gB*wA) * sign(gA - gB)
//   gX = |X-s|^2 - |X-t|^2                          (edge bisector value at A,B)
//   wX = n_c . (X - Pc),  n_c = (Qc-Pc) x (Rc-Pc)   (query-face signed distance)
// factor1 = sign(gA - gB) reuses lp3_feasible_TTS_S_factor1; factor2 below is
// the TTS_S factor2 pattern with the query bisector h replaced by w.  Degree 5.
// NOTE: query is a T face with fixed n_c=(Qc-Pc)x(Rc-Pc); the caller multiplies
// by the face interior-orientation sign sigma_f (Face-orientation convention).
// Inputs: A(_1-_3), B(_4-_6), s(_7-_9), t(_10-_12),
//         Pc(_13-_15), Qc(_16-_18), Rc(_19-_21)
// ---------------------------------------------------------------------------
namespace lp3_feasible_TTS_T_factor2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Bx  = grp::_4;
    constexpr auto By  = grp::_5;
    constexpr auto Bz  = grp::_6;
    constexpr auto sx  = grp::_7;
    constexpr auto sy  = grp::_8;
    constexpr auto sz  = grp::_9;
    constexpr auto tx  = grp::_10;
    constexpr auto ty  = grp::_11;
    constexpr auto tz  = grp::_12;
    constexpr auto Pcx = grp::_13;
    constexpr auto Pcy = grp::_14;
    constexpr auto Pcz = grp::_15;
    constexpr auto Qcx = grp::_16;
    constexpr auto Qcy = grp::_17;
    constexpr auto Qcz = grp::_18;
    constexpr auto Rcx = grp::_19;
    constexpr auto Rcy = grp::_20;
    constexpr auto Rcz = grp::_21;

    // edge bisector values gA, gB = |X-s|^2 - |X-t|^2
    constexpr auto dAsx = Ax - sx;  constexpr auto dAsy = Ay - sy;  constexpr auto dAsz = Az - sz;
    constexpr auto dAtx = Ax - tx;  constexpr auto dAty = Ay - ty;  constexpr auto dAtz = Az - tz;
    constexpr auto dBsx = Bx - sx;  constexpr auto dBsy = By - sy;  constexpr auto dBsz = Bz - sz;
    constexpr auto dBtx = Bx - tx;  constexpr auto dBty = By - ty;  constexpr auto dBtz = Bz - tz;

    constexpr auto gA = dAsx*dAsx + dAsy*dAsy + dAsz*dAsz
                      - dAtx*dAtx - dAty*dAty - dAtz*dAtz;
    constexpr auto gB = dBsx*dBsx + dBsy*dBsy + dBsz*dBsz
                      - dBtx*dBtx - dBty*dBty - dBtz*dBtz;

    // query face normal n_c = (Qc-Pc) x (Rc-Pc)
    constexpr auto dQPx = Qcx - Pcx;  constexpr auto dQPy = Qcy - Pcy;  constexpr auto dQPz = Qcz - Pcz;
    constexpr auto dRPx = Rcx - Pcx;  constexpr auto dRPy = Rcy - Pcy;  constexpr auto dRPz = Rcz - Pcz;
    constexpr auto ncx = dQPy*dRPz - dQPz*dRPy;
    constexpr auto ncy = dQPz*dRPx - dQPx*dRPz;
    constexpr auto ncz = dQPx*dRPy - dQPy*dRPx;

    // query-face signed distances wA = n_c.(A-Pc), wB = n_c.(B-Pc)
    constexpr auto dAPx = Ax - Pcx;  constexpr auto dAPy = Ay - Pcy;  constexpr auto dAPz = Az - Pcz;
    constexpr auto dBPx = Bx - Pcx;  constexpr auto dBPy = By - Pcy;  constexpr auto dBPz = Bz - Pcz;
    constexpr auto wA = ncx*dAPx + ncy*dAPy + ncz*dAPz;
    constexpr auto wB = ncx*dBPx + ncy*dBPy + ncz*dBPz;

    // det[ gA wA ; gB wB ] = gA*wB - wA*gB
    using det2_t = grp::det<decltype(gA), decltype(wA),
                            decltype(gB), decltype(wB)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp3_feasible_TTS_T(double Ax, double Ay, double Az,
                                  double Bx, double By, double Bz,
                                  double sx, double sy, double sz,
                                  double tx, double ty, double tz,
                                  double Pcx, double Pcy, double Pcz,
                                  double Qcx, double Qcy, double Qcz,
                                  double Rcx, double Rcy, double Rcz)
{
    return - lp3_feasible_TTS_T_factor2_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz,
                                                           sx,sy,sz, tx,ty,tz,
                                                           Pcx,Pcy,Pcz, Qcx,Qcy,Qcz, Rcx,Rcy,Rcz)
           * lp3_feasible_TTS_S_factor1_impl::pred{}.apply(Ax,Ay,Az, Bx,By,Bz,
                                                           sx,sy,sz, tx,ty,tz, tx,ty,tz);
}


// ---------------------------------------------------------------------------
// lp3_feasible_TSS_S -- vertex (T,S,S) [Voronoi edge s->t1 ^ s->t2, meeting
// tet face T1=(A,Q,R)] tested against a bisector s->u.
// feasible = sign(Delta_TSS) * sign(Gamma),  Delta_TSS = lp3_D_TSS (factor1).
// Gamma = det[ n_T1 | 0 ; t1-s | c1 ; t2-s | c2 ; u-s | cu ]  (reference A on T1)
//   n_T1 = (Q-A) x (R-A);  cY = |Y-A|^2 - |s-A|^2.
// Rows ordered (n_T1, t1-s, t2-s) to match lp3_D_TSS's Delta exactly.
// Inputs: A(_1-_3), Q(_4-_6), R(_7-_9), s(_10-_12),
//         t1(_13-_15), t2(_16-_18), u(_19-_21)
// ---------------------------------------------------------------------------
namespace lp3_feasible_TSS_S_factor2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Qx  = grp::_4;
    constexpr auto Qy  = grp::_5;
    constexpr auto Qz  = grp::_6;
    constexpr auto Rx  = grp::_7;
    constexpr auto Ry  = grp::_8;
    constexpr auto Rz  = grp::_9;
    constexpr auto sx  = grp::_10;
    constexpr auto sy  = grp::_11;
    constexpr auto sz  = grp::_12;
    constexpr auto t1x = grp::_13;
    constexpr auto t1y = grp::_14;
    constexpr auto t1z = grp::_15;
    constexpr auto t2x = grp::_16;
    constexpr auto t2y = grp::_17;
    constexpr auto t2z = grp::_18;
    constexpr auto ux  = grp::_19;
    constexpr auto uy  = grp::_20;
    constexpr auto uz  = grp::_21;

    // n_T1 = (Q-A) x (R-A)
    constexpr auto dQAx = Qx - Ax;  constexpr auto dQAy = Qy - Ay;  constexpr auto dQAz = Qz - Az;
    constexpr auto dRAx = Rx - Ax;  constexpr auto dRAy = Ry - Ay;  constexpr auto dRAz = Rz - Az;
    constexpr auto nTx = dQAy*dRAz - dQAz*dRAy;
    constexpr auto nTy = dQAz*dRAx - dQAx*dRAz;
    constexpr auto nTz = dQAx*dRAy - dQAy*dRAx;

    // normal rows t1-s, t2-s, u-s
    constexpr auto d1x = t1x - sx;  constexpr auto d1y = t1y - sy;  constexpr auto d1z = t1z - sz;
    constexpr auto d2x = t2x - sx;  constexpr auto d2y = t2y - sy;  constexpr auto d2z = t2z - sz;
    constexpr auto dux = ux  - sx;  constexpr auto duy = uy  - sy;  constexpr auto duz = uz  - sz;

    // offset entries cY = |Y-A|^2 - |s-A|^2
    constexpr auto dsAx  = sx  - Ax;  constexpr auto dsAy  = sy  - Ay;  constexpr auto dsAz  = sz  - Az;
    constexpr auto dt1Ax = t1x - Ax;  constexpr auto dt1Ay = t1y - Ay;  constexpr auto dt1Az = t1z - Az;
    constexpr auto dt2Ax = t2x - Ax;  constexpr auto dt2Ay = t2y - Ay;  constexpr auto dt2Az = t2z - Az;
    constexpr auto duAx  = ux  - Ax;  constexpr auto duAy  = uy  - Ay;  constexpr auto duAz  = uz  - Az;

    constexpr auto sA2 = dsAx*dsAx + dsAy*dsAy + dsAz*dsAz;
    constexpr auto c1  = dt1Ax*dt1Ax + dt1Ay*dt1Ay + dt1Az*dt1Az - sA2;
    constexpr auto c2  = dt2Ax*dt2Ax + dt2Ay*dt2Ay + dt2Az*dt2Az - sA2;
    constexpr auto cu  = duAx*duAx   + duAy*duAy   + duAz*duAz   - sA2;

    // Gamma with rows (t1-s|c1),(t2-s|c2),(n_T1|0),(u-s|cu), cofactor-expanded
    // along the offset column (the n_T1 row is zero and drops out):
    //   Gamma = -c1*det[t2-s;n_T1;u-s] + c2*det[t1-s;n_T1;u-s] + cu*det[t1-s;t2-s;n_T1]
    constexpr auto M1 = grp::det<decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nTx),decltype(nTy),decltype(nTz),
                                 decltype(dux),decltype(duy),decltype(duz)>{};
    constexpr auto M2 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(nTx),decltype(nTy),decltype(nTz),
                                 decltype(dux),decltype(duy),decltype(duz)>{};
    constexpr auto M3 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nTx),decltype(nTy),decltype(nTz)>{};
    constexpr auto expr = c2*M2 + cu*M3 - c1*M1;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp3_feasible_TSS_S(double Ax, double Ay, double Az,
                                  double Qx, double Qy, double Qz,
                                  double Rx, double Ry, double Rz,
                                  double sx, double sy, double sz,
                                  double t1x, double t1y, double t1z,
                                  double t2x, double t2y, double t2z,
                                  double ux, double uy, double uz)
{
    return lp3_D_TSS_impl::pred{}.apply(Ax,Ay,Az, Qx,Qy,Qz, Rx,Ry,Rz,
                                        sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z)
         * lp3_feasible_TSS_S_factor2_impl::pred{}.apply(Ax,Ay,Az, Qx,Qy,Qz, Rx,Ry,Rz,
                                        sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z, ux,uy,uz);
}


// ---------------------------------------------------------------------------
// lp3_feasible_TSS_T -- vertex (T,S,S) [Voronoi edge s->t1 ^ s->t2, meeting the
// vertex-defining face T1] tested against a second, adjacent face T2.
// The two faces share the tet edge A-E; A is the reference on both planes.
//   T1 (vertex-definer) = (A,E,C), n_T1 = (E-A) x (C-A)
//   T2 (query)          = (A,E,D), n_T2 = (E-A) x (D-A)
// feasible = sign(Delta_TSS) * sign(Gamma),  Delta_TSS = lp3_D_TSS(A,E,C,s,t1,t2).
// Gamma (offset column cofactor-expanded, both T rows zero) =
//   -c1*det[t2-s; n_T1; n_T2] + c2*det[t1-s; n_T1; n_T2],  cY=|Y-A|^2-|s-A|^2.
// NOTE: n_T2 uses the fixed order (E-A)x(D-A); this returns feasibility relative
// to that raw orientation. The caller multiplies by the query face's interior
// orientation sign sigma_f (see lp3_predicates.tex, Face-orientation convention).
// Inputs: A(_1-_3), E(_4-_6), C(_7-_9), D(_10-_12), s(_13-_15),
//         t1(_16-_18), t2(_19-_21)
// ---------------------------------------------------------------------------
namespace lp3_feasible_TSS_T_factor2_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto Ex  = grp::_4;
    constexpr auto Ey  = grp::_5;
    constexpr auto Ez  = grp::_6;
    constexpr auto Cx  = grp::_7;
    constexpr auto Cy  = grp::_8;
    constexpr auto Cz  = grp::_9;
    constexpr auto Dx  = grp::_10;
    constexpr auto Dy  = grp::_11;
    constexpr auto Dz  = grp::_12;
    constexpr auto sx  = grp::_13;
    constexpr auto sy  = grp::_14;
    constexpr auto sz  = grp::_15;
    constexpr auto t1x = grp::_16;
    constexpr auto t1y = grp::_17;
    constexpr auto t1z = grp::_18;
    constexpr auto t2x = grp::_19;
    constexpr auto t2y = grp::_20;
    constexpr auto t2z = grp::_21;

    // shared edge direction and the two apex directions from A
    constexpr auto dEAx = Ex - Ax;  constexpr auto dEAy = Ey - Ay;  constexpr auto dEAz = Ez - Az;
    constexpr auto dCAx = Cx - Ax;  constexpr auto dCAy = Cy - Ay;  constexpr auto dCAz = Cz - Az;
    constexpr auto dDAx = Dx - Ax;  constexpr auto dDAy = Dy - Ay;  constexpr auto dDAz = Dz - Az;

    // n_T1 = (E-A) x (C-A),  n_T2 = (E-A) x (D-A)
    constexpr auto nT1x = dEAy*dCAz - dEAz*dCAy;
    constexpr auto nT1y = dEAz*dCAx - dEAx*dCAz;
    constexpr auto nT1z = dEAx*dCAy - dEAy*dCAx;
    constexpr auto nT2x = dEAy*dDAz - dEAz*dDAy;
    constexpr auto nT2y = dEAz*dDAx - dEAx*dDAz;
    constexpr auto nT2z = dEAx*dDAy - dEAy*dDAx;

    // t1-s, t2-s
    constexpr auto d1x = t1x - sx;  constexpr auto d1y = t1y - sy;  constexpr auto d1z = t1z - sz;
    constexpr auto d2x = t2x - sx;  constexpr auto d2y = t2y - sy;  constexpr auto d2z = t2z - sz;

    // cY = |Y-A|^2 - |s-A|^2
    constexpr auto dsAx  = sx  - Ax;  constexpr auto dsAy  = sy  - Ay;  constexpr auto dsAz  = sz  - Az;
    constexpr auto dt1Ax = t1x - Ax;  constexpr auto dt1Ay = t1y - Ay;  constexpr auto dt1Az = t1z - Az;
    constexpr auto dt2Ax = t2x - Ax;  constexpr auto dt2Ay = t2y - Ay;  constexpr auto dt2Az = t2z - Az;
    constexpr auto sA2 = dsAx*dsAx + dsAy*dsAy + dsAz*dsAz;
    constexpr auto c1  = dt1Ax*dt1Ax + dt1Ay*dt1Ay + dt1Az*dt1Az - sA2;
    constexpr auto c2  = dt2Ax*dt2Ax + dt2Ay*dt2Ay + dt2Az*dt2Az - sA2;

    // Gamma = -c1*det[t2-s;n_T1;n_T2] + c2*det[t1-s;n_T1;n_T2]
    constexpr auto M1 = grp::det<decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nT1x),decltype(nT1y),decltype(nT1z),
                                 decltype(nT2x),decltype(nT2y),decltype(nT2z)>{};
    constexpr auto M2 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(nT1x),decltype(nT1y),decltype(nT1z),
                                 decltype(nT2x),decltype(nT2y),decltype(nT2z)>{};
    constexpr auto expr = c2*M2 - c1*M1;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp3_feasible_TSS_T(double Ax, double Ay, double Az,
                                  double Ex, double Ey, double Ez,
                                  double Cx, double Cy, double Cz,
                                  double Dx, double Dy, double Dz,
                                  double sx, double sy, double sz,
                                  double t1x, double t1y, double t1z,
                                  double t2x, double t2y, double t2z)
{
    return lp3_D_TSS_impl::pred{}.apply(Ax,Ay,Az, Ex,Ey,Ez, Cx,Cy,Cz,
                                        sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z)
         * lp3_feasible_TSS_T_factor2_impl::pred{}.apply(Ax,Ay,Az, Ex,Ey,Ez, Cx,Cy,Cz,
                                        Dx,Dy,Dz, sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z);
}


// ---------------------------------------------------------------------------
// lp3_feasible_TSS_T_gen -- vertex (T,S,S) [Voronoi edge s->t1 ^ s->t2, meeting
// its definer face T1=(A,Q1,R1)] tested against an ARBITRARY query face
// T2=(P2,Q2,R2), not assumed to share an edge with T1 (surface edge ordering).
// Reference A on T1 (only T1's offset vanishes); see lp3_predicates.tex,
// Prop "TSS vs an arbitrary face".
//   n_T1=(Q1-A)x(R1-A);  n_T2=(Q2-P2)x(R2-P2);  cY=|Y-A|^2-|s-A|^2;
//   w = n_T2.(P2-A);  Delta_TSS = det[t1-s; t2-s; n_T1] = lp3_D_TSS(A,Q1,R1,...).
//   Gamma = -c1*det[t2-s;n_T1;n_T2] + c2*det[t1-s;n_T1;n_T2] + w*Delta_TSS.
// feasible = sign(Delta_TSS) * sign(Gamma).  Degree 7.  (= lp3_feasible_TSS_T
// when A can be put on both faces, i.e. w=0.)  Orientation of n_T2 is raw; the
// surface ordering uses feasibility*slope, which is invariant to it.
// Inputs: A(_1-_3), Q1(_4-_6), R1(_7-_9), P2(_10-_12), Q2(_13-_15), R2(_16-_18),
//         s(_19-_21), t1(_22-_24), t2(_25-_27)
// ---------------------------------------------------------------------------
namespace lp3_feasible_TSS_T_gen_factor2_impl {
    constexpr auto Ax  = grp::_1;   constexpr auto Ay  = grp::_2;   constexpr auto Az  = grp::_3;
    constexpr auto Q1x = grp::_4;   constexpr auto Q1y = grp::_5;   constexpr auto Q1z = grp::_6;
    constexpr auto R1x = grp::_7;   constexpr auto R1y = grp::_8;   constexpr auto R1z = grp::_9;
    constexpr auto P2x = grp::_10;  constexpr auto P2y = grp::_11;  constexpr auto P2z = grp::_12;
    constexpr auto Q2x = grp::_13;  constexpr auto Q2y = grp::_14;  constexpr auto Q2z = grp::_15;
    constexpr auto R2x = grp::_16;  constexpr auto R2y = grp::_17;  constexpr auto R2z = grp::_18;
    constexpr auto sx  = grp::_19;  constexpr auto sy  = grp::_20;  constexpr auto sz  = grp::_21;
    constexpr auto t1x = grp::_22;  constexpr auto t1y = grp::_23;  constexpr auto t1z = grp::_24;
    constexpr auto t2x = grp::_25;  constexpr auto t2y = grp::_26;  constexpr auto t2z = grp::_27;

    // n_T1 = (Q1-A) x (R1-A)
    constexpr auto dQ1Ax = Q1x-Ax; constexpr auto dQ1Ay = Q1y-Ay; constexpr auto dQ1Az = Q1z-Az;
    constexpr auto dR1Ax = R1x-Ax; constexpr auto dR1Ay = R1y-Ay; constexpr auto dR1Az = R1z-Az;
    constexpr auto nT1x = dQ1Ay*dR1Az - dQ1Az*dR1Ay;
    constexpr auto nT1y = dQ1Az*dR1Ax - dQ1Ax*dR1Az;
    constexpr auto nT1z = dQ1Ax*dR1Ay - dQ1Ay*dR1Ax;
    // n_T2 = (Q2-P2) x (R2-P2)
    constexpr auto dQ2Px = Q2x-P2x; constexpr auto dQ2Py = Q2y-P2y; constexpr auto dQ2Pz = Q2z-P2z;
    constexpr auto dR2Px = R2x-P2x; constexpr auto dR2Py = R2y-P2y; constexpr auto dR2Pz = R2z-P2z;
    constexpr auto nT2x = dQ2Py*dR2Pz - dQ2Pz*dR2Py;
    constexpr auto nT2y = dQ2Pz*dR2Px - dQ2Px*dR2Pz;
    constexpr auto nT2z = dQ2Px*dR2Py - dQ2Py*dR2Px;
    // t1-s, t2-s
    constexpr auto d1x = t1x-sx; constexpr auto d1y = t1y-sy; constexpr auto d1z = t1z-sz;
    constexpr auto d2x = t2x-sx; constexpr auto d2y = t2y-sy; constexpr auto d2z = t2z-sz;
    // cY = |Y-A|^2 - |s-A|^2
    constexpr auto dsAx = sx-Ax;   constexpr auto dsAy = sy-Ay;   constexpr auto dsAz = sz-Az;
    constexpr auto dt1Ax = t1x-Ax; constexpr auto dt1Ay = t1y-Ay; constexpr auto dt1Az = t1z-Az;
    constexpr auto dt2Ax = t2x-Ax; constexpr auto dt2Ay = t2y-Ay; constexpr auto dt2Az = t2z-Az;
    constexpr auto sA2 = dsAx*dsAx + dsAy*dsAy + dsAz*dsAz;
    constexpr auto c1  = dt1Ax*dt1Ax + dt1Ay*dt1Ay + dt1Az*dt1Az - sA2;
    constexpr auto c2  = dt2Ax*dt2Ax + dt2Ay*dt2Ay + dt2Az*dt2Az - sA2;
    // w = n_T2 . (P2 - A).  NOTE: the bisector offsets carry a factor 1/2 that we
    // dropped from c1,c2; to keep the offset COLUMN uniformly scaled the query-face
    // offset must be doubled (2w), i.e. the whole Gamma is scaled by 2 (sign kept).
    constexpr auto dP2Ax = P2x-Ax; constexpr auto dP2Ay = P2y-Ay; constexpr auto dP2Az = P2z-Az;
    constexpr auto w  = nT2x*dP2Ax + nT2y*dP2Ay + nT2z*dP2Az;
    constexpr auto w2 = int_const<2>{} * w;
    // Gamma(x2) = c2*det[t1-s;n_T1;n_T2] - c1*det[t2-s;n_T1;n_T2] + 2w*det[t1-s;t2-s;n_T1]
    constexpr auto M1 = grp::det<decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nT1x),decltype(nT1y),decltype(nT1z),
                                 decltype(nT2x),decltype(nT2y),decltype(nT2z)>{};
    constexpr auto M2 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(nT1x),decltype(nT1y),decltype(nT1z),
                                 decltype(nT2x),decltype(nT2y),decltype(nT2z)>{};
    constexpr auto Dexpr = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                    decltype(d2x),decltype(d2y),decltype(d2z),
                                    decltype(nT1x),decltype(nT1y),decltype(nT1z)>{};
    constexpr auto expr = c2*M2 - c1*M1 + w2*Dexpr;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp3_feasible_TSS_T_gen(double Ax, double Ay, double Az,
                                      double Q1x, double Q1y, double Q1z,
                                      double R1x, double R1y, double R1z,
                                      double P2x, double P2y, double P2z,
                                      double Q2x, double Q2y, double Q2z,
                                      double R2x, double R2y, double R2z,
                                      double sx, double sy, double sz,
                                      double t1x, double t1y, double t1z,
                                      double t2x, double t2y, double t2z)
{
    return lp3_D_TSS_impl::pred{}.apply(Ax,Ay,Az, Q1x,Q1y,Q1z, R1x,R1y,R1z,
                                        sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z)
         * lp3_feasible_TSS_T_gen_factor2_impl::pred{}.apply(
                Ax,Ay,Az, Q1x,Q1y,Q1z, R1x,R1y,R1z,
                P2x,P2y,P2z, Q2x,Q2y,Q2z, R2x,R2y,R2z,
                sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z);
}


// ---------------------------------------------------------------------------
// lp3_feasible_SSS_T -- vertex (S,S,S) [circumcenter of Delaunay tet
// (s,t1,t2,t3)] tested against a tet face T=(A,Q,R).
// feasible = sign(Delta_SSS) * sign(Gamma).
//   factor1 = Delta_SSS = det[t1-s; t2-s; t3-s]  (= orient3d(s,t1,t2,t3),
//             computed locally so its sign convention matches Gamma exactly).
//   n_T = (Q-A) x (R-A);  cY = |Y-A|^2 - |s-A|^2  (reference A on the face).
//   Gamma (offset column cofactor-expanded, the n_T row is zero) =
//     -c1*det[t2-s;t3-s;n_T] + c2*det[t1-s;t3-s;n_T] - c3*det[t1-s;t2-s;n_T].
// NOTE: query is a T face with fixed n_T=(Q-A)x(R-A); the caller multiplies by
// the face interior-orientation sign sigma_f (Face-orientation convention).
// Inputs: s(_1-_3), t1(_4-_6), t2(_7-_9), t3(_10-_12),
//         A(_13-_15), Q(_16-_18), R(_19-_21)
// ---------------------------------------------------------------------------
namespace lp3_feasible_SSS_T_factor1_impl {
    constexpr auto sx  = grp::_1;
    constexpr auto sy  = grp::_2;
    constexpr auto sz  = grp::_3;
    constexpr auto t1x = grp::_4;
    constexpr auto t1y = grp::_5;
    constexpr auto t1z = grp::_6;
    constexpr auto t2x = grp::_7;
    constexpr auto t2y = grp::_8;
    constexpr auto t2z = grp::_9;
    constexpr auto t3x = grp::_10;
    constexpr auto t3y = grp::_11;
    constexpr auto t3z = grp::_12;
    // _13.._21 = A,Q,R unused

    constexpr auto d1x = t1x - sx;  constexpr auto d1y = t1y - sy;  constexpr auto d1z = t1z - sz;
    constexpr auto d2x = t2x - sx;  constexpr auto d2y = t2y - sy;  constexpr auto d2z = t2z - sz;
    constexpr auto d3x = t3x - sx;  constexpr auto d3y = t3y - sy;  constexpr auto d3z = t3z - sz;

    using expr_t = grp::det<
        decltype(d1x), decltype(d1y), decltype(d1z),
        decltype(d2x), decltype(d2y), decltype(d2z),
        decltype(d3x), decltype(d3y), decltype(d3z)
    >;
    constexpr auto expr = expr_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

namespace lp3_feasible_SSS_T_factor2_impl {
    constexpr auto sx  = grp::_1;
    constexpr auto sy  = grp::_2;
    constexpr auto sz  = grp::_3;
    constexpr auto t1x = grp::_4;
    constexpr auto t1y = grp::_5;
    constexpr auto t1z = grp::_6;
    constexpr auto t2x = grp::_7;
    constexpr auto t2y = grp::_8;
    constexpr auto t2z = grp::_9;
    constexpr auto t3x = grp::_10;
    constexpr auto t3y = grp::_11;
    constexpr auto t3z = grp::_12;
    constexpr auto Ax  = grp::_13;
    constexpr auto Ay  = grp::_14;
    constexpr auto Az  = grp::_15;
    constexpr auto Qx  = grp::_16;
    constexpr auto Qy  = grp::_17;
    constexpr auto Qz  = grp::_18;
    constexpr auto Rx  = grp::_19;
    constexpr auto Ry  = grp::_20;
    constexpr auto Rz  = grp::_21;

    // n_T = (Q-A) x (R-A)
    constexpr auto dQAx = Qx - Ax;  constexpr auto dQAy = Qy - Ay;  constexpr auto dQAz = Qz - Az;
    constexpr auto dRAx = Rx - Ax;  constexpr auto dRAy = Ry - Ay;  constexpr auto dRAz = Rz - Az;
    constexpr auto nTx = dQAy*dRAz - dQAz*dRAy;
    constexpr auto nTy = dQAz*dRAx - dQAx*dRAz;
    constexpr auto nTz = dQAx*dRAy - dQAy*dRAx;

    // t_i - s
    constexpr auto d1x = t1x - sx;  constexpr auto d1y = t1y - sy;  constexpr auto d1z = t1z - sz;
    constexpr auto d2x = t2x - sx;  constexpr auto d2y = t2y - sy;  constexpr auto d2z = t2z - sz;
    constexpr auto d3x = t3x - sx;  constexpr auto d3y = t3y - sy;  constexpr auto d3z = t3z - sz;

    // cY = |Y-A|^2 - |s-A|^2
    constexpr auto dsAx  = sx  - Ax;  constexpr auto dsAy  = sy  - Ay;  constexpr auto dsAz  = sz  - Az;
    constexpr auto dt1Ax = t1x - Ax;  constexpr auto dt1Ay = t1y - Ay;  constexpr auto dt1Az = t1z - Az;
    constexpr auto dt2Ax = t2x - Ax;  constexpr auto dt2Ay = t2y - Ay;  constexpr auto dt2Az = t2z - Az;
    constexpr auto dt3Ax = t3x - Ax;  constexpr auto dt3Ay = t3y - Ay;  constexpr auto dt3Az = t3z - Az;
    constexpr auto sA2 = dsAx*dsAx + dsAy*dsAy + dsAz*dsAz;
    constexpr auto c1  = dt1Ax*dt1Ax + dt1Ay*dt1Ay + dt1Az*dt1Az - sA2;
    constexpr auto c2  = dt2Ax*dt2Ax + dt2Ay*dt2Ay + dt2Az*dt2Az - sA2;
    constexpr auto c3  = dt3Ax*dt3Ax + dt3Ay*dt3Ay + dt3Az*dt3Az - sA2;

    // Gamma = -c1*det[d2;d3;nT] + c2*det[d1;d3;nT] - c3*det[d1;d2;nT]
    constexpr auto M1 = grp::det<decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(d3x),decltype(d3y),decltype(d3z),
                                 decltype(nTx),decltype(nTy),decltype(nTz)>{};
    constexpr auto M2 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(d3x),decltype(d3y),decltype(d3z),
                                 decltype(nTx),decltype(nTy),decltype(nTz)>{};
    constexpr auto M3 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nTx),decltype(nTy),decltype(nTz)>{};
    constexpr auto expr = c2*M2 - c1*M1 - c3*M3;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int lp3_feasible_SSS_T(double sx, double sy, double sz,
                                  double t1x, double t1y, double t1z,
                                  double t2x, double t2y, double t2z,
                                  double t3x, double t3y, double t3z,
                                  double Ax, double Ay, double Az,
                                  double Qx, double Qy, double Qz,
                                  double Rx, double Ry, double Rz)
{
    return lp3_feasible_SSS_T_factor1_impl::pred{}.apply(
                sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z, t3x,t3y,t3z,
                Ax,Ay,Az, Qx,Qy,Qz, Rx,Ry,Rz)
         * lp3_feasible_SSS_T_factor2_impl::pred{}.apply(
                sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z, t3x,t3y,t3z,
                Ax,Ay,Az, Qx,Qy,Qz, Rx,Ry,Rz);
}


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

// robust 2D orientation (Shewchuk), for the collinearity guard.
extern "C" int orient2d(double, double, double, double, double, double);

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
using _size_lp3_D_TSS     = show_stage_d_size<lp3_D_TSS_impl::exact::results_size>*;
using _size_lp3_feasible_TTS_S_f1 = show_stage_d_size<lp3_feasible_TTS_S_factor1_impl::exact::results_size>*;
using _size_lp3_feasible_TTS_S_f2 = show_stage_d_size<lp3_feasible_TTS_S_factor2_impl::exact::results_size>*;
using _size_lp3_feasible_TSS_S_f2 = show_stage_d_size<lp3_feasible_TSS_S_factor2_impl::exact::results_size>*;
using _size_lp3_feasible_TSS_T_f2 = show_stage_d_size<lp3_feasible_TSS_T_factor2_impl::exact::results_size>*;
using _size_lp3_feasible_SSS_T_f1 = show_stage_d_size<lp3_feasible_SSS_T_factor1_impl::exact::results_size>*;
using _size_lp3_feasible_SSS_T_f2 = show_stage_d_size<lp3_feasible_SSS_T_factor2_impl::exact::results_size>*;
using _size_incircle3d = show_stage_d_size<incircle3d_impl::exact::results_size>*;
using _size_orient3d_dd = show_stage_d_size<orient3d_dd_impl::exact::results_size>*;
// lp_feasible_T1_S_T2 reuses lp_D_T1_S_impl and lp_feasible_T1_T2_S_impl -- no new size entry
// lp_feasible_T0_S_T2 reuses lp_D_T0_S_impl and lp_feasible_T0_T2_S_impl -- no new size entry
#endif
