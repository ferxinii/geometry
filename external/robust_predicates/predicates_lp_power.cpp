#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expressions.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_d.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/staged_predicate.hpp"

#include "predicates_lp.h"   /* lp3_slope_TTS, lp3_D_TSS: the slope-only
                              * factors are weight-invariant (a_S untouched),
                              * so they are reused UNWEIGHTED from here rather
                              * than re-derived -- see geometric_tests.tex,
                              * appendix:LP, the a_S-untouched proof. */

namespace grp = boost::geometry::detail::generic_robust_predicates;

// Compile-time integer constant leaf (see predicates_lp_n3.cpp: same role,
// duplicated here since that file's anonymous helper isn't shared via a
// header -- this is the only cross-file duplication, everything else reuses
// the existing exported lp3_slope_TTS/lp3_D_TSS functions as black boxes).
template <int N>
struct int_const : grp::static_constant_interface<double> {
    static constexpr double value = static_cast<double>(N);
    static constexpr bool non_negative = (N >= 0);
};

// LIBRARY BUG WORKAROUND (expansion_arithmetic.hpp): `<expansion> - <scalar
// leaf>` lowers to expansion_minus(InIter,InIter,Real,...), which forwards to
// an expansion_plus<...,e_negate,f_negate,...> instantiation that does not
// exist for a scalar right-hand operand (confirmed by isolated repro: fails
// to compile with "no matching function for call to 'expansion_plus'" at
// exactly this expansion-minus-scalar shape, for ANY expansion length -- not
// a degree/complexity limit). `<scalar> - <scalar>` and `<expansion> +
// <scalar>` both go through different, working code paths. Every weight
// substitution below is therefore written as a scalar-minus-scalar
// difference computed FIRST (e.g. `dwt = wt - ws`), then added once to the
// expansion (`gX + dwt`) -- never `gX - ws + wt`. Mathematically identical,
// same exactness guarantees (the scalar subtraction of two raw leaves is the
// same primitive operation already used throughout predicates_lp_n3.cpp,
// e.g. `dAsx = Ax - sx`); only the operation ORDER changes to dodge the bug.


// ===========================================================================
// power_lp2_feasible_T0_T1_S -- weighted analogue of lp2_feasible_T0_T1_S.
// sign(|ti-A|^2 - |s-A|^2 + ws - wt).
// Derivation: unweighted feasible(TTT;S) = -sign(|V-s|^2-|V-t|^2)
// (geometric_tests.tex, prop:tet-case-table, proof (a)). Weighted power
// bisector substitutes |x-s|^2-w_s <= |x-t|^2-w_t for the favor-s condition,
// i.e. the same value shifted by -w_s+w_t before the sign/comparison:
//   feasible^w = -sign(|V-s|^2-w_s - |V-t|^2+w_t)
//             =  sign(|V-t|^2 - |V-s|^2 + ws - wt)
//             =  sign(expr + ws - wt),  expr = the existing unweighted code's
// |ti-A|^2-|s-A|^2 (V=A, t=ti). At ws=wt=0 this is exactly the unweighted
// expr, so sign agreement at zero weight is definitional, not incidental.
// Inputs: A(_1-_3), s(_4-_6), ti(_7-_9), ws(_10), wt(_11)
// ===========================================================================
namespace power_lp2_feasible_T0_T1_S_impl {
    constexpr auto Ax  = grp::_1;
    constexpr auto Ay  = grp::_2;
    constexpr auto Az  = grp::_3;
    constexpr auto sx  = grp::_4;
    constexpr auto sy  = grp::_5;
    constexpr auto sz  = grp::_6;
    constexpr auto tix = grp::_7;
    constexpr auto tiy = grp::_8;
    constexpr auto tiz = grp::_9;
    constexpr auto ws  = grp::_10;
    constexpr auto wt  = grp::_11;

    constexpr auto dsx  = sx  - Ax;
    constexpr auto dsy  = sy  - Ay;
    constexpr auto dsz  = sz  - Az;
    constexpr auto dtix = tix - Ax;
    constexpr auto dtiy = tiy - Ay;
    constexpr auto dtiz = tiz - Az;

    constexpr auto dw = ws - wt;
    constexpr auto expr = dtix*dtix + dtiy*dtiy + dtiz*dtiz
                        - dsx*dsx - dsy*dsy - dsz*dsz
                        + dw;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp2_feasible_T0_T1_S(double Ax, double Ay, double Az,
                                          double sx, double sy, double sz,
                                          double tix, double tiy, double tiz,
                                          double ws, double wt)
{
    return power_lp2_feasible_T0_T1_S_impl::pred{}.apply(Ax,Ay,Az, sx,sy,sz,
                                                         tix,tiy,tiz, ws,wt);
}


// ===========================================================================
// power_lp3_feasible_TTS_S / power_lp3_det_TTS -- shared weighted factor2.
// feasible = -sign(gA^w*hB^w - gB^w*hA^w) * sign(gA-gB), same structure as
// the unweighted lp3_feasible_TTS_S, with the CRITICAL-RULE leaf substitution
// (geometric_tests.tex a_S-untouched proof): gX^w = gX - ws + wt,
// hX^w = hX - ws + wu.  This is NOT "old value plus a constant" -- expanding
// gA^w*hB^w - gB^w*hA^w produces cross terms:
//   (gA-ws+wt)(hB-ws+wu) - (gB-ws+wt)(hA-ws+wu)
//     = (gA*hB - gB*hA) + (wu-ws)(gA-gB) + (wt-ws)(hB-hA) + ... (further
//       cross terms in ws that do NOT cancel in general), so gX^w, hX^w are
//       substituted at the leaf and the det<> template expansion below
//       produces the correct expanded polynomial automatically -- no
//       hand-derived top-level correction term is used.
// power_lp3_det_TTS reuses this exact same weighted product with (t1,t2) in
// place of (t,u) (g1 playing the role of g, g2 the role of h), mirroring how
// the unweighted lp3_det_TTS reuses lp3_feasible_TTS_S_factor2_impl.
// Inputs: A(_1-_3), B(_4-_6), s(_7-_9), t(_10-_12), u(_13-_15),
//         ws(_16), wt(_17), wu(_18)
// ===========================================================================
namespace power_lp3_feasible_TTS_S_factor2_impl {
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
    constexpr auto ws = grp::_16;
    constexpr auto wt = grp::_17;
    constexpr auto wu = grp::_18;

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

    // Leaf substitution (a_S untouched, only the bisector OFFSET shifts).
    // dwt,dwu computed scalar-minus-scalar first, then added once (see the
    // library-bug workaround note above).
    constexpr auto dwt = wt - ws;
    constexpr auto dwu = wu - ws;
    constexpr auto gAw = gA + dwt;
    constexpr auto gBw = gB + dwt;
    constexpr auto hAw = hA + dwu;
    constexpr auto hBw = hB + dwu;

    // det[ gAw hAw ; gBw hBw ] = gAw*hBw - hAw*gBw -- template expands the
    // cross terms in ws,wt,wu automatically.
    using det2_t = grp::det<decltype(gAw), decltype(hAw),
                            decltype(gBw), decltype(hBw)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp3_feasible_TTS_S(double Ax, double Ay, double Az,
                                        double Bx, double By, double Bz,
                                        double sx, double sy, double sz,
                                        double tx, double ty, double tz,
                                        double ux, double uy, double uz,
                                        double ws, double wt, double wu)
{
    return - power_lp3_feasible_TTS_S_factor2_impl::pred{}.apply(
                    Ax,Ay,Az, Bx,By,Bz, sx,sy,sz, tx,ty,tz, ux,uy,uz, ws,wt,wu)
           * lp3_slope_TTS(Ax,Ay,Az, Bx,By,Bz, sx,sy,sz, tx,ty,tz);
}

extern "C" int power_lp3_det_TTS(double Ax, double Ay, double Az,
                                 double Bx, double By, double Bz,
                                 double sx, double sy, double sz,
                                 double t1x, double t1y, double t1z,
                                 double t2x, double t2y, double t2z,
                                 double ws, double wt1, double wt2)
{
    return power_lp3_feasible_TTS_S_factor2_impl::pred{}.apply(
                Ax,Ay,Az, Bx,By,Bz, sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z,
                ws,wt1,wt2);
}


// ===========================================================================
// power_lp3_feasible_TTS_T -- weighted analogue of lp3_feasible_TTS_T.
// feasible = -sign(gA^w*wB - gB^w*wA) * sign(gA-gB).  gX^w = gX-ws+wt as in
// TTS_S; wX (query T-face signed distance, purely geometric) is
// weight-independent and untouched -- same leaf-substitution-then-expand
// rule as TTS_S, just with h replaced by w.
// Inputs: A(_1-_3), B(_4-_6), s(_7-_9), t(_10-_12),
//         Pc(_13-_15), Qc(_16-_18), Rc(_19-_21), ws(_22), wt(_23)
// ===========================================================================
namespace power_lp3_feasible_TTS_T_factor2_impl {
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
    constexpr auto ws  = grp::_22;
    constexpr auto wt  = grp::_23;

    constexpr auto dAsx = Ax - sx;  constexpr auto dAsy = Ay - sy;  constexpr auto dAsz = Az - sz;
    constexpr auto dAtx = Ax - tx;  constexpr auto dAty = Ay - ty;  constexpr auto dAtz = Az - tz;
    constexpr auto dBsx = Bx - sx;  constexpr auto dBsy = By - sy;  constexpr auto dBsz = Bz - sz;
    constexpr auto dBtx = Bx - tx;  constexpr auto dBty = By - ty;  constexpr auto dBtz = Bz - tz;

    constexpr auto gA = dAsx*dAsx + dAsy*dAsy + dAsz*dAsz
                      - dAtx*dAtx - dAty*dAty - dAtz*dAtz;
    constexpr auto gB = dBsx*dBsx + dBsy*dBsy + dBsz*dBsz
                      - dBtx*dBtx - dBty*dBty - dBtz*dBtz;
    constexpr auto dwt = wt - ws;
    constexpr auto gAw = gA + dwt;
    constexpr auto gBw = gB + dwt;

    constexpr auto dQPx = Qcx - Pcx;  constexpr auto dQPy = Qcy - Pcy;  constexpr auto dQPz = Qcz - Pcz;
    constexpr auto dRPx = Rcx - Pcx;  constexpr auto dRPy = Rcy - Pcy;  constexpr auto dRPz = Rcz - Pcz;
    constexpr auto ncx = dQPy*dRPz - dQPz*dRPy;
    constexpr auto ncy = dQPz*dRPx - dQPx*dRPz;
    constexpr auto ncz = dQPx*dRPy - dQPy*dRPx;

    constexpr auto dAPx = Ax - Pcx;  constexpr auto dAPy = Ay - Pcy;  constexpr auto dAPz = Az - Pcz;
    constexpr auto dBPx = Bx - Pcx;  constexpr auto dBPy = By - Pcy;  constexpr auto dBPz = Bz - Pcz;
    constexpr auto wA = ncx*dAPx + ncy*dAPy + ncz*dAPz;
    constexpr auto wB = ncx*dBPx + ncy*dBPy + ncz*dBPz;

    // det[ gAw wA ; gBw wB ] = gAw*wB - wA*gBw
    using det2_t = grp::det<decltype(gAw), decltype(wA),
                            decltype(gBw), decltype(wB)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp3_feasible_TTS_T(double Ax, double Ay, double Az,
                                        double Bx, double By, double Bz,
                                        double sx, double sy, double sz,
                                        double tx, double ty, double tz,
                                        double Pcx, double Pcy, double Pcz,
                                        double Qcx, double Qcy, double Qcz,
                                        double Rcx, double Rcy, double Rcz,
                                        double ws, double wt)
{
    return - power_lp3_feasible_TTS_T_factor2_impl::pred{}.apply(
                    Ax,Ay,Az, Bx,By,Bz, sx,sy,sz, tx,ty,tz,
                    Pcx,Pcy,Pcz, Qcx,Qcy,Qcz, Rcx,Rcy,Rcz, ws,wt)
           * lp3_slope_TTS(Ax,Ay,Az, Bx,By,Bz, sx,sy,sz, tx,ty,tz);
}


// ===========================================================================
// power_lp3_feasible_TSS_S -- weighted analogue of lp3_feasible_TSS_S.
// feasible = sign(Delta_TSS) * sign(Gamma^w).  Delta_TSS (the slope,
// weight-invariant) is reused unweighted via lp3_D_TSS.  Gamma is LINEAR in
// c1,c2,cu (a sum, each offset appearing to the first power) -- per the
// critical-rule split this is the direct-substitution case, no cross terms:
//   ci^w = ci + ws - w_ti  (resp. cu^w = cu + ws - wu), the same M1,M2,M3
// geometric minors (pure orientation on t1-s,t2-s,u-s,n_T1: no weight
// dependence) reused unweighted.
// Inputs: A(_1-_3), Q(_4-_6), R(_7-_9), s(_10-_12),
//         t1(_13-_15), t2(_16-_18), u(_19-_21),
//         ws(_22), wt1(_23), wt2(_24), wu(_25)
// ===========================================================================
namespace power_lp3_feasible_TSS_S_factor2_impl {
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
    constexpr auto ws  = grp::_22;
    constexpr auto wt1 = grp::_23;
    constexpr auto wt2 = grp::_24;
    constexpr auto wu  = grp::_25;

    constexpr auto dQAx = Qx - Ax;  constexpr auto dQAy = Qy - Ay;  constexpr auto dQAz = Qz - Az;
    constexpr auto dRAx = Rx - Ax;  constexpr auto dRAy = Ry - Ay;  constexpr auto dRAz = Rz - Az;
    constexpr auto nTx = dQAy*dRAz - dQAz*dRAy;
    constexpr auto nTy = dQAz*dRAx - dQAx*dRAz;
    constexpr auto nTz = dQAx*dRAy - dQAy*dRAx;

    constexpr auto d1x = t1x - sx;  constexpr auto d1y = t1y - sy;  constexpr auto d1z = t1z - sz;
    constexpr auto d2x = t2x - sx;  constexpr auto d2y = t2y - sy;  constexpr auto d2z = t2z - sz;
    constexpr auto dux = ux  - sx;  constexpr auto duy = uy  - sy;  constexpr auto duz = uz  - sz;

    constexpr auto dsAx  = sx  - Ax;  constexpr auto dsAy  = sy  - Ay;  constexpr auto dsAz  = sz  - Az;
    constexpr auto dt1Ax = t1x - Ax;  constexpr auto dt1Ay = t1y - Ay;  constexpr auto dt1Az = t1z - Az;
    constexpr auto dt2Ax = t2x - Ax;  constexpr auto dt2Ay = t2y - Ay;  constexpr auto dt2Az = t2z - Az;
    constexpr auto duAx  = ux  - Ax;  constexpr auto duAy  = uy  - Ay;  constexpr auto duAz  = uz  - Az;

    constexpr auto sA2 = dsAx*dsAx + dsAy*dsAy + dsAz*dsAz;
    constexpr auto c1  = dt1Ax*dt1Ax + dt1Ay*dt1Ay + dt1Az*dt1Az - sA2;
    constexpr auto c2  = dt2Ax*dt2Ax + dt2Ay*dt2Ay + dt2Az*dt2Az - sA2;
    constexpr auto cu  = duAx*duAx   + duAy*duAy   + duAz*duAz   - sA2;

    // Direct substitution (linear-in-c group): ci^w = ci + (ws - w_ti).
    constexpr auto dw1 = ws - wt1;
    constexpr auto dw2 = ws - wt2;
    constexpr auto dwu = ws - wu;
    constexpr auto c1w = c1 + dw1;
    constexpr auto c2w = c2 + dw2;
    constexpr auto cuw = cu + dwu;

    constexpr auto M1 = grp::det<decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nTx),decltype(nTy),decltype(nTz),
                                 decltype(dux),decltype(duy),decltype(duz)>{};
    constexpr auto M2 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(nTx),decltype(nTy),decltype(nTz),
                                 decltype(dux),decltype(duy),decltype(duz)>{};
    constexpr auto M3 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nTx),decltype(nTy),decltype(nTz)>{};
    constexpr auto expr = c2w*M2 + cuw*M3 - c1w*M1;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp3_feasible_TSS_S(double Ax, double Ay, double Az,
                                        double Qx, double Qy, double Qz,
                                        double Rx, double Ry, double Rz,
                                        double sx, double sy, double sz,
                                        double t1x, double t1y, double t1z,
                                        double t2x, double t2y, double t2z,
                                        double ux, double uy, double uz,
                                        double ws, double wt1, double wt2, double wu)
{
    return lp3_D_TSS(Ax,Ay,Az, Qx,Qy,Qz, Rx,Ry,Rz,
                     sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z)
         * power_lp3_feasible_TSS_S_factor2_impl::pred{}.apply(
                Ax,Ay,Az, Qx,Qy,Qz, Rx,Ry,Rz, sx,sy,sz,
                t1x,t1y,t1z, t2x,t2y,t2z, ux,uy,uz, ws,wt1,wt2,wu);
}


// ===========================================================================
// power_lp3_feasible_TSS_T -- weighted analogue of lp3_feasible_TSS_T.
// Same c1^w, c2^w substitution as TSS_S; both T-face minors M1,M2 are pure
// orientation (n_T1,n_T2 built only from A,E,C,D), untouched.
// Inputs: A(_1-_3), E(_4-_6), C(_7-_9), D(_10-_12), s(_13-_15),
//         t1(_16-_18), t2(_19-_21), ws(_22), wt1(_23), wt2(_24)
// ===========================================================================
namespace power_lp3_feasible_TSS_T_factor2_impl {
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
    constexpr auto ws  = grp::_22;
    constexpr auto wt1 = grp::_23;
    constexpr auto wt2 = grp::_24;

    constexpr auto dEAx = Ex - Ax;  constexpr auto dEAy = Ey - Ay;  constexpr auto dEAz = Ez - Az;
    constexpr auto dCAx = Cx - Ax;  constexpr auto dCAy = Cy - Ay;  constexpr auto dCAz = Cz - Az;
    constexpr auto dDAx = Dx - Ax;  constexpr auto dDAy = Dy - Ay;  constexpr auto dDAz = Dz - Az;

    constexpr auto nT1x = dEAy*dCAz - dEAz*dCAy;
    constexpr auto nT1y = dEAz*dCAx - dEAx*dCAz;
    constexpr auto nT1z = dEAx*dCAy - dEAy*dCAx;
    constexpr auto nT2x = dEAy*dDAz - dEAz*dDAy;
    constexpr auto nT2y = dEAz*dDAx - dEAx*dDAz;
    constexpr auto nT2z = dEAx*dDAy - dEAy*dDAx;

    constexpr auto d1x = t1x - sx;  constexpr auto d1y = t1y - sy;  constexpr auto d1z = t1z - sz;
    constexpr auto d2x = t2x - sx;  constexpr auto d2y = t2y - sy;  constexpr auto d2z = t2z - sz;

    constexpr auto dsAx  = sx  - Ax;  constexpr auto dsAy  = sy  - Ay;  constexpr auto dsAz  = sz  - Az;
    constexpr auto dt1Ax = t1x - Ax;  constexpr auto dt1Ay = t1y - Ay;  constexpr auto dt1Az = t1z - Az;
    constexpr auto dt2Ax = t2x - Ax;  constexpr auto dt2Ay = t2y - Ay;  constexpr auto dt2Az = t2z - Az;
    constexpr auto sA2 = dsAx*dsAx + dsAy*dsAy + dsAz*dsAz;
    constexpr auto c1  = dt1Ax*dt1Ax + dt1Ay*dt1Ay + dt1Az*dt1Az - sA2;
    constexpr auto c2  = dt2Ax*dt2Ax + dt2Ay*dt2Ay + dt2Az*dt2Az - sA2;

    constexpr auto dw1 = ws - wt1;
    constexpr auto dw2 = ws - wt2;
    constexpr auto c1w = c1 + dw1;
    constexpr auto c2w = c2 + dw2;

    constexpr auto M1 = grp::det<decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nT1x),decltype(nT1y),decltype(nT1z),
                                 decltype(nT2x),decltype(nT2y),decltype(nT2z)>{};
    constexpr auto M2 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(nT1x),decltype(nT1y),decltype(nT1z),
                                 decltype(nT2x),decltype(nT2y),decltype(nT2z)>{};
    constexpr auto expr = c2w*M2 - c1w*M1;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp3_feasible_TSS_T(double Ax, double Ay, double Az,
                                        double Ex, double Ey, double Ez,
                                        double Cx, double Cy, double Cz,
                                        double Dx, double Dy, double Dz,
                                        double sx, double sy, double sz,
                                        double t1x, double t1y, double t1z,
                                        double t2x, double t2y, double t2z,
                                        double ws, double wt1, double wt2)
{
    return lp3_D_TSS(Ax,Ay,Az, Ex,Ey,Ez, Cx,Cy,Cz,
                     sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z)
         * power_lp3_feasible_TSS_T_factor2_impl::pred{}.apply(
                Ax,Ay,Az, Ex,Ey,Ez, Cx,Cy,Cz, Dx,Dy,Dz, sx,sy,sz,
                t1x,t1y,t1z, t2x,t2y,t2z, ws,wt1,wt2);
}


// ===========================================================================
// power_lp3_feasible_TSS_T_gen -- weighted analogue of lp3_feasible_TSS_T_gen
// (arbitrary query face T2, not sharing an edge with T1).  Same c1^w,c2^w
// substitution; the T2-vs-A offset term w = n_T2.(P2-A) (and its x2 scaling,
// int_const<2>) is weight-independent -- purely a T-face-vs-reference-point
// distance, untouched.
// Inputs: A(_1-_3), Q1(_4-_6), R1(_7-_9), P2(_10-_12), Q2(_13-_15), R2(_16-_18),
//         s(_19-_21), t1(_22-_24), t2(_25-_27), ws(_28), wt1(_29), wt2(_30)
// ===========================================================================
namespace power_lp3_feasible_TSS_T_gen_factor2_impl {
    constexpr auto Ax  = grp::_1;   constexpr auto Ay  = grp::_2;   constexpr auto Az  = grp::_3;
    constexpr auto Q1x = grp::_4;   constexpr auto Q1y = grp::_5;   constexpr auto Q1z = grp::_6;
    constexpr auto R1x = grp::_7;   constexpr auto R1y = grp::_8;   constexpr auto R1z = grp::_9;
    constexpr auto P2x = grp::_10;  constexpr auto P2y = grp::_11;  constexpr auto P2z = grp::_12;
    constexpr auto Q2x = grp::_13;  constexpr auto Q2y = grp::_14;  constexpr auto Q2z = grp::_15;
    constexpr auto R2x = grp::_16;  constexpr auto R2y = grp::_17;  constexpr auto R2z = grp::_18;
    constexpr auto sx  = grp::_19;  constexpr auto sy  = grp::_20;  constexpr auto sz  = grp::_21;
    constexpr auto t1x = grp::_22;  constexpr auto t1y = grp::_23;  constexpr auto t1z = grp::_24;
    constexpr auto t2x = grp::_25;  constexpr auto t2y = grp::_26;  constexpr auto t2z = grp::_27;
    constexpr auto ws  = grp::_28;  constexpr auto wt1 = grp::_29;  constexpr auto wt2 = grp::_30;

    constexpr auto dQ1Ax = Q1x-Ax; constexpr auto dQ1Ay = Q1y-Ay; constexpr auto dQ1Az = Q1z-Az;
    constexpr auto dR1Ax = R1x-Ax; constexpr auto dR1Ay = R1y-Ay; constexpr auto dR1Az = R1z-Az;
    constexpr auto nT1x = dQ1Ay*dR1Az - dQ1Az*dR1Ay;
    constexpr auto nT1y = dQ1Az*dR1Ax - dQ1Ax*dR1Az;
    constexpr auto nT1z = dQ1Ax*dR1Ay - dQ1Ay*dR1Ax;
    constexpr auto dQ2Px = Q2x-P2x; constexpr auto dQ2Py = Q2y-P2y; constexpr auto dQ2Pz = Q2z-P2z;
    constexpr auto dR2Px = R2x-P2x; constexpr auto dR2Py = R2y-P2y; constexpr auto dR2Pz = R2z-P2z;
    constexpr auto nT2x = dQ2Py*dR2Pz - dQ2Pz*dR2Py;
    constexpr auto nT2y = dQ2Pz*dR2Px - dQ2Px*dR2Pz;
    constexpr auto nT2z = dQ2Px*dR2Py - dQ2Py*dR2Px;
    constexpr auto d1x = t1x-sx; constexpr auto d1y = t1y-sy; constexpr auto d1z = t1z-sz;
    constexpr auto d2x = t2x-sx; constexpr auto d2y = t2y-sy; constexpr auto d2z = t2z-sz;
    constexpr auto dsAx = sx-Ax;   constexpr auto dsAy = sy-Ay;   constexpr auto dsAz = sz-Az;
    constexpr auto dt1Ax = t1x-Ax; constexpr auto dt1Ay = t1y-Ay; constexpr auto dt1Az = t1z-Az;
    constexpr auto dt2Ax = t2x-Ax; constexpr auto dt2Ay = t2y-Ay; constexpr auto dt2Az = t2z-Az;
    constexpr auto sA2 = dsAx*dsAx + dsAy*dsAy + dsAz*dsAz;
    constexpr auto c1  = dt1Ax*dt1Ax + dt1Ay*dt1Ay + dt1Az*dt1Az - sA2;
    constexpr auto c2  = dt2Ax*dt2Ax + dt2Ay*dt2Ay + dt2Az*dt2Az - sA2;

    constexpr auto dw1 = ws - wt1;
    constexpr auto dw2 = ws - wt2;
    constexpr auto c1w = c1 + dw1;
    constexpr auto c2w = c2 + dw2;

    // Weight-independent T2-vs-reference-point term (see lp3_feasible_TSS_T_gen).
    constexpr auto dP2Ax = P2x-Ax; constexpr auto dP2Ay = P2y-Ay; constexpr auto dP2Az = P2z-Az;
    constexpr auto w  = nT2x*dP2Ax + nT2y*dP2Ay + nT2z*dP2Az;
    constexpr auto w2 = int_const<2>{} * w;
    constexpr auto M1 = grp::det<decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nT1x),decltype(nT1y),decltype(nT1z),
                                 decltype(nT2x),decltype(nT2y),decltype(nT2z)>{};
    constexpr auto M2 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(nT1x),decltype(nT1y),decltype(nT1z),
                                 decltype(nT2x),decltype(nT2y),decltype(nT2z)>{};
    constexpr auto Dexpr = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                    decltype(d2x),decltype(d2y),decltype(d2z),
                                    decltype(nT1x),decltype(nT1y),decltype(nT1z)>{};
    constexpr auto expr = c2w*M2 - c1w*M1 + w2*Dexpr;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp3_feasible_TSS_T_gen(double Ax, double Ay, double Az,
                                            double Q1x, double Q1y, double Q1z,
                                            double R1x, double R1y, double R1z,
                                            double P2x, double P2y, double P2z,
                                            double Q2x, double Q2y, double Q2z,
                                            double R2x, double R2y, double R2z,
                                            double sx, double sy, double sz,
                                            double t1x, double t1y, double t1z,
                                            double t2x, double t2y, double t2z,
                                            double ws, double wt1, double wt2)
{
    return lp3_D_TSS(Ax,Ay,Az, Q1x,Q1y,Q1z, R1x,R1y,R1z,
                     sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z)
         * power_lp3_feasible_TSS_T_gen_factor2_impl::pred{}.apply(
                Ax,Ay,Az, Q1x,Q1y,Q1z, R1x,R1y,R1z,
                P2x,P2y,P2z, Q2x,Q2y,Q2z, R2x,R2y,R2z,
                sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z, ws,wt1,wt2);
}


// ===========================================================================
// power_lp3_feasible_SSS_T -- weighted analogue of lp3_feasible_SSS_T.
// feasible = sign(Delta_SSS) * sign(Gamma^w).  Delta_SSS = det[t1-s;t2-s;t3-s]
// is weight-invariant (a_S untouched) -- recomputed locally here, unweighted,
// duplicating lp3_feasible_SSS_T_factor1_impl exactly (that namespace lives
// in a different translation unit and isn't exported standalone, so it is
// re-derived here rather than reused -- everything downstream of it, the
// det<> exact-arithmetic pattern, is identical). Gamma is LINEAR in c1,c2,c3:
// direct substitution ci^w = ci + ws - w_ti, same M1,M2,M3 minors (pure
// orientation on t1-s,t2-s,t3-s,n_T) reused unweighted.
// Inputs: s(_1-_3), t1(_4-_6), t2(_7-_9), t3(_10-_12),
//         A(_13-_15), Q(_16-_18), R(_19-_21),
//         ws(_22), wt1(_23), wt2(_24), wt3(_25)
// ===========================================================================
namespace power_lp3_feasible_SSS_T_factor1_impl {
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
    // _13.._25 = A,Q,R,ws,wt1,wt2,wt3 unused (weight-invariant slope).

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

namespace power_lp3_feasible_SSS_T_factor2_impl {
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
    constexpr auto ws  = grp::_22;
    constexpr auto wt1 = grp::_23;
    constexpr auto wt2 = grp::_24;
    constexpr auto wt3 = grp::_25;

    constexpr auto dQAx = Qx - Ax;  constexpr auto dQAy = Qy - Ay;  constexpr auto dQAz = Qz - Az;
    constexpr auto dRAx = Rx - Ax;  constexpr auto dRAy = Ry - Ay;  constexpr auto dRAz = Rz - Az;
    constexpr auto nTx = dQAy*dRAz - dQAz*dRAy;
    constexpr auto nTy = dQAz*dRAx - dQAx*dRAz;
    constexpr auto nTz = dQAx*dRAy - dQAy*dRAx;

    constexpr auto d1x = t1x - sx;  constexpr auto d1y = t1y - sy;  constexpr auto d1z = t1z - sz;
    constexpr auto d2x = t2x - sx;  constexpr auto d2y = t2y - sy;  constexpr auto d2z = t2z - sz;
    constexpr auto d3x = t3x - sx;  constexpr auto d3y = t3y - sy;  constexpr auto d3z = t3z - sz;

    constexpr auto dsAx  = sx  - Ax;  constexpr auto dsAy  = sy  - Ay;  constexpr auto dsAz  = sz  - Az;
    constexpr auto dt1Ax = t1x - Ax;  constexpr auto dt1Ay = t1y - Ay;  constexpr auto dt1Az = t1z - Az;
    constexpr auto dt2Ax = t2x - Ax;  constexpr auto dt2Ay = t2y - Ay;  constexpr auto dt2Az = t2z - Az;
    constexpr auto dt3Ax = t3x - Ax;  constexpr auto dt3Ay = t3y - Ay;  constexpr auto dt3Az = t3z - Az;
    constexpr auto sA2 = dsAx*dsAx + dsAy*dsAy + dsAz*dsAz;
    constexpr auto c1  = dt1Ax*dt1Ax + dt1Ay*dt1Ay + dt1Az*dt1Az - sA2;
    constexpr auto c2  = dt2Ax*dt2Ax + dt2Ay*dt2Ay + dt2Az*dt2Az - sA2;
    constexpr auto c3  = dt3Ax*dt3Ax + dt3Ay*dt3Ay + dt3Az*dt3Az - sA2;

    constexpr auto dw1 = ws - wt1;
    constexpr auto dw2 = ws - wt2;
    constexpr auto dw3 = ws - wt3;
    constexpr auto c1w = c1 + dw1;
    constexpr auto c2w = c2 + dw2;
    constexpr auto c3w = c3 + dw3;

    constexpr auto M1 = grp::det<decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(d3x),decltype(d3y),decltype(d3z),
                                 decltype(nTx),decltype(nTy),decltype(nTz)>{};
    constexpr auto M2 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(d3x),decltype(d3y),decltype(d3z),
                                 decltype(nTx),decltype(nTy),decltype(nTz)>{};
    constexpr auto M3 = grp::det<decltype(d1x),decltype(d1y),decltype(d1z),
                                 decltype(d2x),decltype(d2y),decltype(d2z),
                                 decltype(nTx),decltype(nTy),decltype(nTz)>{};
    constexpr auto expr = c2w*M2 - c1w*M1 - c3w*M3;

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp3_feasible_SSS_T(double sx, double sy, double sz,
                                        double t1x, double t1y, double t1z,
                                        double t2x, double t2y, double t2z,
                                        double t3x, double t3y, double t3z,
                                        double Ax, double Ay, double Az,
                                        double Qx, double Qy, double Qz,
                                        double Rx, double Ry, double Rz,
                                        double ws, double wt1, double wt2, double wt3)
{
    return power_lp3_feasible_SSS_T_factor1_impl::pred{}.apply(
                sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z, t3x,t3y,t3z,
                Ax,Ay,Az, Qx,Qy,Qz, Rx,Ry,Rz, ws,wt1,wt2,wt3)
         * power_lp3_feasible_SSS_T_factor2_impl::pred{}.apply(
                sx,sy,sz, t1x,t1y,t1z, t2x,t2y,t2z, t3x,t3y,t3z,
                Ax,Ay,Az, Qx,Qy,Qz, Rx,Ry,Rz, ws,wt1,wt2,wt3);
}


// ===========================================================================
// power_lp2_feasible_S_S_T0 -- weighted analogue of lp2_feasible_S_S_T0.
// feasible = sign(a_j*b_k-a_k*b_j) * sign(b_k*c_j2^w - c_k2^w*b_j), where
// a_m=(tm-s).(B-A), b_m=(tm-s).(C-A) are pure direction dot-products
// (weight-invariant, a_S untouched) and reused unweighted via lp2_det2 for
// factor1; c_m^w = c_m + (ws - w_tm), c_m = |tm-A|^2-|s-A|^2 is the LINEAR
// leaf substitution (each c appears to the first power, multiplied only by
// the weight-invariant b, so direct substitution is exact -- no cross terms,
// same class as power_lp3_feasible_TSS_S's Gamma).
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9), s(_10-_12), tj(_13-_15), tk(_16-_18),
//         ws(_19), wtj(_20), wtk(_21)
// ===========================================================================
namespace power_lp2_feasible_S_S_T0_factor2_impl {
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
    constexpr auto ws  = grp::_19;
    constexpr auto wtj = grp::_20;
    constexpr auto wtk = grp::_21;

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

    constexpr auto dwj = ws - wtj;
    constexpr auto dwk = ws - wtk;
    constexpr auto c_jw = c_j + dwj;
    constexpr auto c_kw = c_k + dwk;

    // -(b_j*c_k2^w - c_j2^w*b_k)
    using det2_t = grp::det<decltype(b_k), decltype(c_kw),
                            decltype(b_j), decltype(c_jw)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp2_feasible_S_S_T0(double Ax, double Ay, double Az,
                                         double Bx, double By, double Bz,
                                         double Cx, double Cy, double Cz,
                                         double sx, double sy, double sz,
                                         double tjx, double tjy, double tjz,
                                         double tkx, double tky, double tkz,
                                         double ws, double wtj, double wtk)
{
    return lp2_det2(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                    sx,sy,sz, tjx,tjy,tjz, tkx,tky,tkz)
         * power_lp2_feasible_S_S_T0_factor2_impl::pred{}.apply(
                Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz, sx,sy,sz,
                tjx,tjy,tjz, tkx,tky,tkz, ws,wtj,wtk);
}


// ===========================================================================
// power_lp2_feasible_S_S_T1 -- weighted analogue of lp2_feasible_S_S_T1.
// feasible = sign(a_j*b_k-a_k*b_j) * sign(c_j2^w*a_k - a_j*c_k2^w).  a_m
// (weight-invariant) reused unweighted; c_m^w = c_m + (ws-w_tm) direct
// substitution, same linear-in-c class as S_S_T0.
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9, unused), s(_10-_12), tj(_13-_15),
//         tk(_16-_18), ws(_19), wtj(_20), wtk(_21)
// ===========================================================================
namespace power_lp2_feasible_S_S_T1_factor2_impl {
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
    constexpr auto ws  = grp::_19;
    constexpr auto wtj = grp::_20;
    constexpr auto wtk = grp::_21;

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

    constexpr auto dwj = ws - wtj;
    constexpr auto dwk = ws - wtk;
    constexpr auto c_jw = c_j + dwj;
    constexpr auto c_kw = c_k + dwk;

    // -(c_j2^w*a_k - a_j*c_k2^w)
    using det2_t = grp::det<decltype(c_kw), decltype(a_k),
                            decltype(c_jw), decltype(a_j)>;
    constexpr auto expr = det2_t{};

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp2_feasible_S_S_T1(double Ax, double Ay, double Az,
                                         double Bx, double By, double Bz,
                                         double Cx, double Cy, double Cz,
                                         double sx, double sy, double sz,
                                         double tjx, double tjy, double tjz,
                                         double tkx, double tky, double tkz,
                                         double ws, double wtj, double wtk)
{
    return lp2_det2(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                    sx,sy,sz, tjx,tjy,tjz, tkx,tky,tkz)
         * power_lp2_feasible_S_S_T1_factor2_impl::pred{}.apply(
                Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz, sx,sy,sz,
                tjx,tjy,tjz, tkx,tky,tkz, ws,wtj,wtk);
}


// ===========================================================================
// power_lp2_feasible_S_S_T2 -- weighted analogue of lp2_feasible_S_S_T2.
// feasible = sign(a_j*b_k-a_k*b_j) * sign(c_j2^w*b_k - b_j*c_k2^w
//   + a_j*c_k2^w - c_j2^w*a_k - 2*a_j*b_k + 2*b_j*a_k).  a_m,b_m
// (weight-invariant) reused unweighted; c_m^w = c_m + (ws-w_tm) direct
// substitution -- expr is linear in c_j, c_k separately (never c_j*c_k), so
// substituting the leaf is exact, same class as S_S_T0/T1.
// Inputs: A(_1-_3), B(_4-_6), C(_7-_9), s(_10-_12), tj(_13-_15), tk(_16-_18),
//         ws(_19), wtj(_20), wtk(_21)
// ===========================================================================
namespace power_lp2_feasible_S_S_T2_factor2_impl {
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
    constexpr auto ws  = grp::_19;
    constexpr auto wtj = grp::_20;
    constexpr auto wtk = grp::_21;

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

    constexpr auto dwj = ws - wtj;
    constexpr auto dwk = ws - wtk;
    constexpr auto c_jw = c_j + dwj;
    constexpr auto c_kw = c_k + dwk;

    // -(c_j2^w*b_k - b_j*c_k2^w + a_j*c_k2^w - c_j2^w*a_k - 2*a_j*b_k + 2*b_j*a_k)
    constexpr auto expr = b_j*c_kw - c_jw*b_k - a_j*c_kw + c_jw*a_k
                        + (a_j*b_k + a_j*b_k) - (b_j*a_k + b_j*a_k);

    using semi_static = grp::forward_error_semi_static<expr, double, grp::robust_rules<true>>;
    using exact       = grp::stage_d<expr, double>;
    using pred        = grp::staged_predicate<semi_static, exact>;
}

extern "C" int power_lp2_feasible_S_S_T2(double Ax, double Ay, double Az,
                                         double Bx, double By, double Bz,
                                         double Cx, double Cy, double Cz,
                                         double sx, double sy, double sz,
                                         double tjx, double tjy, double tjz,
                                         double tkx, double tky, double tkz,
                                         double ws, double wtj, double wtk)
{
    return lp2_det2(Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz,
                    sx,sy,sz, tjx,tjy,tjz, tkx,tky,tkz)
         * power_lp2_feasible_S_S_T2_factor2_impl::pred{}.apply(
                Ax,Ay,Az, Bx,By,Bz, Cx,Cy,Cz, sx,sy,sz,
                tjx,tjy,tjz, tkx,tky,tkz, ws,wtj,wtk);
}
