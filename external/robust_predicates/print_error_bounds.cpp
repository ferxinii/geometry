/* Made with help from Claude and taking inspiration from 
 * util/render_error_bound.hpp in the paper's repository */

#include <array>
#include <string>
#include <sstream>
#include <limits>
#include <iomanip>
#include <iostream>
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp"
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/forward_error_bound.hpp"
#include "robust_predicates.cpp"


// ---------------------------------------------------------------------------
// ----------- GENERIC PRINTING FUNCTIONS ------------------------------------
// ---------------------------------------------------------------------------

using namespace boost::geometry::detail::generic_robust_predicates;

template <typename T>
constexpr T constexpr_abs(T x) { return x < 0 ? -x : x; }

template <typename Expression>
constexpr int underflow_guard_count()
{   
    /* detect underflow_guard_constant nodes (to print k * uN) */

    if constexpr (Expression::is_leaf && Expression::argn == 0)
    {
        using ct = typename Expression::value_type;
        constexpr ct uN = std::numeric_limits<ct>::min();
        constexpr ct val = Expression::value;
        if (val <= 0) return 0;
        // check if val is an integer multiple of uN
        constexpr int k = (int)(val / uN + 0.5);
        if (constexpr_abs(val - k * uN) < uN * 0.01)  // 10% margin?
            return k;
    }
    return 0;
}

template <typename Expression, std::size_t N>
std::string to_latex(const std::array<const char*, N>& dict)
{
    /* symbolic expression printer
       takes a dict mapping argument indices to variable names */

    std::stringstream out;

    if constexpr (Expression::is_leaf) {
        if constexpr (Expression::argn == 0) {
            constexpr int k = underflow_guard_count<Expression>();
            if (k > 0) {
                if (k == 1) out << "u_N";
                else out << k << "u_N";
            } else {
                out << std::setprecision(std::numeric_limits<double>::digits10)
                    << Expression::value;
            }
        } else {
            if constexpr (Expression::argn < N) out << dict[Expression::argn];
            else out << "_" << Expression::argn;
        }
    }
    else if constexpr (Expression::operator_type == operator_types::abs) {
        // skip abs around a square — x^2 >= 0 always so |x^2| = x^2
        if constexpr (Expression::child::operator_type == operator_types::product &&
                      std::is_same_v<typename Expression::child::left,
                                     typename Expression::child::right>)
            out << to_latex<typename Expression::child>(dict);
        else
            out << "\\left|" << to_latex<typename Expression::child>(dict) << "\\right|";
    }

    else if constexpr (Expression::operator_type == operator_types::sum) {
        if constexpr (Expression::left::operator_type == operator_types::difference)
            out << "\\left(" << to_latex<typename Expression::left>(dict) << "\\right)";
        else out << to_latex<typename Expression::left>(dict);

        out << " \\oplus ";

        if constexpr (Expression::right::operator_type == operator_types::difference)
            out << "\\left(" << to_latex<typename Expression::right>(dict) << "\\right)";
        else out << to_latex<typename Expression::right>(dict);
    }

    else if constexpr (Expression::operator_type == operator_types::difference) {
        if constexpr (Expression::left::operator_type == operator_types::difference)
            out << "\\left(" << to_latex<typename Expression::left>(dict) << "\\right)";
        else out << to_latex<typename Expression::left>(dict);

        out << " \\ominus ";

        if constexpr (   Expression::right::operator_type == operator_types::difference
                      || Expression::right::operator_type == operator_types::sum)
            out << "\\left(" << to_latex<typename Expression::right>(dict) << "\\right)";
        else out << to_latex<typename Expression::right>(dict);
    }

    else if constexpr (Expression::operator_type == operator_types::product) {
        constexpr int power = right_power<typename Expression::left,
                                          typename Expression::right>();
        if constexpr (power == -1)
        {
            if constexpr (   Expression::left::operator_type == operator_types::sum
                          || Expression::left::operator_type == operator_types::difference)
                out << "\\left(" << to_latex<typename Expression::left>(dict) << "\\right)";
            else out << to_latex<typename Expression::left>(dict);

            out << " \\odot ";

            if constexpr (   Expression::right::operator_type == operator_types::sum
                          || Expression::right::operator_type == operator_types::difference)
                out << "\\left(" << to_latex<typename Expression::right>(dict) << "\\right)";
            else
                out << to_latex<typename Expression::right>(dict);
        } else {
            if constexpr (   Expression::left::operator_type == operator_types::sum
                          || Expression::left::operator_type == operator_types::difference)
                out << "\\left(" << to_latex<typename Expression::left>(dict)
                    << "\\right)^{" << (power + 1) << "}";
            else
                out << to_latex<typename Expression::left>(dict)
                    << "^{" << (power + 1) << "}";
        }
    }

    return out.str();
}

template <typename Expression, std::size_t N>
void print_error_bound(const char* name,
                       const std::array<const char*, N>& dict)
{   
    /* print full error bound: coefficient * magnitude
       with epsilon coefficients shown symbolically */

    using eb = forward_error_bound<Expression, double, robust_rules<true>>;
    using magnitude = typename eb::error_bound::magnitude;
    constexpr auto a = eb::error_bound::a;

    std::cout << "=== " << name << " ===\n";

    // print magnitude expression symbolically
    std::cout << "E = \\left(";
    if (a[0] != 0) std::cout << a[0] << "\\varepsilon";
    if (a[1] != 0) {
        if (a[1] > 0) std::cout << " + " << a[1] << "\\varepsilon^2";
        else          std::cout << " - " << -a[1] << "\\varepsilon^2";
    }
    if (a[2] != 0) {
        if (a[2] > 0) std::cout << " + " << a[2] << "\\varepsilon^3";
        else          std::cout << " - " << -a[2] << "\\varepsilon^3";
    }
    std::cout << ", " << to_latex<magnitude>(dict) << "\\right)\n\n";
}




// ---------------------------------------------------------------------------
// --------- PRINT IMPLEMENTED PREDICATES ------------------------------------
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// dictionaries
// ---------------------------------------------------------------------------
constexpr std::array<const char*, 7> powertest1d_dict {
    "",
    "x_a", "w_a", "x_b", "w_b", "x_c", "w_c"
};

constexpr std::array<const char*, 13> powertest2d_dict {
    "",
    "a_x", "a_y", "w_a",
    "b_x", "b_y", "w_b",
    "c_x", "c_y", "w_c",
    "d_x", "d_y", "w_d"
};

constexpr std::array<const char*, 21> powertest3d_dict {
    "",
    "a_x", "a_y", "a_z", "w_a",
    "b_x", "b_y", "b_z", "w_b",
    "c_x", "c_y", "c_z", "w_c",
    "d_x", "d_y", "d_z", "w_d",
    "e_x", "e_y", "e_z", "w_e"
};

constexpr std::array<const char*, 7> orient2d_dict {
    "",
    "a_x", "a_y", "b_x", "b_y", "c_x", "c_y"
};

constexpr std::array<const char*, 13> orient3d_dict {
    "",
    "a_x", "a_y", "a_z",
    "b_x", "b_y", "b_z",
    "c_x", "c_y", "c_z",
    "d_x", "d_y", "d_z"
};

constexpr std::array<const char*, 9> incircle_dict {
    "",
    "a_x", "a_y", "b_x", "b_y", "c_x", "c_y", "d_x", "d_y"
};

constexpr std::array<const char*, 16> insphere_dict {
    "",
    "a_x", "a_y", "a_z",
    "b_x", "b_y", "b_z",
    "c_x", "c_y", "c_z",
    "d_x", "d_y", "d_z",
    "e_x", "e_y", "e_z"
};





// ---------------------------------------------------------------------------
// add a raw evaluator to print_expressions.cpp
#include <array>
#include "boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_eval.hpp"

// compare insphere vs powertest3d with zero weights
namespace compare_impl {
    // reproduce insphere matrix using powertest3d arguments but zero weights
    constexpr auto ax=grp::_1;
    constexpr auto ay=grp::_2;
    constexpr auto az=grp::_3;
    constexpr auto bx=grp::_5;
    constexpr auto by=grp::_6;
    constexpr auto bz=grp::_7;
    constexpr auto cx=grp::_9;
    constexpr auto cy=grp::_10;
    constexpr auto cz=grp::_11;
    constexpr auto dx=grp::_13;
    constexpr auto dy=grp::_14;
    constexpr auto dz=grp::_15;
    constexpr auto ex=grp::_17;
    constexpr auto ey=grp::_18;
    constexpr auto ez=grp::_19;

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

    using expr_t = grp::det <
        decltype(dax), decltype(day), decltype(daz), decltype(la),
        decltype(dbx), decltype(dby), decltype(dbz), decltype(lb),
        decltype(dcx), decltype(dcy), decltype(dcz), decltype(lc),
        decltype(ddx), decltype(ddy), decltype(ddz), decltype(ld)
    >;
    constexpr auto expr = expr_t{};
}

template <typename Expr, std::size_t N>
double eval_expr(const std::array<double, N>& args) {
    return grp::evaluate_expression<Expr>(args);
}

int main()
{

    // test points
    std::array<double, 15> insphere_args {{
        1,0,0,   // pa
       -1,0,0,   // pb
        0,1,0,   // pc
        0,0,1,   // pd
        0,0,0    // pe_in
    }};

    std::array<double, 20> powertest_args {{
        1,0,0,0,   // pa, wa=0
       -1,0,0,0,   // pb, wb=0
        0,1,0,0,   // pc, wc=0
        0,0,1,0,   // pd, wd=0
        0,0,0,0    // pe, we=0
    }};

    double insphere_val = eval_expr<decltype(insphere_impl::expr)>(insphere_args);
    double powertest_val = eval_expr<decltype(compare_impl::expr)>(powertest_args);

    std::cout << "insphere value:  " << insphere_val  << "\n";
    std::cout << "powertest value: " << powertest_val << "\n";
    std::cout << "ratio: " << insphere_val / powertest_val << "\n";


    std::cout << "insphere expression:\n";
    std::cout << print_expression<decltype(grp::insphere{})>() << "\n";
    std::cout << "insphere_weighted expression:\n";
    std::cout << print_expression<decltype(powertest3d_impl::expr)>() << "\n";

    print_error_bound<decltype(orient2d_impl::expr)>
        ("orient2d", orient2d_dict);

    print_error_bound<decltype(orient3d_impl::expr)>
        ("orient3d", orient3d_dict);

    print_error_bound<decltype(incircle_impl::expr)>
        ("incircle", incircle_dict);

    print_error_bound<decltype(insphere_impl::expr)>
        ("insphere", insphere_dict);

    print_error_bound<decltype(powertest1d_impl::expr)>
        ("powertest1d", powertest1d_dict);

    print_error_bound<decltype(powertest2d_impl::expr)>
        ("powertest2d", powertest2d_dict);

    print_error_bound<decltype(powertest3d_impl::expr)>
        ("powertest3d", powertest3d_dict);

    return 0;
}

