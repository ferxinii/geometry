#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FORWARD_ERROR_BOUND_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FORWARD_ERROR_BOUND_HPP

#include <array>
#include <limits>
#include <type_traits>

#include <boost/mp11/utility.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/simple_orient2d.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/stage_a_error_bound.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/semi_static_filter.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename Expression, typename CalculationType, typename Rules>
struct forward_error_bound;


template <typename CalculationType, int c = 1>
struct underflow_guard_constant : public static_constant_interface<CalculationType>
{
    // 2 * u_N
    static constexpr CalculationType value =
        ( c * std::numeric_limits<CalculationType>::min() )
#ifdef BOOST_GEOMETRY_GRP_AVOID_DENORM
        / ( std::numeric_limits<CalculationType>::epsilon() / 2.0 )
#endif
        ;
    static constexpr bool non_negative = true;
    static constexpr int count = c;
};

template <typename Exp>
constexpr int count_underflow_guards()
{
    if constexpr ( Exp::operator_type == operator_types::sum )
    {
        using r = typename Exp::right;
        if constexpr ( requires { r::count; typename r::value_type; } )
        {
            if constexpr (std::is_same_v<r, underflow_guard_constant<typename r::value_type, r::count>>)
                return r::count;
        }
    }
    return 0;
}

template <typename Exp>
constexpr auto strip_underflow_guards()
{
    if constexpr ( Exp::operator_type == operator_types::sum )
    {
        using r = typename Exp::right;
        if constexpr ( requires { r::count; typename r::value_type; } )
        {
            if constexpr (std::is_same_v<r, underflow_guard_constant<typename r::value_type, r::count>>)
            {
                return typename Exp::left{};
            }
        }
    }
    else
    {
        return Exp{};
    }
}

template
<
    bool UnderflowProtection = true
>
struct ozaki_simple_fp_lemma_31
{
    template <typename Expression, typename, typename Rules>
    static constexpr bool applicable()
    {
        if constexpr( Expression::operator_type == operator_types::product )
        {
            if constexpr(    (   Expression::left::operator_type == operator_types::sum
                              || Expression::left::operator_type == operator_types::difference)
                          && (   Expression::right::operator_type == operator_types::sum
                              || Expression::right::operator_type == operator_types::difference) )
            {
                return Expression::left::left::is_leaf && Expression::left::right::is_leaf
                    && Expression::right::left::is_leaf && Expression::right::right::is_leaf;
            }
        }
        return false;
    }

    template <typename Expression, typename CalculationType, typename>
    struct error_bound
    {
        using magnitude = 
            mp11::mp_if_c
                <
                    UnderflowProtection,
                    sum
                        <
                            abs<Expression>,
                            underflow_guard_constant<CalculationType>
                        >,
                    abs<Expression>
                >;

        static constexpr std::array<long, 3> a {3, -(phi<CalculationType> - 14), 0};
    };
};

template <bool UnderflowProtection = true>
struct inexact_leaves
{
    template <typename Expression, typename, typename>
    static constexpr bool applicable()
    {
        return Expression::is_leaf;
    }

    template <typename Expression, typename CalculationType, typename>
    struct error_bound
    {
    private:
        template <typename Constant>
        struct abs_constant : public Constant
        {
            static constexpr auto value = std::abs(Constant::value);
        };

        static constexpr auto magnitude_impl()
        {
            if constexpr (Expression::argn == 0)
            {
                return abs_constant<Expression>{};
            }
            else
            {
                return mp11::mp_if_c
                        <
                            UnderflowProtection,
                            sum
                                <
                                    abs<Expression>,
                                    underflow_guard_constant<CalculationType>
                                >,
                            abs<Expression>
                        >{};
            }
        }
    public:
        using magnitude = decltype(magnitude_impl());
        static constexpr std::array<long, 3> a {1, 0, 0};
    };
};

struct exact_leaves
{
    template <typename Expression, typename, typename>
    static constexpr bool applicable()
    {
        return Expression::is_leaf;
    }

    template <typename Expression, typename, typename>
    struct error_bound
    {
    private:
        template <typename Constant>
        struct abs_constant : public Constant
        {
            static constexpr auto value = std::abs(Constant::value);
        };

        static constexpr auto magnitude_impl()
        {
            if constexpr (Expression::argn == 0)
            {
                return abs_constant<Expression>{};
            }
            else
            {
                return abs<Expression>{};
            }
        }
    public:
        using magnitude = decltype(magnitude_impl());
        static constexpr std::array<long, 3> a {0, 0, 0};
    };
};

struct exact_leaves_sumdiff
{
    template <typename Expression, typename, typename Rules>
    static constexpr bool applicable()
    {
        if constexpr(    Expression::operator_type == operator_types::sum 
                      || Expression::operator_type == operator_types::difference )
        {
            return Expression::left::is_leaf && Expression::right::is_leaf;
        }
        return false;
    }

    template <typename Expression, typename, typename>
    struct error_bound
    {
        using magnitude = abs<Expression>;
        static constexpr std::array<long, 3> a {1, 0, 0};
    };
};

template
<
    bool UnderflowProtection = true
>
struct exact_leaves_product
{
    template <typename Expression, typename, typename Rules>
    static constexpr bool applicable()
    {
        if constexpr( Expression::operator_type == operator_types::product )
        {
            return Expression::left::is_leaf && Expression::right::is_leaf;
        }
        return false;
    }

    template <typename Expression, typename CalculationType, typename>
    struct error_bound
    {
        using magnitude =
            mp11::mp_if_c
                <
                    UnderflowProtection,
                    sum
                        <
                            abs<Expression>,
                            underflow_guard_constant<CalculationType>
                        >,
                    abs<Expression>
                >;
        static constexpr std::array<long, 3> a {1, 0, 0};
    };
};

template <typename LEB, typename REB, typename CalculationType>
constexpr auto collapse_underflow_guards()
{
    constexpr int ugs =   count_underflow_guards<typename LEB::magnitude>()
                        + count_underflow_guards<typename REB::magnitude>();
    if constexpr (ugs == 0)
    {
        return sum<typename LEB::magnitude, typename REB::magnitude>{};
    }
    else
    {
        return sum
            <
                sum
                    <
                        decltype(strip_underflow_guards<typename LEB::magnitude>()),
                        decltype(strip_underflow_guards<typename REB::magnitude>())
                    >,
                underflow_guard_constant<CalculationType, ugs>
            >{};
    }
}

struct inexacts_sumdiff
{
    template <typename Expression, typename, typename>
    static constexpr bool applicable()
    {
        return    Expression::operator_type == operator_types::sum
               || Expression::operator_type == operator_types::difference;
    }

    template <typename Expression, typename CalculationType, typename Rules>
    struct error_bound
    {
    private:
        using leb = typename forward_error_bound<typename Expression::left, CalculationType, Rules>::error_bound;
        using reb = typename forward_error_bound<typename Expression::right, CalculationType, Rules>::error_bound;
        static constexpr std::array<long, 3> max_a = coeff_max(leb::a, reb::a);
    public:
        using magnitude = decltype(collapse_underflow_guards<leb, reb, CalculationType>());
        static constexpr std::array<long, 3> a = coeff_inc_first(coeff_mult_by_1_plus_eps(max_a));
    };
};

template
<
    bool UnderflowProtection = true
>
struct pow2_product
{
private:
    template <typename Expression>
    static constexpr auto pow2_exp()
    {
        if constexpr (Expression::is_leaf)
        {
            if constexpr (Expression::argn == 0)
            {
                auto val = Expression::value;
                if (val == 1)
                    return std::pair {true, 0};
                if (val < 0)
                    val = -val;
                int exponent = 0;
                if (val < 1)
                {
                    while (val < 1)
                    {
                        val *= 2;
                        --exponent;
                    }
                    return std::pair {val == 1, val == 1 ? exponent : 0};
                }
                if (val > 1)
                {
                    while (val > 1)
                    {
                        exponent++;
                        if( (val / 2) * 2 != val )
                            return std::pair {false, 0};
                        val /= 2;
                    }
                    return std::pair {val == 1, val == 1 ? exponent : 0};
                }
            }
        }
        return std::pair {false, 0};
    }
public:
    template <typename Expression, typename, typename>
    static constexpr bool applicable()
    {
        if constexpr (Expression::operator_type == operator_types::product)
        {
            constexpr auto lexp = pow2_exp<typename Expression::left>();
            constexpr auto rexp = pow2_exp<typename Expression::right>();
            return    (   lexp.first
                       && (!UnderflowProtection || lexp.second >= 0 ))
                   || (   rexp.first
                       && (!UnderflowProtection || rexp.second >= 0 ));
        }
        return false;
    }

    template <typename Expression, typename CalculationType, typename Rules>
    struct error_bound
    {
    private:
        using leb = typename forward_error_bound
            <
                typename Expression::left,
                CalculationType,
                Rules
            >::error_bound;
        using reb = typename forward_error_bound
            <
                typename Expression::right,
                CalculationType,
                Rules
            >::error_bound;
    public:
        using magnitude =
                product<typename leb::magnitude, typename reb::magnitude>;
        static constexpr std::array<long, 3> a = coeff_max(leb::a, reb::a);
    };
};

template
<
    bool UnderflowProtection = true
>
struct inexacts_product
{
    template <typename Expression, typename, typename>
    static constexpr bool applicable()
    {
        return Expression::operator_type == operator_types::product;
    }

    template <typename Expression, typename CalculationType, typename Rules>
    struct error_bound
    {
    private:
        using leb = typename forward_error_bound
            <
                typename Expression::left,
                CalculationType,
                Rules
            >::error_bound;
        using reb = typename forward_error_bound
            <
                typename Expression::right,
                CalculationType,
                Rules
            >::error_bound;
        static constexpr std::array<long, 3> a_prod = coeff_product(leb::a, reb::a);
    public:
        using magnitude = mp11::mp_if_c
            <
                UnderflowProtection,
                sum
                    <
                        product
                            <
                                typename leb::magnitude,
                                typename reb::magnitude
                            >,
                        underflow_guard_constant<CalculationType>
                    >,
                product<typename leb::magnitude, typename reb::magnitude>
            >;
        static constexpr std::array<long, 3> a = coeff_inc_first(coeff_mult_by_1_plus_eps(a_prod));
    };
};

template <typename Expression, typename CalculationType, typename Rules>
struct applicable_to
{
    template <typename Rule>
    using fn = mp11::mp_bool<Rule::template applicable<Expression, CalculationType, Rules>()>;
};

template <typename Expression, typename CalculationType, typename Rules>
struct forward_error_bound
{
    using rule = mp11::mp_at
        <
            Rules,
            mp11::mp_find_if_q
                <
                    Rules,
                    applicable_to<Expression, CalculationType, Rules>
                >
        >;

    using error_bound = typename rule::template error_bound
        <
            Expression,
            CalculationType,
            Rules
        >;
};

template <typename Expression, typename CalculationType, typename Rules>
struct forward_error_condition_sumdiff
{
private:
    using leb = typename forward_error_bound<typename Expression::left, CalculationType, Rules>::error_bound;
    using reb = typename forward_error_bound<typename Expression::right, CalculationType, Rules>::error_bound;
    static constexpr auto max_a = coeff_max(leb::a, reb::a);
    static constexpr std::array<long, 3> a =
        coeff_mult_by_1_plus_eps(coeff_mult_by_1_plus_eps(coeff_div_by_1_minus_eps(max_a)));
    static constexpr long c = round_to_next_2_pow<a[0]>();
    static constexpr long eps_square_coeff =
        c == 0 ? 0 : (a[2] > 0 ? c * ((a[1] + 1) / c + 1) : c * (a[1] / c + 1));
public:
    using magnitude = decltype(collapse_underflow_guards<leb, reb, CalculationType>());//sum<typename leb::magnitude, typename reb::magnitude>;
    static constexpr std::array<long, 2> coefficients {a[0], eps_square_coeff};
};

template
<
    typename Expression,
    typename CalculationType,
    typename Rules
>
struct forward_error_bound_expression_impl
{
private:
    using ct = CalculationType;
    using fe_cond = forward_error_condition_sumdiff<Expression, CalculationType, Rules>;
    static constexpr ct eps = std::numeric_limits<ct>::epsilon() / 2.0;
    struct constant : public static_constant_interface<ct>
    {
        static constexpr ct value =
              fe_cond::coefficients[0] * eps
            + fe_cond::coefficients[1] * eps * eps;
        static constexpr bool non_negative = true;
    };
public:
    using type = product<constant, typename fe_cond::magnitude>;
};

template <bool UnderflowProtection>
using robust_rules =
    mp11::mp_list
        <
            exact_leaves,
            exact_leaves_sumdiff,
            pow2_product<UnderflowProtection>,
            exact_leaves_product<UnderflowProtection>,
            inexacts_sumdiff,
            ozaki_simple_fp_lemma_31<UnderflowProtection>,
            inexacts_product<UnderflowProtection>
        >;

template <bool UnderflowProtection>
using rounded_input_rules =
    mp11::mp_list
        <
            inexact_leaves<UnderflowProtection>,
            pow2_product<UnderflowProtection>,
            inexacts_sumdiff,
            inexacts_product<UnderflowProtection>
        >;

template
<
    auto Expression,
    typename CalculationType,
    typename Rules = robust_rules<true>
>
using forward_error_bound_expression =
    typename forward_error_bound_expression_impl<decltype(Expression), CalculationType, Rules>::type;

template
<
    auto Expression,
    typename CalculationType = double,
    typename Rules = robust_rules<true>
>
using forward_error_semi_static =
    semi_static_filter
        <
            Expression,
            forward_error_bound_expression<Expression, CalculationType, Rules>{},
            CalculationType
        >;

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_FORWARD_ERROR_BOUND_HPP
