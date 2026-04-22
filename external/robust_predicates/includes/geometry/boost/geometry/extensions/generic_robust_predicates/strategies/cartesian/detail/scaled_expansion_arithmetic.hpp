// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2021 Tinko Bartels, Berlin, Germany.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SCALED_EXPANSION_ARITHMETIC_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SCALED_EXPANSION_ARITHMETIC_HPP

//TODO remove
#include <iostream>

#include <cassert>
#include <cmath>
#include <utility>
#include <limits>
#include <algorithm>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expansion_arithmetic.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename FPT>
constexpr FPT ce_exp2(auto exp)
{
    if(exp == 0) return FPT(1);
    const FPT factor = exp > 0 ? FPT(2) : FPT(1./2.);
    exp *= (exp > 0 ? 1 : -1);
    FPT a(1.);
    for(int i = 0; i < exp; ++i)
        a *= factor;
    return a;
}

template <typename FPT>
constexpr int max_exp = (std::numeric_limits<FPT>::max_exponent - 1) / 2;

template <typename FPT, int degree = 1>
constexpr int min_exp =
    (std::numeric_limits<FPT>::min_exponent + 1 + std::numeric_limits<FPT>::digits * degree) / 2 + 1;

template <typename FPT>
constexpr FPT scaling_upper = ce_exp2<FPT>(max_exp<FPT>);

template <typename FPT, int degree = 1>
constexpr FPT scaling_lower = ce_exp2<FPT>(min_exp<FPT, degree>);

template <typename FPT>
constexpr FPT scaling_upper_tail = scaling_upper<FPT> * std::numeric_limits<FPT>::epsilon();

template <typename FPT, int degree = 1>
constexpr FPT scaling_factor = ce_exp2<FPT>(max_exp<FPT> - min_exp<FPT, degree>);

constexpr int exp_end = std::numeric_limits<int>::max();

namespace debug_scaled_expansion
{

template
<
    int degree = 1,
    typename ExpIter
>
inline bool nonoverlapping(ExpIter e_exps_begin,
                           ExpIter e_exps_end)
{
    using FPT = std::remove_cvref_t<decltype(*e_exps_begin->first)>;
    std::vector<FPT> comps;
    for(auto it = e_exps_begin; it->second != exp_end; ++it)
    {
        if(it->first >= std::next(it)->first || it->second >= std::next(it)->second)
            return false;
        comps.clear();
        auto e_exp = it->second;
        for (auto it2 = e_exps_begin; it2->second != exp_end; ++it2)
        {
            auto exp = it2->second - e_exp;
            auto sf = 1.;
            if (exp < -1) sf = 0;
            else if (exp == -1) sf /= scaling_factor<FPT, degree>;
            else if (exp ==  0) sf = 1.;
            else if (exp ==  1) sf *= scaling_factor<FPT, degree>;
            else if (exp >   1) sf = 0;
            auto e_it = it2->first;
            auto next = std::next(it2)->first;
            for (; e_it != next; ++e_it)
            {
                if (*e_it > scaling_upper<FPT>) return false;
                if (*e_it <= scaling_lower<FPT, degree>) return false;
                comps.push_back( (*e_it) * sf );
            }
        }
        if ( !debug_expansion::expansion_nonoverlapping(comps.cbegin(), comps.cend()) )
            return false;
    }
    return true;
}

template
<
    typename ExpIter
>
inline void print(ExpIter e_exps_begin, ExpIter e_exps_end)
{
    auto e_begin = e_exps_begin->first;
    std::cout << "Printing expansion: \n";
    auto e_exps_it = e_exps_begin;
    for (; e_exps_it->second != exp_end; ++e_exps_it)
    {
        std::cout << "SF ^ " << e_exps_it->second << ": \t";
        for (auto e_it = e_exps_it->first; e_it != std::next(e_exps_it)->first; ++e_it)
        {
            std::cout << *e_it << ", ";
        }
        std::cout << "\n";
    }
    std::cout << "total length: " << std::distance(e_begin, e_exps_it->first) << "\n\n";
}

}

template
<
    int degree = 1,
    typename FPT
>
constexpr auto scale_number(FPT a)
{
    if (std::abs(a) != 0 && std::abs(a) <= scaling_lower<FPT, degree>)
        return std::pair {a * scaling_factor<FPT, degree>, -1};
    if (std::abs(a) > scaling_upper<FPT>)
        return std::pair {a / scaling_factor<FPT, degree>, 1};
    return std::pair {a, 0};
}

template
<
    bool NegateE = false,
    bool NegateB = false,
    typename FPT,
    typename ExpInIter,
    typename ExpOutIter
>
constexpr ExpOutIter scaled_merge_expansion(ExpInIter e_exps_begin,
                                            ExpInIter,
                                            FPT b,
                                            ExpOutIter h_exps_begin,
                                            ExpOutIter)
{
    auto h_begin = h_exps_begin->first;
    auto h_it = h_begin;
    if constexpr(NegateB != NegateE)
        b = -b;
    const auto [b_n, b_exp] = scale_number(b);
    auto e_exps_it = e_exps_begin;
    auto h_exps_it = h_exps_begin;
    while (e_exps_it->second < b_exp && e_exps_it->second != exp_end)
    {
        //std::cout << "adding lower segment exp " << e_exps_it->second << "\n";
        *h_exps_it++ = std::pair{h_it, e_exps_it->second};
        auto next = std::next(e_exps_it)->first;
        h_it = std::copy(e_exps_it->first, next, h_it);
        ++e_exps_it;
    }
    *h_exps_it++ = std::pair{h_it, b_exp};
    //std::cout << "adding segment exp " << b_exp << "\n";
    auto hh_begin = h_it;
    h_it = insert_ze<true, false>(h_it, b_n);
    if (e_exps_it->second == b_exp)
    {
        auto next = std::next(e_exps_it)->first;
        h_it = std::copy(e_exps_it->first, next, h_it);
        //std::cout << "sorting range: " << std::distance(h_begin, hh_begin) << " to " 
        //          << std::distance(h_begin, h_it) << "\n";
        std::sort(hh_begin, h_it, [](const FPT& a, const FPT& b)
                {
                    return std::abs(a) < std::abs(b);
                });
        ++e_exps_it;
    }
    while (e_exps_it->second != exp_end)
    {
        //std::cout << "adding upper segment exp " << e_exps_it->second << "\n";
        *h_exps_it++ = std::pair{h_it, e_exps_it->second};
        auto next = std::next(e_exps_it)->first;
        auto h_old_it = h_it;
        h_it = std::copy(e_exps_it->first, next, h_it);
        //std::cout << "copied " << std::distance(e_exps_begin->first, e_exps_it->first) 
        //          << " to " << std::distance(e_exps_begin->first, next) << " over " 
        //          << std::distance(h_begin, h_old_it) << " to "
        //          << std::distance(h_begin, h_it) << "\n";
        ++e_exps_it;
    }
    if constexpr (NegateE)
    {
        for (auto hh_it = h_begin; hh_it != h_it; ++hh_it)
        {
            *hh_it = -*hh_it;
        }
    }
    *h_exps_it++ = std::pair {h_it, exp_end};
    return h_exps_it;
}

template
<
    bool NegateE = false,
    bool NegateF = false,
    typename FPT,
    typename ExpInIter1,
    typename ExpInIter2,
    typename ExpOutIter
>
constexpr auto scaled_merge_expansion(ExpInIter1 e_exps_begin,
                                      ExpInIter1,
                                      ExpInIter2 f_exps_begin,
                                      ExpInIter2,
                                      ExpOutIter h_exps_begin,
                                      ExpOutIter)
{
    auto e_exps_it = e_exps_begin;
    auto f_exps_it = f_exps_begin;
    auto h_exps_it = h_exps_begin;
    auto h_begin = h_exps_begin->first;
    auto h_it = h_begin;

    while (e_exps_it->second != exp_end && f_exps_it->second != exp_end)
    {
        if (e_exps_it->second < f_exps_it->second)
        {
            auto h_old_it = h_it;
            *h_exps_it++ = std::pair {h_it, e_exps_it->second};
            auto next = std::next(e_exps_it)->first;
            h_it = std::copy(e_exps_it->first, next, h_it);
            ++e_exps_it;
            if constexpr (NegateE)
            {
                for (auto hh_it = h_old_it; hh_it != h_it; ++hh_it)
                    *hh_it = -*hh_it;
            }
        }
        else if (f_exps_it->second < e_exps_it->second)
        {
            auto h_old_it = h_it;
            *h_exps_it++ = std::pair {h_it, f_exps_it->second};
            auto next = std::next(f_exps_it)->first;
            h_it = std::copy(f_exps_it->first, next, h_it);
            ++f_exps_it;
            if constexpr (NegateF)
            {
                for (auto hh_it = h_old_it; hh_it != h_it; ++hh_it)
                    *hh_it = -*hh_it;
            }
        }
        else // f_exps_it->second == e_exps_it->second
        {
            auto h_old_it = h_it;
            *h_exps_it++ = std::pair {h_it, e_exps_it->second};
            auto next = std::next(e_exps_it)->first;
            h_it = std::copy(e_exps_it->first, next, h_it);
            ++e_exps_it;
            if constexpr (NegateE)
            {
                for (auto hh_it = h_old_it; hh_it != h_it; ++hh_it)
                    *hh_it = -*hh_it;
            }
            auto h_old_it2 = h_it;
            next = std::next(f_exps_it)->first;
            h_it = std::copy(f_exps_it->first, next, h_it);
            ++f_exps_it;
            if constexpr (NegateF)
            {
                for (auto hh_it = h_old_it2; hh_it != h_it; ++hh_it)
                    *hh_it = -*hh_it;
            }
            std::inplace_merge(h_old_it, h_old_it2, h_it, abs_comp{});
        }
    }
    while (e_exps_it->second != exp_end)
    {
        auto h_old_it = h_it;
        *h_exps_it++ = std::pair {h_it, e_exps_it->second};
        auto next = std::next(e_exps_it)->first;
        h_it = std::copy(e_exps_it->first, next, h_it);
        ++e_exps_it;
        if constexpr (NegateE)
        {
            for (auto hh_it = h_old_it; hh_it != h_it; ++hh_it)
                *hh_it = -*hh_it;
        }
    }
    while (f_exps_it->second != exp_end)
    {
        auto h_old_it = h_it;
        *h_exps_it++ = std::pair {h_it, f_exps_it->second};
        auto next = std::next(f_exps_it)->first;
        h_it = std::copy(f_exps_it->first, next, h_it);
        ++f_exps_it;
        if constexpr (NegateF)
        {
            for (auto hh_it = h_old_it; hh_it != h_it; ++hh_it)
                *hh_it = -*hh_it;
        }
    }
    return std::pair {h_it, h_exps_it};
}
/*
template
<
    bool NegateE = false,
    bool NegateF = false,
    typename FPT,
    typename ExpInIter1,
    typename ExpInIter2,
    typename ExpOutIter
>
constexpr ExpOutIter scaled_merge_expansion_inplace(ExpInIter1 e_exps_begin,
                                              ExpInIter1 e_exps_end,
                                              ExpInIter2 f_exps_begin,
                                              ExpInIter2 f_exps_end,
                                              ExpOutIter h_exps_begin,
                                              ExpOutIter h_exps_end)
{
    auto e_exps_it = e_exps_begin;
    auto f_exps_it = f_exps_begin;
    auto h_exps_it = h_exps_begin;
    auto h_begin = h_exps_begin->first;
    auto h_it = h_begin;
    auto e_end = e_exps_begin->first;

    assert(e_exps_begin->first == h_exps_begin->first);
    assert(e_exps_begin == h_exps_begin);
    if constexpr (NegateE)
        for (auto e_exps_it2 = e_exps_begin; e_exps_it2->second != exp_end; ++e_exps_it2)
        {
            e_end = std::next(e_exps_it2)->first;
            for (auto e_it = e_exps_it2->first; e_it != e_end; ++e_it)
                *e_it = -*e_it;
        }
    if constexpr (NegateF)
        for (auto f_exps_it2 = f_exps_begin; f_exps_it2->second != exp_end; ++f_exps_it2)
            for (auto f_it = f_exps_it2->first; f_it != std::next(f_exps_it2)->first; ++f_it)
                *f_it = -*f_it;

    while (e_exps_it->second != exp_end && f_exps_it->second != exp_end)
    {
        if (e_exps_it->second < f_exps_it->second)
        {
            h_it = ++e_exps_it->first;
        }
        else if (f_exps_it->second < e_exps_it->second)
        {
            auto f_exp = *f_exps_it;
            auto next = std::next(f_exps_it)->first;
            auto diff = std::distance(f_exps_it->first, next);
            for (auto ee_exps_it = e_exps_it; ee_exps_it != e_exps_end; ++ee_exps_it)
                ee_exps_it->first += diff;
            std::rotate(e_end, f_exps_it->first, next);
            std::rotate(h_it, e_end, e_end + diff);
            e_end += diff;
            std::rotate(e_exps_it, e_exps_end, e_exps_end + 1);
            *e_exps_it = f_exp;
            ++e_exps_it;
            ++e_exps_end;
            ++f_exps_it;
            h_it += diff;
        }
        else // f_exps_it->second == e_exps_it->second
        {
            auto h_old_it = h_it;
            h_it = ++e_exps_it->first;
            auto f_exp = *f_exps_it;
            auto next = std::next(f_exps_it)->first;
            auto diff = std::distance(f_exps_it->first, next);
            for (auto ee_exps_it = e_exps_it; ee_exps_it != e_exps_end; ++ee_exps_it)
                ee_exps_it->first += diff;
            std::rotate(e_end, f_exps_it->first, next);
            std::rotate(h_it, e_end, e_end + diff);
            e_end += diff;
            std::rotate(e_exps_it, e_exps_end, e_exps_end + 1);
            *e_exps_it = f_exp;
            ++e_exps_it;
            ++e_exps_end;
            ++f_exps_it;
            std::inplace_merge(h_old_it, h_it, h_it + diff, abs_comp{});
            h_it += diff;
        }
    }
    h_exps_it = e_exps_end;
    h_it = e_end;
    while (f_exps_it != f_exps_end)
    {
        auto h_old_it = h_it;
        *h_exps_it++ = std::pair {h_it, f_exps_it->second};
        auto next = std::next(f_exps_it)->first;
        std::rotate(h_it, f_exps_it->first, next);
        auto diff = std::distance(f_exps_it->first, next);
        h_it += diff;
        ++f_exps_it;
    }
    return h_exps_it;
}
*/
template <typename ExpInIter, typename ExpOutIter, typename FPOutIter>
constexpr ExpOutIter relocate_scaled_expansion(ExpInIter e_exps_begin,
                                               ExpInIter e_exps_end,
                                               ExpOutIter h_exps_begin,
                                               FPOutIter h_begin)
{
    auto h_exps_end = std::copy(e_exps_begin, e_exps_end, h_exps_begin);
    const auto e_begin = e_exps_begin->first;
    auto dist = std::distance(e_begin, h_begin);
    auto e_end = e_begin;
    auto h_exps_it = h_exps_begin;
    auto e_exps_it = e_exps_begin;
    for(; e_exps_it->second != exp_end; ++e_exps_it)
        (h_exps_it++)->first += dist;
    h_exps_it->first += dist;
    e_end = e_exps_it->first;
    auto h_end = std::copy(e_begin, e_end, h_begin);
    assert(h_exps_it->first == h_end);
    return ++h_exps_it;
}

template <int degree = 1, typename ExpIter>
constexpr ExpIter renormalize_scaled_expansion(ExpIter h_exps_begin,
                                               ExpIter h_exps_end)
{
    using FPT = std::remove_cvref_t<decltype(*h_exps_begin->first)>;
    FPT Q = 0;
    auto Q_exp = h_exps_begin->second - 2;
    const auto h_begin = h_exps_begin->first;
    auto h_it = h_begin;
    int i = 0;
/*    for (auto it = h_exps_begin; it->second != exp_end; ++it)
        std::cout << "segment " << i 
                  << " exp: " << it->second
                  << " from " << std::distance(h_begin, it->first)
                  << " to " << std::distance(h_begin, std::next(it)->first)
                  << "\n";*/
    auto h_exps_it = h_exps_begin;
    for (; h_exps_it->second != exp_end; ++h_exps_it)
    {
        auto h_exp = h_exps_it->second;
        auto next = std::next(h_exps_it)->first;
        auto g_it = h_exps_it->first;
        /*std::cout << "processing segment " << std::distance(h_exps_begin, h_exps_it)
                  << ", carrying " << Q
                  << " * (SF ^ " << Q_exp << "), "
                  << "h_exp: " << h_exp
                  << "; from: " << std::distance(h_begin, g_it)
                  << " to: " << std::distance(h_begin, next) << "\n";*/
        if (    std::abs(Q) >= scaling_upper<FPT> * std::numeric_limits<FPT>::epsilon()
             && Q_exp == h_exp - 1 )
        {
            Q /= scaling_factor<FPT, degree>;
            const auto Q_new = Q + *g_it;
            auto h_new = fast_two_sum_tail(*g_it, Q, Q_new);
            h_it = insert_ze<true>(h_it, h_new * scaling_factor<FPT, degree>);
            Q = Q_new;
        }
        else
        {
            h_it = insert_ze<true>(h_it, Q);
            Q = *g_it;
        }
        Q_exp = h_exp;
        h_exps_it->first = h_it;
        ++g_it;/*
        std::cout << "Going into loop "
                  << " carrying " << Q << "  * (SF ^ " << Q_exp << ")\n";*/
        for (; g_it != next; ++g_it)
        {
            auto Q_new = Q + *g_it;
            auto h_new = fast_two_sum_tail(*g_it, Q, Q_new);/*
            std::cout << "" << Q_new << "\t = "
                      << "" << Q << " + \t"
                      << "" << *g_it << "; R: \t"
                      << "" << h_new << "\n";*/
            h_it = insert_ze<true>(h_it, h_new);
            Q = Q_new;
        }
 //       debug_scaled_expansion::print(h_exps_begin, h_it);
    }
    h_it = insert_ze<true>(h_it, Q);
    h_exps_it->first = h_it;
    for (auto h_exps_it = h_exps_begin; h_exps_it->second != exp_end; ++h_exps_it)
    {
        auto next = std::next(h_exps_it)->first;
        for (auto hh_it = h_exps_it->first; hh_it != next; ++hh_it)
        {
            if (std::abs(*hh_it) <= scaling_lower<FPT, degree>)
            {
                *hh_it *= scaling_factor<FPT, degree>;
                h_exps_it->first++;
            }
            if (std::abs(*hh_it) > scaling_upper<FPT>)
            {
                *hh_it /= scaling_factor<FPT, degree>;
                if (h_exps_it + 1 != h_exps_end)
                    (h_exps_it + 1)->first--;
                else
                {
                    *h_exps_end++ = std::pair {hh_it, h_exps_it->second + 1};
                }
            }
        }
    }
    if (h_exps_begin->first != h_begin)
    {
        std::rotate(h_exps_begin, h_exps_end, h_exps_end + 1);
        h_exps_end++;
        *h_exps_begin = std::pair {h_begin, (h_exps_begin+1)->second - 1};
    }
    for (auto h_exps_it = h_exps_begin;
         h_exps_it != h_exps_end - 1 && h_exps_it != h_exps_end;
         ++h_exps_it)
    {
        if(h_exps_it->first == (h_exps_it + 1)->first)
        {
            std::rotate(h_exps_it, h_exps_it + 1, h_exps_end);
            h_exps_end--;
        }
    }
    return h_exps_end;
}

template
<
    int degree = 1,
    bool NegateE = false,
    bool NegateB = false,
    typename FPT,
    typename ExpInIter,
    typename ExpOutIter,
    typename FPOutIter
>
constexpr auto scaled_scale_expansion(ExpInIter e_exps_begin,
                                      ExpInIter,
                                      FPT b,
                                      ExpOutIter h_exps_begin,
                                      FPOutIter h_begin)
{
    auto h_it = h_begin;
    auto e_exps_it = e_exps_begin;
    auto h_exps_it = h_exps_begin;
    const auto [b_n, b_exp] = scale_number(b);

    FPT h_new = 0;
    FPT Q = 0;
    auto Q_exp = e_exps_begin->second - 3;

    *h_exps_it++ = std::pair {h_it, e_exps_it->second + b_exp};
    for (; e_exps_it->second != exp_end; ++e_exps_it)
    {
        auto e_it = e_exps_it->first;
        const auto T = *e_it * b_n;
        const auto t = two_product_tail(*e_it, b_n, T);
        auto Q_new_exp = e_exps_it->second + b_exp;
        if (   Q_exp == Q_new_exp - 1
            && std::abs(Q) >= scaling_upper<FPT> * std::numeric_limits<FPT>::epsilon())
        {
            Q /= scaling_factor<FPT, degree>;
            auto Q_new = Q + t;
            auto h_new = two_sum_tail(Q, t, Q_new);
            h_it = insert_ze<true>(h_it, h_new * scaling_factor<FPT, degree>);
            Q = Q_new;
        }
        else
        {
            h_it = insert_ze<true>(h_it, Q);
            Q = t;
        }
        Q_exp = Q_new_exp;
        for (;;);
    }
}

template <typename Expression>
constexpr int degree()
{
    if constexpr (Expression::operator_type == operator_types::product)
    {
        if constexpr (std::is_same_v<typename Expression::left, typename Expression::right>)
            return 2 * degree<typename Expression::left>();
        else
            return   degree<typename Expression::left>()
                   + degree<typename Expression::right>();
    }
    else if constexpr (   Expression::operator_type == operator_types::sum
                       || Expression::operator_type == operator_types::difference)
    {
        return std::max(degree<typename Expression::left>(),
                        degree<typename Expression::right>());
    }
    else if constexpr (Expression::is_leaf)
    {
        return 1;
    }
    return -999;
}

template <typename Expression>
constexpr int summands()
{
    if constexpr (Expression::operator_type == operator_types::product)
        return   summands<typename Expression::left>()
               * summands<typename Expression::right>();
    else if constexpr (   Expression::operator_type == operator_types::sum
                       || Expression::operator_type == operator_types::difference)
        return   summands<typename Expression::left>()
               + summands<typename Expression::right>();
    else if constexpr (Expression::is_leaf)
        return 1;
    return -999;
}

constexpr std::size_t nonzeroes(auto decomp)
{
    std::size_t i = 0;
    for (const auto& d : decomp)
        if(d.second != 0)
            ++i;
    return i;
}

constexpr std::size_t nofac = std::numeric_limits<std::size_t>::max();

constexpr void aggregate(auto& decomp)
{
    int i = 0;
    for(int j = 1; j < decomp.size(); ++j)
    {
        if( decomp[i].first == decomp[j].first )
        {
            decomp[i].second += decomp[j].second;
            decomp[j].second = 0;
        }
        else
        {
            if(decomp[i].second != 0)
                ++i;
            decomp[i] = decomp[j];
        }
    }
    for(++i; i < decomp.size(); ++i)
        decomp[i].second = 0;
}

template <typename Expression, typename Real>
constexpr auto decompose()
{
    constexpr auto deg = degree<Expression>();
    if constexpr (Expression::is_leaf)
    {
        std::array<std::pair<std::array<std::pair<std::size_t, Real>, deg>, int>, 1> out;
        if constexpr (Expression::argn == 0)
        {
            out[0].first[0] = std::pair { std::size_t{0}, Real(Expression::value) };
            static_assert( Real(Expression::value) == Expression::value );
        }
        else
            out[0].first[0] = std::pair { Expression::argn, Real{} };
        out[0].second = 1;
        return out;
    }
    else
    {
        int i = 0;
        if constexpr (Expression::operator_type == operator_types::product)
        {
            if constexpr (std::is_same_v<typename Expression::left, typename Expression::right>)
            {
                constexpr auto left = decompose<typename Expression::left, Real>();
                constexpr std::size_t leftsize = nonzeroes(left);
                std::array
                    <
                        std::pair<std::array<std::pair<std::size_t, Real>, deg>, int>,
                        leftsize * (leftsize + 1) / 2
                    > out;
                for (int li = 0; li < leftsize; ++li)
                {
                    for (int j = 0; j < deg / 2; ++j)
                        out[i].first[ 2 * j ] = out[i].first[ 2 * j + 1 ] = left[li].first[j];
                    out[i].second = left[li].second * left[li].second;
                    ++i;
                }
                for(int li = 0; li < leftsize; ++li)
                {
                    for (int li2 = li + 1; li2 < leftsize; ++li2)
                    {
                        std::merge(left[li].first.begin(), left[li].first.end(),
                                   left[li2].first.begin(), left[li2].first.end(),
                                   out[i].first.begin());
                        out[i].second = 2 * left[li].second * left[li2].second; 
                        ++i;
                    }
                }
                std::sort(out.begin(), out.end());
                aggregate(out);
                return out;
            }
            else
            {
                constexpr auto left = decompose<typename Expression::left, Real>();
                constexpr std::size_t leftsize = nonzeroes(left);
                constexpr auto right = decompose<typename Expression::right, Real>();
                constexpr std::size_t rightsize = nonzeroes(right);
                std::array
                    <
                        std::pair<std::array<std::pair<std::size_t, Real>, deg>, int>,
                        leftsize * rightsize
                    > out;
                for(int li = 0; li < leftsize; ++li)
                    for (auto ri = 0; ri < rightsize; ++ri)
                    {
                        std::merge(left[li].first.begin(), left[li].first.end(),
                                   right[ri].first.begin(), right[ri].first.end(),
                                   out[i].first.begin());
                        out[i].second = left[li].second * right[ri].second;
                        ++i;
                    }
                std::sort(out.begin(), out.end());
                aggregate(out);
                return out;
            }
        }
        else if constexpr (   Expression::operator_type == operator_types::sum
                           || Expression::operator_type == operator_types::difference)
        {
            constexpr auto left = decompose<typename Expression::left, Real>();
            constexpr std::size_t leftsize = nonzeroes(left);
            constexpr auto right = decompose<typename Expression::right, Real>();
            constexpr std::size_t rightsize = nonzeroes(right);
            std::array
                <
                    std::pair<std::array<std::pair<std::size_t, Real>, deg>, int>,
                    leftsize + rightsize
                > out;
            constexpr bool negate_r = Expression::operator_type == operator_types::difference;
            auto out_temp = out;
            for (int li = 0; li < leftsize; ++li)
            {
                for (int j = 0; j < deg; ++j)
                    if (j < left[li].first.size())
                        out_temp[i].first[j] = left[li].first[j];
                    else
                        out_temp[i].first[j] = std::pair {nofac, Real{}};
                out_temp[i].second = left[li].second;
                ++i;
            }
            for (int ri = 0; ri < rightsize; ++ri)
            {
                for (int j = 0; j < deg; ++j)
                    if (j < right[ri].first.size())
                        out_temp[i].first[j] = right[ri].first[j];
                    else
                        out_temp[i].first[j] = std::pair {nofac, Real{}};
                out_temp[i].second = negate_r ? -right[ri].second : right[ri].second;
                ++i;
            }
            std::merge(out_temp.begin(), out_temp.begin() + leftsize,
                       out_temp.begin() + leftsize, out_temp.begin() + leftsize + rightsize,
                       out.begin()); //inplace_merge not constexpr
            aggregate(out);
            return out;
        }
    }
}

template <typename Expression, typename Real = double>
int scaled_sign(const auto& input)
{
    constexpr auto summands = decompose<Expression, Real>();
    constexpr auto deg = degree<Expression>();
    constexpr auto s = nonzeroes(summands);

    constexpr auto scaling_factor_exp = max_exp<Real> - min_exp<Real, deg>;

    constexpr auto exp_upper_bound =
        (std::numeric_limits<Real>::max_exponent * deg - max_exp<Real>) / scaling_factor_exp;
    constexpr auto exp_lower_bound =
        ((std::numeric_limits<Real>::min_exponent + 1 - std::numeric_limits<Real>::digits) * deg - min_exp<Real>) / scaling_factor_exp - 1;
    constexpr auto levels = exp_upper_bound - exp_lower_bound + 1;
    constexpr auto exp_index_offset = -exp_lower_bound;

    std::array<Real, ce_exp2<std::size_t>(deg)> tmp;
    std::fill(tmp.begin(), tmp.end(), Real(0));
    std::array<std::array<Real, deg + 1>, s> products;
    for (auto& p : products)
        std::fill(p.begin(), p.end(), Real(0));
    std::array<Real, s * (deg + 1)> all_comps;
    std::fill(all_comps.begin(), all_comps.end(), Real(0));
    std::array<int, s> scales;
    std::array<typename std::array<Real, deg + 1>::iterator, s> lower_scale_ends, products_end;
    std::fill(scales.begin(), scales.end(), 0);

    std::array<int, levels> level_comps;
    std::array<typename std::array<Real, s * (deg + 1)>::iterator, levels> level_ends, level_begins;
    std::fill(level_comps.begin(), level_comps.end(), 0);

    for (int i = 0; i < s; ++i)
    {
        bool zero = false;
        for (const auto& fac : summands[i].first)
            if (fac.first != 0 && fac.first != nofac && input[fac.first] == 0)
            {
                zero = true;
                break;
            }
        if (zero)
            continue;
        products[i][0] = Real(summands[i].second);
        products_end[i] = products[i].begin() + 1;
        for (int j = 0; j < summands[i].first.size() && summands[i].first[j].first != nofac; ++j)
        {
            const Real b = input[summands[i].first[j].first];
            const auto [b_n, b_exp] = scale_number<deg>(b);
            scales[i] += b_exp;
            auto tmp_end = scale_expansion<true>(products[i].begin(),
                                           products_end[i],
                                           b_n,
                                           tmp.begin(),
                                           tmp.end());
            products_end[i] = vec_sum_err_branch(tmp.begin(), tmp_end, products[i].begin());
            if (std::abs(products[i].front()) > scaling_upper<Real>)
            {
                ++scales[i];
                for (auto prod_it = products[i].begin(); prod_it != products_end[i]; ++prod_it)
                    *prod_it /= scaling_factor<Real, deg>;
            }
            else if (   std::abs(products[i].front()) != 0
                     && std::abs(products[i].front()) <= scaling_lower<Real, deg>)
            {
                --scales[i];
                for (auto prod_it = products[i].begin(); prod_it != products_end[i]; ++prod_it)
                    *prod_it *= scaling_factor<Real, deg>;
            }
            std::reverse(products[i].begin(), products_end[i]);
        }
        auto prod_length = std::distance(products[i].begin(), products_end[i]);
        lower_scale_ends[i] = products[i].begin();
        while( *lower_scale_ends[i] != 0 && std::abs(*lower_scale_ends[i]) <= scaling_lower<Real, deg>)
            lower_scale_ends[i]++;
        const auto lower_length = std::distance(products[i].begin(), lower_scale_ends[i]);
        level_comps[exp_index_offset + scales[i]] += (prod_length - lower_length);
        level_comps[exp_index_offset + scales[i] - 1] += lower_length;
    }
    level_begins[0] = level_ends[0] = all_comps.begin();
    for(int i = 0; i < levels - 1; ++i)
    {
        level_begins[i + 1] = level_ends[i + 1] = level_begins[i] + level_comps[i];
    }
    for(int i = 0; i < s; ++i)
    {
         level_ends[scales[i] + exp_index_offset - 1] = 
             std::copy(products[i].begin(), lower_scale_ends[i],
                       level_ends[scales[i] + exp_index_offset - 1]);
         level_ends[scales[i] + exp_index_offset] = 
             std::copy(lower_scale_ends[i], products_end[i],
                       level_ends[scales[i] + exp_index_offset]);

    }
    for(int i = 0; i < levels; ++i)
    {
        std::sort(level_begins[i], level_ends[i],
                  [](Real a, Real b) {return std::abs(a) > std::abs(b);});
    }
    Real remainder_fac = (Real(1.) - std::numeric_limits<Real>::epsilon() * 4) / s;
    Real Q = 0, q = 0;
    int level = levels - 1;
    for (int i = level; i >= 0; --i)
    {
        Q *= scaling_factor<Real, deg>;
        q *= scaling_factor<Real, deg>;
        for (auto c_it = level_begins[i]; c_it != level_ends[i]; ++c_it )
        {
            if( std::abs(*c_it) < std::abs(Q) * remainder_fac )
                return (Q > 0) - (Q < 0);
            auto Q_new = Q + *c_it;
            auto tail = two_sum_tail(Q, *c_it, Q_new);
            auto q_new = q + tail;
            assert(two_sum_tail(q, tail, q_new) == Real(0));
            auto Q_new2 = Q_new + q_new;
            auto q_new2 = two_sum_tail(Q_new, q_new, Q_new2);
            Q = Q_new2;
            q = q_new2;
        }
    }
    return (Q > 0) - (Q < 0);
}

constexpr std::size_t count_unique_summands(const auto decomp)
{
    if(decomp.size() == 0) return 0;
    int i = 0;
    int j = 0;
    auto cur = decomp[0].first;
    for(const auto summand : decomp)
    {
        if (cur == summand.first)
        {
            if (summand.second)
                --j;
            else
                ++j;
        }
        else
        {
            if (j != 0)
                ++i;
            cur = summand.first;
            if (summand.second)
                j = -1;
            else
                j = 1;
        }
    }
    if (j != 0) ++i;
    return i;
}

template <typename Expression, typename Real>
constexpr auto decompose_reduced()
{
    constexpr auto deg = degree<Expression>();
    constexpr auto decomp = decompose<Expression, Real>();
    constexpr auto unique = count_unique_summands(decomp);
    std::array
        <
            std::pair<std::array<std::pair<std::size_t, Real>, deg>, int>,
            unique
        > out;

    if(decomp.size() == 0) return out;
    int i = 0;
    int j = 0;
    auto cur = decomp[0].first;
    for(const auto summand : decomp)
    {
        if (cur == summand.first)
        {
            if (summand.second)
                --j;
            else
                ++j;
        }
        else
        {
            if (j != 0)
                out[i++] = std::pair{cur, j};
            cur = summand.first;
            if (summand.second)
                j = -1;
            else
                j = 1;
        }
    }
    if (j != 0)
        out[i++] = std::pair{cur, j};
    return out;
}

}} // detail::generic_robust_predicates

}} // boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_SCALED_EXPANSION_ARITHMETIC_HPP

