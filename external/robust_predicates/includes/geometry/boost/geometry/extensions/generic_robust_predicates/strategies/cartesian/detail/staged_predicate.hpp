// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGED_PREDICATE_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGED_PREDICATE_HPP

#ifdef CGAL_PROFILE
#include <iostream>
#include <boost/core/demangle.hpp>
#include <vector>
#include <random>
#endif

#include <type_traits>
#include <array>
#include <tuple>

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>
#include <boost/mp11/set.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{

template <typename Stage>
using is_stateful = boost::mp11::mp_bool<Stage::stateful>;

template <typename Stage>
using is_updatable = boost::mp11::mp_bool<Stage::updates>;

#ifdef CGAL_PROFILE

#ifndef CGAL_PROFILE_RESERVOIR
#define CGAL_PROFILE_RESERVOIR 10
#endif


template <typename Stages, typename FailArr>
struct staged_predicate_profiler
{

    static staged_predicate_profiler& instance()
    {
        static staged_predicate_profiler x;
        return x;
    }

    std::array<long, std::tuple_size_v<Stages>> successes = {};
    std::array<long, std::tuple_size_v<Stages>> failures = {};
    std::array<std::vector<FailArr>, std::tuple_size_v<Stages>> failure_records;
    std::mt19937 gen = std::mt19937(std::random_device{}());

    void record_fail(FailArr a, int s)
    {
        auto& f = failure_records[s];
        constexpr int k = CGAL_PROFILE_RESERVOIR;
        auto i = failures[s]++;
        if(i < k)
            f.push_back(a);
        else
        {
            std::uniform_int_distribution<> distrib(0, i - 1);
            const int j = distrib(gen);
            if(j < k)
            {
                f[j] = a;
            }
        }
    }

    ~staged_predicate_profiler()
    {
        long sum = 0;
        std::cerr << "Staged predicate chain: " /*<< boost::core::demangle( typeid(Stages).name() )*/ << "\n";
        std::size_t const N = std::tuple_size_v<Stages>;
        mp11::mp_for_each<mp11::mp_iota_c<N>>( [&]( auto I ){
            const int i = I.value;
            using stage = mp11::mp_at<Stages, decltype(I)>;
            std::cerr << i << " " << successes[i] << " " << stage_label<stage> << "\n";
            sum += successes[i];
#ifdef CGAL_PROFILE_PRINT_FAILURES
            for (const auto& failure : failure_records[i])
            {
                for (const auto& f : failure)
                    std::cerr << f << ", ";
                std::cerr << "\n";
            }
#endif
        });
        std::cerr << "total: " << sum << "\n";
    }
};
#endif

template <std::size_t I, typename Stages, typename ...Reals>
inline int staged_apply(const Stages& stages, const Reals&... args)
{
#ifdef CGAL_PROFILE
    using real = boost::mp11::mp_front<boost::mp11::mp_list<Reals...>>;
    using input_arr = std::array<real, sizeof...(Reals)>;
    using profiler = staged_predicate_profiler<Stages, input_arr>;
#endif
    if constexpr(I == std::tuple_size_v<Stages>)
    {
        return sign_uncertain;
    }
    else
    {
        int result = std::get<I>(stages).apply(args...);
        if( result != sign_uncertain )
        {
#ifdef CGAL_PROFILE
            profiler::instance().successes[I]++;
#endif
            return result;
        }
        else
        {
#ifdef CGAL_PROFILE
            input_arr arr {{static_cast<real>(args)...}};
            profiler::instance().record_fail(arr, I);
#endif
            return staged_apply<I + 1>(stages, args...);
        }
    }
}

template
<
    typename ...Stages
>
struct staged_predicate
{
private:
    using stages = std::tuple<Stages...>;
    stages m_stages;
public:
    static constexpr bool stateful =
        boost::mp11::mp_any_of<stages, is_stateful>::value;
    static constexpr bool updates =
        boost::mp11::mp_any_of<stages, is_updatable>::value;

    staged_predicate() = default;

    staged_predicate(Stages... s) : m_stages(s...) {}

    template <typename ...Reals>
    inline void update(const Reals&... args)
    {
        std::size_t const N = sizeof...(Stages);
        mp11::mp_for_each<mp11::mp_iota_c<N>>( [&]( auto I ){
            if constexpr(mp11::mp_at<stages, decltype(I)>::updates)
            {
                std::get<decltype(I)::value>(m_stages).update(args...);
            }
        });
    }

    template <typename ...Reals>
    inline int apply(const Reals&... args) const
    {
        //std::size_t const N = sizeof...(Stages);
        return staged_apply<0>(m_stages, args...);
    }
};

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGED_PREDICATE_HPP
