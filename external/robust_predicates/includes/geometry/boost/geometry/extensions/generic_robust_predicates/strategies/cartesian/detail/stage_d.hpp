// Boost.Geometry (aka GGL, Generic Geometry Library)

// Copyright (c) 2020 Tinko Bartels, Berlin, Germany.

// Contributed and/or modified by Tinko Bartels,
//   as part of Google Summer of Code 2020 program.

// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_D_HPP
#define BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_D_HPP

#include <cstddef>
#include <array>
#include <algorithm>

#include <boost/mp11/integral.hpp>
#include <boost/mp11/list.hpp>
#include <boost/mp11/function.hpp>
#include <boost/mp11/algorithm.hpp>

#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expression_tree.hpp>
#include <boost/geometry/extensions/generic_robust_predicates/strategies/cartesian/detail/expansion_eval.hpp>

namespace boost { namespace geometry
{

namespace detail { namespace generic_robust_predicates
{
template <
    auto Expression,
    typename Real = double,
    template <int> class ZEPolicy = default_zero_elimination_policy,
    template <int, int> class FEPolicy = default_fast_expansion_sum_policy
>
// FERNANDO: Modified this!
struct stage_d
{
    static constexpr bool stateful = false;
    static constexpr bool updates = false;

    using expression_t = decltype(Expression);
    using stack             = typename boost::mp11::mp_unique<post_order<expression_t>>;
    using evals             = typename boost::mp11::mp_remove_if<stack, is_leaf>;
    using sizes_pre         = boost::mp11::mp_transform <
                                  expansion_size,
                                  boost::mp11::mp_pop_back<evals>>;
    using sizes             = boost::mp11::mp_push_back <
                                  sizes_pre,
                                  boost::mp11::mp_size_t <
                                      final_expansion_size <
                                          expression_t,
                                          expansion_size<typename expression_t::left>::value,
                                          expansion_size<typename expression_t::right>::value
                                      >()>>;
    using accumulated_sizes = boost::mp11::mp_push_front <
                                  boost::mp11::mp_partial_sum <
                                      sizes,
                                      boost::mp11::mp_size_t<0>,
                                      boost::mp11::mp_plus>,
                                  boost::mp11::mp_size_t<0>>;

    static constexpr std::size_t results_size =
        boost::mp11::mp_back<accumulated_sizes>::value;

    // Threshold above which we use heap storage. Sized conservatively to
    // stay well within typical 8 MB stack limits and avoid cache pressure.
    // 4096 doubles = 32 KB, comfortably on stack; anything larger goes heap.
    static constexpr std::size_t stack_threshold = 4096;

    template <typename ...Reals>
    static inline int apply(const Reals&... args)
    {
        auto run = [&](Real* rbegin, Real* rend) -> int {
            auto most_significant = eval_expansions <
                evals,
                sizes,
                accumulated_sizes,
                Real*,
                Real,
                false,
                ZEPolicy,
                FEPolicy
            >(rbegin, rend, args...) - 1;

            if (*most_significant == 0) return  0;
            if (*most_significant  > 0) return  1;
            return -1;
        };

        if constexpr (results_size > stack_threshold) {
            // Thread-local buffer: allocated once per thread, reused forever.
            // If this stage_d is instantiated with multiple expression sizes,
            // each instantiation has its own thread_local buffer because
            // thread_local variables at function scope are per-instantiation.
            static thread_local std::vector<Real> buf;
            if (buf.size() < results_size)
                buf.resize(results_size);
            return run(buf.data(), buf.data() + results_size);
        } else {
            std::array<Real, results_size> buf;
            return run(buf.data(), buf.data() + results_size);
        }
    }
};


// template <
//     auto Expression,
//     typename Real = double,
//     template <int> class ZEPolicy = default_zero_elimination_policy,
//     template <int, int> class FEPolicy = default_fast_expansion_sum_policy
// >
// struct stage_d
// {
//     static constexpr bool stateful = false;
//     static constexpr bool updates = false;
//
//     using expression_t = decltype(Expression);
//
//     // Move these from inside apply() to here:
//     using stack = typename boost::mp11::mp_unique<post_order<expression_t>>;
//     using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
//     using sizes_pre = boost::mp11::mp_transform<
//         expansion_size,
//         boost::mp11::mp_pop_back<evals>
//     >;
//     using sizes = boost::mp11::mp_push_back<
//         sizes_pre,
//         boost::mp11::mp_size_t<
//             final_expansion_size<
//                 expression_t,
//                 expansion_size<typename expression_t::left>::value,
//                 expansion_size<typename expression_t::right>::value
//             >()
//         >
//     >;
//     using accumulated_sizes = boost::mp11::mp_push_front<
//         boost::mp11::mp_partial_sum<
//             sizes,
//             boost::mp11::mp_size_t<0>,
//             boost::mp11::mp_plus
//         >,
//         boost::mp11::mp_size_t<0>
//     >;
//
//     // The total number of doubles in the results buffer:
//     static constexpr std::size_t results_size =
//         boost::mp11::mp_back<accumulated_sizes>::value;
//
//     template <typename ...Reals>
//     static inline int apply(const Reals&... args)
//     {
//         // Remove the duplicate local using-declarations that were here before.
//         // results_size is now available directly:
//         // std::vector<Real> results(results_size);  // heap allocation (your previous fix)
//     
//         std::vector<Real> results(results_size);
//         Real* rbegin = results.data();          // plain Real*, not __wrap_iter
//         Real* rend   = rbegin + results_size;
//
//         auto most_significant = eval_expansions <
//                 evals,
//                 sizes,
//                 accumulated_sizes,
//                 Real*,                          // explicit Iter type
//                 Real,
//                 false,
//                 ZEPolicy,
//                 FEPolicy
//             >(rbegin, rend, args...) - 1;
//
//         if (*most_significant == 0) return  0;
//         if (*most_significant  > 0) return  1;
//         return -1;
//
//         // auto most_significant = eval_expansions<
//         //         evals,
//         //         sizes,
//         //         accumulated_sizes,
//         //         decltype(results.begin()),
//         //         Real,
//         //         false,
//         //         ZEPolicy,
//         //         FEPolicy
//         //     >(results.begin(), results.end(), args...) - 1;
//         //
//         // if (*most_significant == 0) return 0;
//         // else if (*most_significant > 0) return 1;
//         // else return -1;
//     }
// };


// template
// <
//     auto Expression,
//     typename Real = double,
//     template <int> class ZEPolicy = default_zero_elimination_policy,
//     template <int, int> class FEPolicy = default_fast_expansion_sum_policy
// >
// struct stage_d
// {
//     static constexpr bool stateful = false;
//     static constexpr bool updates = false;
//
//     using expression_t = decltype(Expression);
//
//     template <typename ...Reals>
//     static inline int apply(const Reals&... args)
//     {
//         using stack = typename boost::mp11::mp_unique<post_order<expression_t>>;
//         using evals = typename boost::mp11::mp_remove_if<stack, is_leaf>;
//         using sizes_pre = boost::mp11::mp_transform
//             <
//                 expansion_size,
//                 boost::mp11::mp_pop_back<evals>
//             >;
//         using sizes = boost::mp11::mp_push_back
//             <
//                 sizes_pre,
//                 boost::mp11::mp_size_t
//                     <
//                         final_expansion_size
//                             <
//                                 expression_t,
//                                 expansion_size<typename expression_t::left>::value,
//                                 expansion_size<typename expression_t::right>::value
//                             >()
//                     >
//             >;
//         using accumulated_sizes = boost::mp11::mp_push_front
//             <
//                 boost::mp11::mp_partial_sum
//                     <
//                         sizes,
//                         boost::mp11::mp_size_t<0>,
//                         boost::mp11::mp_plus
//                     >,
//                 boost::mp11::mp_size_t<0>
//             >;
//
//         using result_array =
//             std::array<Real, boost::mp11::mp_back<accumulated_sizes>::value>;
//         result_array results;
//
//         auto most_significant = eval_expansions
//             <
//                 evals,
//                 sizes,
//                 accumulated_sizes,
//                 decltype(results.begin()),
//                 Real,
//                 false,
//                 ZEPolicy,
//                 FEPolicy
//             >(results.begin(), results.end(), args...) - 1;
//         if( *most_significant == 0)
//         {
//             return 0;
//         }
//         else if( *most_significant > 0 )
//         {
//             return 1;
//         }
//         else
//         {
//             return -1;
//         }
//     }
// };

template
<
    auto E,
    typename CT,
    template <int> class ZEP,
    template <int, int> class FEP
>
std::string stage_label<stage_d<E, CT, ZEP, FEP >> =
          std::string("stage_d < Expression = ")
        + expression_label<E>
        + std::string(", CalculationType = ")
        + boost::core::demangle( typeid( CT ).name() )
        + std::string(", Zero elimination policy = ")
        + boost::core::demangle( typeid( ZEP<2> ).name() )
        + std::string(", Fast expansion sum policy = ")
        + boost::core::demangle( typeid( FEP<2, 2> ).name() )
        + std::string(" >");

}} // namespace detail::generic_robust_predicates

}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_EXTENSIONS_GENERIC_ROBUST_PREDICATES_STRATEGIES_CARTESIAN_DETAIL_STAGE_D_HPP
