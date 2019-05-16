//
// Created by Jan Groschaft on 4.4.19.
//

/*
 * Module for splitting network into multiple segments. Current implementation divides the network into multiple segments
 * of similar cost, where the cost of a segment increases with it's height and the number of active vertices in it.
 */

#ifndef MAXFLOW_PARTITIONING_H
#define MAXFLOW_PARTITIONING_H

#include <memory>
#include <vector>
#include <cmath>

namespace partitioning
{
    namespace detail
    {
        template <typename T>
        auto get_cost ( T total_active, T height ) noexcept
        {
            auto log = log2 ( height + 1 );
            return height * log * log * total_active;
        }

        template <typename label_info, typename T>
        auto get_next_partition ( const std::unique_ptr<label_info[]> & labels, const double cost_per_segment, const T low, const T highest_active ) noexcept
        {
            auto total_active = T { 0 };
            for ( T i = low + 1; i <= highest_active; ++i )
            {
                total_active += labels[i] . active_vertices . size ();
                double cost = get_cost ( total_active, i - low );
                if ( cost >= cost_per_segment )
                    return std::make_tuple ( i + 1, cost );
            }
            return std::make_tuple ( highest_active + 1, get_cost ( total_active, highest_active - low ) );
        }


        template <typename label_info, typename T>
        auto calculate_partitions ( const std::unique_ptr<label_info[]> & labels, const double cost_per_segment, const T highest_active,
                                    std::vector<T> & partitions, const std::size_t thread_count ) noexcept
        {
            static constexpr T min_segment_height = T { 1 } << 2;
            T low = 0;
            double total_cost = 0;

            partitions . clear ();
            partitions . emplace_back ( low );

            while ( partitions . size () < thread_count )
            {
                auto[high, cost] = get_next_partition ( labels, cost_per_segment, low, highest_active );
                if ( high >= highest_active - 1 )
                    break;
                total_cost += cost;
                high = std::max ( high, low + min_segment_height );
                partitions . emplace_back ( high );
                low = high;
            }

            partitions . emplace_back ( highest_active + 1 );
            double cost; //we have to use std::ignore, because otherwise we will get unused variable warning for the first member...
            std::tie ( std::ignore, cost ) = get_next_partition ( labels, std::numeric_limits<double>::max (), low, highest_active );
            total_cost += cost;
            return total_cost;
        }
    }

    template <typename label_info, typename T>
    std::vector<T> get_partitions ( const std::unique_ptr<label_info[]> & labels, const T highest_active,
                                    std::size_t thread_count, const T min_active_per_thread ) noexcept
    {
        T total_active_cnt = 0;
        for ( std::size_t i = 0; i <= highest_active; ++i )
            total_active_cnt += labels[i] . active_vertices . size ();

        //ensure there are enough active vertices, adjust thread count if not
        if ( total_active_cnt < thread_count * min_active_per_thread )
            thread_count = std::max ( T { 1 }, total_active_cnt / min_active_per_thread );

        std::vector<T> partitions, result = { 0, highest_active + 1 };
        partitions . reserve ( thread_count );

        if ( thread_count == 1 )
            return result;

        const double max_cost = detail::get_cost ( total_active_cnt, highest_active + 1 );
        const double min_cost = 1;
        double lowest_cost = max_cost;
        double low = min_cost, high = max_cost;
        double mid;

        while ( low <= high )
        {
            mid = ( low + high ) / 2;

            auto total_cost = detail::calculate_partitions ( labels, mid, highest_active, partitions, thread_count );
            if ( total_cost <= lowest_cost )
            {
                lowest_cost = total_cost;
                std::swap ( result, partitions );
                high = std::ceil ( mid - 1 );
            } else
                low = std::floor ( mid + 1 );
        }

        return result;
    }
}


#endif //MAXFLOW_PARTITIONING_H
