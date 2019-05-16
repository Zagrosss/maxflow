//
// Created by Jan Groschaft on 31.3.19.
//

/*
 * Parallel implementation of Ahuja-Orlin's algorithm, divides the network into multiple segments.
 */

#ifndef MAXFLOW_AHUJA_ORLIN_SEGMENT_H
#define MAXFLOW_AHUJA_ORLIN_SEGMENT_H


#include "../../common_types.h"
#include "../../data_structures/linked_list.h"
#include "../../data_structures/thread_local_buffer_pool.h"
#include "partitioning.h"
#include <memory>
#include <cassert>
#include <chrono>
#include <cstring>
#include <omp.h>
#include <cmath>
#include <algorithm>
#include <atomic>

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

namespace ahuja_orlin_segment
{
    template <template <class> typename vector, typename T, typename U>
    class max_flow_instance
    {
        struct alignas (CACHE_LINE_SIZE) vertex
        {
            vertex * next = nullptr;
            vertex * prev = nullptr;
            U excess { 0 };
            T label;
            T original_label;
            std::atomic_flag discovered = ATOMIC_FLAG_INIT;
        };

        struct label_info
        {
            data_structures::linked_list<vertex> active_vertices { };
            data_structures::linked_list<vertex> inactive_vertices { };

            void reset ( ) noexcept
            {
                active_vertices . clear ();
                inactive_vertices . clear ();
            }
        };

        vector<vector<cached_edge<T, U>>> _residual_network;
        std::unique_ptr<label_info[]> _labels;
        std::unique_ptr<vertex[]> _vertices;
        data_structures::thread_local_buffer_pool<T> _pool;
        std::unique_ptr<T[]> _q;
        std::unique_ptr<label_info[]> _thread_local_labels;
        T _source, _sink, _highest_active { 0 }, _highest_vertex { 0 };
        U _max_cap;
        std::size_t _thread_count, _original_relabel_threshold { 0 };
        const std::size_t _max_thread_count;
        int64_t _min_cpu_time_per_phase { 0 };
        std::atomic<std::size_t> _relabel_threshold { 0 };


    public:
        max_flow_instance ( vector<vector<cached_edge<T, U>>> graph, T source, T sink,
                            std::size_t thread_count = static_cast<size_t>(omp_get_max_threads ()) )
                :
                _residual_network ( std::move ( graph ) ),
                _labels ( std::make_unique<label_info[]> ( _residual_network . size () + 1 ) ),
                _vertices ( std::make_unique<vertex[]> ( _residual_network . size () ) ),
                _pool ( data_structures::thread_local_buffer_pool<T> { thread_count, _residual_network . size () } ),
                _q ( std::make_unique<T[]> ( _residual_network . size () ) ),
                _thread_local_labels ( std::make_unique<label_info[]> ( thread_count ) ),
                _source ( source ), _sink ( sink ), _thread_count ( thread_count ),
                _max_thread_count ( thread_count )
        {
            omp_set_num_threads ( static_cast<int> ( _max_thread_count ) );
            init ();
        }

        U find_max_flow ( )
        {
            global_relabel ( _max_cap );
            auto K = static_cast<U> ( std::ceil ( std::log2 ( _max_cap ) ));
            for ( U k = 0; k <= K; ++k )
            {
                auto delta = static_cast<U> ( std::pow ( 2, K - k ));
                set_active_inactive ( delta );
                while ( _highest_active != 0 )
                {
                    parallel_phase ( delta );
                    if ( !is_any_active ( delta ) )
                    {
                        set_original_labels ();
                        break;
                    }
                    global_relabel ( delta );
                }
            }
            return _vertices[_sink] . excess;
        }

        void preflow_to_flow ( )
        {
            std::swap ( _source, _sink );
            _highest_vertex = _residual_network . size ();
            find_max_flow ();
            std::swap ( _source, _sink );
            #ifdef DEBUG
            for ( std::size_t i = 0; i < _residual_network . size(); ++i )
                if ( i != _source && i != _sink )
                    if ( _vertices[i] . excess > 0 )
                        std::cerr << "Excess violation: vertex " << i << ", excess " << _vertices[i] . excess << '\n';
            #endif
        }

        auto steal_network ( )
        {
            return std::move ( _residual_network );
        }

    private:
        static constexpr T ALPHA = 6, BETA = 12;
        static constexpr double GLOBAL_RELABEL_FREQ = 1;
        static constexpr T min_active_per_thread = 10;

        void init ( ) noexcept
        {
            _max_cap = 0;
            #pragma omp parallel for schedule(static) reduction(max:_max_cap)
            for ( std::size_t i = 0; i < _residual_network[_source] . size (); ++i )
            {
                auto & edge = _residual_network[_source][i];
                _max_cap = std::max ( _max_cap, edge . r_capacity );
                _vertices[edge . dst_vertex] . excess = edge . r_capacity;
                edge . reverse_r_capacity += edge . r_capacity;
                _residual_network[edge . dst_vertex][edge . reverse_edge_index] . r_capacity += edge . r_capacity;
                _residual_network[edge . dst_vertex][edge . reverse_edge_index] . reverse_r_capacity -= edge . r_capacity;
                edge . r_capacity = 0;
            }

            std::size_t m = 0;
            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
                m += _residual_network[i] . size ();
            _original_relabel_threshold = ( _residual_network . size () * ALPHA + m / 2 ) * 1;
        }


        void set_active_inactive ( const U delta ) noexcept
        {
            for ( std::size_t i = 0; i <= _highest_vertex; ++i )
                _labels[i] . reset ();

            omp_set_num_threads ( static_cast<int> ( _max_thread_count ) );

            _highest_active = _highest_vertex = 0;
            const auto delta_half = delta / 2;
            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
            {
                if ( i == _source || i == _sink || _vertices[i] . label == _residual_network . size () )
                    continue;

                if ( _vertices[i] . excess > delta_half )
                {
                    _highest_active = std::max ( _highest_active, _vertices[i] . label );
                    _labels[_vertices[i] . label] . active_vertices . push ( &_vertices[i] );
                } else
                    _labels[_vertices[i] . label] . inactive_vertices . push ( &_vertices[i] );
                _highest_vertex = std::max ( _highest_vertex, _vertices[i] . label );
            }
        }

        bool is_any_active ( const U delta ) const noexcept
        {
            const auto delta_half = delta / 2;
            bool res = false;
            omp_set_num_threads ( _max_thread_count );
            #pragma omp parallel for schedule(static) reduction(||:res)
            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
                if ( _vertices[i] . excess > delta_half && _vertices[i] . label < _residual_network . size () &&
                     i != _sink )
                    res = true;
            omp_set_num_threads ( _thread_count );
            return res;
        }

        void set_original_labels ( ) noexcept
        {
            omp_set_num_threads ( static_cast<int> ( _max_thread_count ) );
            #pragma omp parallel for schedule(static)
            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
                _vertices[i] . original_label = _vertices[i] . label;
        }

        struct thread_local_data
        {
            const U delta;
            int64_t & cpu_time;
            const T low;
            const T high;
            T lowest_active;
            T & highest_active;
            T & highest_vertex;
            T relabel_progress;
        };


        void parallel_phase ( const U delta )
        {
            for ( ;; )
            {
                _relabel_threshold = _original_relabel_threshold;
                const auto partitions = partitioning::get_partitions ( _labels, _highest_active, _thread_count, min_active_per_thread );
                const auto actual_thread_cnt = partitions . size () - 1;
                omp_set_num_threads ( static_cast<int> ( actual_thread_cnt ) );

                int64_t cpu_time = 0;
                if ( actual_thread_cnt == 1 )
                {
                    push_relabel ( thread_local_data { delta, cpu_time, 0,
                                                       static_cast<T> ( _residual_network . size () ),
                                                       1, _highest_active, _highest_vertex, 0 } );
                    _thread_count = std::min ( _thread_count * 2, _max_thread_count );
                    return;
                }

                T highest_active = 0, highest_vertex = 0;
                #pragma omp parallel for schedule(static) reduction(+:cpu_time) reduction(max:highest_active) reduction(max:highest_vertex)
                for ( std::size_t i = 0; i < actual_thread_cnt; ++i )
                {
                    T low = partitions[i], high = partitions[i + 1];
                    highest_active = highest_vertex = high - 1;
                    push_relabel ( thread_local_data { delta, cpu_time, low, high, low, highest_active, highest_vertex,
                                                       0 } );
                    _relabel_threshold -= _original_relabel_threshold / actual_thread_cnt;
                    //add back vertices that are still active but couldn't have been relabeled to higher partition
                    _labels[high -
                            1] . active_vertices . append_list ( _thread_local_labels[omp_get_thread_num ()] . active_vertices );
                    if ( !_labels[high - 1] . active_vertices . empty () )
                        highest_active = highest_vertex = high - 1;
                }
                _highest_active = highest_active;
                _highest_vertex = std::max ( _highest_vertex, highest_vertex );

                if ( cpu_time > _min_cpu_time_per_phase )
                {
                    _thread_count = std::min ( _thread_count * 2, _max_thread_count );
                    return;
                }

                _thread_count = std::max ( actual_thread_cnt / 2, std::size_t { 1 } );
                set_original_labels ();
            }
        }

        void push_relabel ( thread_local_data data ) noexcept
        {
            auto start = std::chrono::high_resolution_clock::now ();
            while ( data . lowest_active <= data . highest_vertex )
            {
                if ( _labels[data . lowest_active] . active_vertices . empty () || data . lowest_active == data . low )
                {
                    ++data . lowest_active;
                    continue;
                }

                auto vertex = get_vertex_idx ( _labels[data . lowest_active] . active_vertices . front () );
                process ( vertex, data . lowest_active, data );

                if ( data . relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold )
                    break;
            }
            auto end = std::chrono::high_resolution_clock::now ();
            data . cpu_time += std::chrono::duration_cast<std::chrono::milliseconds> ( end - start ) . count ();
        }


        inline T get_vertex_idx ( vertex * n ) const noexcept
        {
            return std::distance ( _vertices . get (), n );
        }

        inline void process ( const T vertex, T label, thread_local_data & data ) noexcept
        {
            if ( push ( vertex, label, data ) )
                return;
            relabel ( vertex, label, data );
        }


        //original labels have to be set before the parallel phase starts and they don't change until the next one
        inline bool same_thread ( T original_label, T low, T high ) const noexcept
        {
            return original_label >= low && original_label < high;
        }


        inline bool push ( const T vertex, const T label, thread_local_data & data ) noexcept
        {
            const auto target_label = label - 1;
            for ( auto & edge : _residual_network[vertex] )
            {
                if ( edge . r_capacity > 0 &&
                     same_thread ( _vertices[edge . dst_vertex] . original_label, data . low, data . high )
                     && _vertices[edge . dst_vertex] . label == target_label )
                {
                    const auto old_target_excess = _vertices[edge . dst_vertex] . excess;

                    auto flow = std::min ( _vertices[vertex] . excess, edge . r_capacity );
                    if ( edge . dst_vertex != _sink && target_label != data . low )
                        flow = std::min ( flow, data . delta - _vertices[edge . dst_vertex] . excess );
                    _vertices[vertex] . excess -= flow;
                    _vertices[edge . dst_vertex] . excess += flow;
                    edge . r_capacity -= flow;
                    edge . reverse_r_capacity += flow;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . reverse_r_capacity -= flow;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . r_capacity += flow;

                    bool ret = false;
                    if ( _vertices[vertex] . excess <= data . delta / 2 )
                    {
                        _labels[label] . active_vertices . remove ( &_vertices[vertex] );
                        _labels[label] . inactive_vertices . push ( &_vertices[vertex] );
                        ret = true;
                    }

                    if ( edge . dst_vertex == _source || edge . dst_vertex == _sink )
                    {
                        if ( ret ) return true;
                        else continue;
                    }

                    //since we don't process vertices in data . low, we must ensure that we don't push them multiple times
                    if ( _vertices[edge . dst_vertex] . excess > data . delta / 2 &&
                         old_target_excess <= data . delta / 2 )
                    {
                        _labels[target_label] . inactive_vertices . remove ( &_vertices[edge . dst_vertex] );
                        _labels[target_label] . active_vertices . push ( &_vertices[edge . dst_vertex] );
                        --data . lowest_active;
                        ret = true;
                    }
                    if ( ret )
                        return true;
                }
            }
            return false;
        }


        inline void relabel ( const T vertex, const T current_label, thread_local_data & data ) noexcept
        {
            data . relabel_progress += BETA;
            const auto new_label = calculate_new_label ( vertex, data );
            _labels[current_label] . active_vertices . remove ( &_vertices[vertex] );
            _vertices[vertex] . label = std::min ( new_label, data . high - 1 );
            if ( new_label < data . high )
            {
                data . highest_vertex = std::max ( data . highest_vertex, new_label );
                data . highest_active = std::max ( data . highest_active, new_label );
                _labels[new_label] . active_vertices . push ( &_vertices[vertex] );
            } else if ( data . high < _residual_network . size () )
                //vertex is still active, but we cannot relabel it to another partition, so we remember it and add it back to active vertices at the end of this phase
                _thread_local_labels[omp_get_thread_num ()] . active_vertices . push ( &_vertices[vertex] );


            if ( _labels[current_label] . active_vertices . empty () &&
                 _labels[current_label] . inactive_vertices . empty () &&
                 current_label != data . high - 1 )
            {
                gap_relabel ( current_label, data );
                _vertices[vertex] . label = _residual_network . size ();
            }
        }


        inline T calculate_new_label ( const T vertex, thread_local_data & data ) noexcept
        {
            T increase_to = data . high - 1;
            for ( auto & edge :  _residual_network[vertex] )
            {
                if ( edge . r_capacity == 0 ||
                     !same_thread ( _vertices[edge . dst_vertex] . original_label, data . low, data . high ) )
                    continue;
                increase_to = std::min ( increase_to, _vertices[edge . dst_vertex] . label );
            }
            data . relabel_progress += _residual_network[vertex] . size ();
            return increase_to + 1;
        }


        void global_relabel ( const U delta ) noexcept
        {
            auto start = std::chrono::high_resolution_clock::now ();
            omp_set_num_threads ( static_cast<int> ( _max_thread_count ) );
            const auto not_reached = _residual_network . size ();

            #pragma omp parallel for schedule(static)
            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
            {
                _vertices[i] . discovered . clear ( std::memory_order_relaxed );
                _vertices[i] . label = _vertices[i] . original_label = not_reached;
            }

            #pragma omp parallel for schedule(static)
            for ( std::size_t i = 0; i <= _highest_vertex; ++i )
                _labels[i] . reset ();

            _vertices[_sink] . label = _vertices[_sink] . original_label = 0;
            _vertices[_sink] . discovered . test_and_set ( std::memory_order_relaxed );
            _highest_active = 0;

            _q[0] = _sink;
            std::size_t current_queue_size = 1;
            T current_distance = 0;
            const auto delta_half = delta / 2;

            while ( current_queue_size > 0 )
            {
                #pragma omp parallel for schedule(static)
                for ( std::size_t i = 0; i < current_queue_size; ++i )
                {
                    auto thr_id = omp_get_thread_num ();
                    auto current_vertex = _q[i];

                    for ( auto edge : _residual_network[current_vertex] )
                    {
                        if ( edge . reverse_r_capacity > 0 )
                        {
                            if ( !_vertices[edge . dst_vertex] . discovered . test_and_set ( std::memory_order_relaxed ) )
                            {
                                _vertices[edge . dst_vertex] . label = current_distance + 1;
                                _vertices[edge . dst_vertex] . original_label = current_distance + 1;
                                _pool . push_back ( edge . dst_vertex, static_cast<std::size_t>(thr_id) );

                                auto * node = &_vertices[edge . dst_vertex];
                                if ( _vertices[edge . dst_vertex] . excess > delta_half )
                                    _thread_local_labels[thr_id] . active_vertices . push ( node );
                                else
                                    _thread_local_labels[thr_id] . inactive_vertices . push ( node );
                            }
                        }
                    }
                }

                current_queue_size = _pool . swap_data ( _q );
                ++current_distance;

                for ( std::size_t i = 0; i < _max_thread_count; ++i ) //append together all thread_local info
                {
                    _labels[current_distance] . active_vertices . append_list ( _thread_local_labels[i] . active_vertices );
                    _labels[current_distance] . inactive_vertices . append_list (
                            _thread_local_labels[i] . inactive_vertices );
                }

                if ( !_labels[current_distance] . active_vertices . empty () )
                    _highest_active = current_distance;
            }

            _highest_vertex = current_distance - 1;

            omp_set_num_threads ( static_cast<int> ( _thread_count ) );
            auto end = std::chrono::high_resolution_clock::now ();

            _min_cpu_time_per_phase = std::chrono::duration_cast<std::chrono::milliseconds> ( end - start ) . count ();
        }


        //gap heuristic restricted to single segment
        void gap_relabel ( const T gap_height, const thread_local_data & data ) noexcept
        {
            for ( auto current_height = gap_height + 1; current_height <= data . highest_vertex; ++current_height )
            {
                while ( !_labels[current_height] . active_vertices . empty () )
                {
                    auto * ptr = _labels[current_height] . active_vertices . pop ();
                    auto vertex_idx = get_vertex_idx ( ptr );
                    _vertices[vertex_idx] . label = _residual_network . size ();
                }
                while ( !_labels[current_height] . inactive_vertices . empty () )
                {
                    auto * ptr = _labels[current_height] . inactive_vertices . pop ();
                    auto vertex_idx = get_vertex_idx ( ptr );
                    _vertices[vertex_idx] . label = _residual_network . size ();
                }
            }
            data . highest_vertex = data . highest_active = std::max ( gap_height - 1, data . low );
        }
    };
}

#endif //MAXFLOW_AHUJA_ORLIN_SEGMENT_H
