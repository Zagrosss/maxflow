//
// Created by Jan Groschaft on 12/11/18.
//

/*
 * Implementation of Goldberg-Tarjan's parallel push-relabel algorithm. Description can be found in
 * Goldberg, Andrew and Tarjan, Robert, A New Approach to the Maximum-Flow Problem, J. ACM, 1988.
 *
 * This implementation is also based on detailed pseudocode presented in
 * Baumstark, Niklas, Speeding up Maximum Flow Computations on Shared-Memory Platforms, KIT, Karlsruhe, 2014.
 */

#ifndef MAXFLOW_PARALLEL_PUSH_RELABEL_H
#define MAXFLOW_PARALLEL_PUSH_RELABEL_H

#include <memory>
#include <chrono>
#include <iostream>
#include <atomic>
#include <omp.h>
#include <algorithm>
#include "../../common_types.h"
#include "../../data_structures/queue.h"
#include "../../data_structures/thread_local_buffer_pool.h"

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

namespace parallel_push_relabel
{
    template <template <class> typename vector, typename T, typename U>
    class max_flow_instance
    {
        struct alignas (CACHE_LINE_SIZE) vertex
        {
            U excess { 0 };
            std::atomic<U> new_excess { 0 };
            T label;
            T new_label;
            std::atomic_flag discovered = ATOMIC_FLAG_INIT;
        };

        vector<vector<cached_edge<T, U>>> _residual_network;
        std::unique_ptr<vertex[]> _vertices;
        std::unique_ptr<T[]> _active { };
        data_structures::thread_local_buffer_pool<T> _pool;
        T _source, _sink, _relabel_threshold, _active_cnt;
        std::size_t _relabel_progress;
        const T _thread_count;
    public:
        max_flow_instance ( vector<vector<cached_edge<T, U>>> graph, T source, T sink,
                            std::size_t thread_count = static_cast<size_t>(omp_get_max_threads ()) )
                :
                _residual_network ( std::move ( graph ) ),
                _vertices ( std::make_unique<vertex[]> ( _residual_network . size () ) ),
                _active ( std::make_unique<T[]> ( _residual_network . size () ) ),
                _pool ( data_structures::thread_local_buffer_pool<T> ( thread_count, _residual_network . size () ) ),
                _source ( source ), _sink ( sink ), _active_cnt ( 0 ), _relabel_progress ( 0 ),
                _thread_count ( thread_count )
        {
            omp_set_num_threads ( static_cast<int> ( _thread_count ) );
            init ();
        }

        uint64_t _phase_cnt = 0;
        uint64_t _push_cnt = 0;
        uint64_t _global_update_cnt = 0;

        U find_max_flow ( ) noexcept
        {
            find_max_flow_inner ();

            #ifdef DEBUG
            std::cout << "global updates:\t" << _global_update_cnt << std::endl;
            std::cout << "phase cnt: " << _phase_cnt << std::endl;
            std::cout << "pushes: " << _push_cnt << std::endl;
            #endif
            return _vertices[_sink] . new_excess + _vertices[_sink] . excess;
        }

        void preflow_to_flow ( )
        {
            std::swap ( _source, _sink );
            find_max_flow_inner ();
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
        static constexpr double GLOBAL_RELABEL_FREQ = 0.5;

        void init ( ) noexcept
        {
            #pragma omp parallel for schedule(static)
            for ( std::size_t i = 0; i < _residual_network[_source] . size (); ++i )
            {
                auto & edge = _residual_network[_source][i];
                _vertices[edge . dst_vertex] . excess = edge . r_capacity;
                edge . reverse_r_capacity += edge . r_capacity;
                _residual_network[edge . dst_vertex][edge . reverse_edge_index] . r_capacity += edge . r_capacity;
                _residual_network[edge . dst_vertex][edge . reverse_edge_index] . reverse_r_capacity -= edge . r_capacity;
                edge . r_capacity = 0;

            }

            T m = 0;
            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
                m += _residual_network[i] . size ();
            _relabel_threshold = _residual_network . size () * ALPHA + m / 2;
        }


        void find_max_flow_inner ( )
        {
            global_relabel ();

            for ( ;; )
            {
                if ( _active_cnt == 0 )
                    return;

                ++_phase_cnt;
                uint64_t push_cnt_per_phase = 0;

                #pragma omp parallel
                {
                    #pragma omp for schedule(static) reduction(+:push_cnt_per_phase)
                    for ( T i = 0; i < _active_cnt; ++i )
                    {
                        auto thr_id = omp_get_thread_num ();
                        auto vertex = _active[i];
                        if ( _vertices[vertex] . label == _residual_network . size () )
                            continue;
                        push ( vertex, _vertices[vertex] . label, thr_id, push_cnt_per_phase );
                    }
                    //stage 2
                    #pragma omp for schedule(static) reduction(+:_relabel_progress)
                    for ( T i = 0; i < _active_cnt; ++i )
                    {
                        auto thr_id = omp_get_thread_num ();
                        auto vertex = _active[i];
                        relabel ( vertex, thr_id, _relabel_progress );
                    }
                    //stage 3
                    #pragma omp for schedule(static)
                    for ( T i = 0; i < _active_cnt; ++i )
                    {
                        auto vertex = _active[i];
                        _vertices[vertex] . label = _vertices[vertex] . new_label;
                        _vertices[vertex] . discovered . clear ( std::memory_order_relaxed );
                    }
                    //stage 4
                    #pragma omp single
                    _active_cnt = _pool . swap_data ( _active );

                    #pragma omp for schedule(static)
                    for ( T i = 0; i < _active_cnt; ++i )
                    {
                        auto vertex = _active[i];
                        _vertices[vertex] . excess += _vertices[vertex] . new_excess . load ( std::memory_order_relaxed );
                        _vertices[vertex] . new_excess . store ( 0, std::memory_order_relaxed );
                        _vertices[vertex] . discovered . clear ( std::memory_order_relaxed );
                    }
                }

                if ( _relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold || push_cnt_per_phase == 0 )
                {
                    _relabel_progress = 0;
                    global_relabel ();
                }

                _push_cnt += push_cnt_per_phase;
            }
        }

        inline void push ( const T vertex, const T label, int thr_id, uint64_t & push_cnt ) noexcept
        {
            const auto target_label = label - 1;
            for ( auto & edge : _residual_network[vertex] )
            {
                if ( edge . r_capacity > 0 && _vertices[edge . dst_vertex] . label == target_label )
                {
                    auto flow = std::min ( _vertices[vertex] . excess, edge . r_capacity );
                    if ( edge . dst_vertex != _source && edge . dst_vertex != _sink )
                        if ( !_vertices[edge . dst_vertex] . discovered . test_and_set ( std::memory_order_relaxed ) )
                            _pool . push_back ( edge . dst_vertex, static_cast<size_t>(thr_id) );
                    ++push_cnt;
                    _vertices[vertex] . excess -= flow;
                    _vertices[edge . dst_vertex] . new_excess . fetch_add ( flow, std::memory_order_relaxed );
                    edge . r_capacity -= flow;
                    edge . reverse_r_capacity += flow;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . reverse_r_capacity -= flow;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . r_capacity += flow;
                    if ( _vertices[vertex] . excess == 0 )
                        return;
                }
            }
        }

        inline void relabel ( const T vertex, const int thr_id, std::size_t & relabel_progress ) noexcept
        {
            if ( _vertices[vertex] . excess > 0 || _vertices[vertex] . label == _residual_network . size () )
            {
                relabel_progress += BETA;
                _vertices[vertex] . new_label = calculate_new_label ( vertex );
                relabel_progress += _residual_network[vertex] . size ();
                if ( _vertices[vertex] . new_label == _residual_network . size () )
                {
                    _vertices[vertex] . excess += _vertices[vertex] . new_excess;
                    _vertices[vertex] . new_excess = 0;
                    return;
                }

                if ( !_vertices[vertex] . discovered . test_and_set ( std::memory_order_relaxed ) )
                    _pool . push_back ( vertex, static_cast<size_t>(thr_id) );
            } else
                _vertices[vertex] . new_label = _vertices[vertex] . label;
        }

        inline T calculate_new_label ( const T vertex ) noexcept
        {
            T increase_to = _residual_network . size () - 1;
            for ( auto & edge : _residual_network[vertex] )
            {
                if ( edge . r_capacity == 0 )
                    continue;

                increase_to = std::min ( increase_to, _vertices[edge . dst_vertex] . label );
            }
            return increase_to + 1;
        }


        void global_relabel ( ) noexcept
        {
            ++_global_update_cnt;
            const auto not_reached = _residual_network . size ();

            #pragma omp parallel for schedule(static)
            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
                _vertices[i] . label = not_reached;

            _vertices[_sink] . label = 0;
            _vertices[_sink] . discovered . test_and_set ();
            assert ( _pool . empty () );
            _active[0] = _sink;
            std::size_t current_queue_size = 1;
            T current_distance = 0;

            while ( current_queue_size > 0 )
            {
                #pragma omp parallel for schedule(static)
                for ( std::size_t i = 0; i < current_queue_size; ++i )
                {
                    auto thr_id = omp_get_thread_num ();
                    auto current_vertex = _active[i];

                    for ( auto edge : _residual_network[current_vertex] )
                    {
                        if ( edge . reverse_r_capacity > 0 )
                        {
                            if ( !_vertices[edge . dst_vertex] . discovered . test_and_set ( std::memory_order_relaxed ) )
                            {
                                _vertices[edge . dst_vertex] . label = current_distance + 1;
                                _pool . push_back ( edge . dst_vertex, static_cast<std::size_t>(thr_id) );
                            }
                        }
                    }
                }
                current_queue_size = _pool . swap_data ( _active );
                ++current_distance;
            }


            #pragma omp parallel for schedule(static)
            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
            {
                auto thr_id = omp_get_thread_num ();
                if ( _vertices[i] . label != not_reached && _vertices[i] . excess > 0 && i != _sink )
                    _pool . push_back ( i, static_cast<size_t>(thr_id) );
                _vertices[i] . discovered . clear ( std::memory_order_relaxed );
            }

            _active_cnt = _pool . swap_data ( _active );
        }
    };
}


#endif //MAXFLOW_PARALLEL_PUSH_RELABEL_H
