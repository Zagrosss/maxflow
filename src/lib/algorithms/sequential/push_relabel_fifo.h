//
// Created by Jan Groschaft on 16.11.18.
//

/*
 * Push-relabel, FIFO active vertex selection.
 */

#ifndef MAXFLOW_PUSH_RELABEL_FIFO_H
#define MAXFLOW_PUSH_RELABEL_FIFO_H

#include "../../common_types.h"
#include "../../data_structures/queue.h"
#include "../../data_structures/linked_list.h"
#include "../../data_structures/circular_queue.h"
#include <memory>
#include <iostream>
#include <queue>
#include <chrono>

namespace push_relabel_fifo
{
    template <template <class> typename vector, typename T, typename U>
    class max_flow_instance
    {
        struct vertex
        {
            U excess { 0 };
            T label;
        };
        using pair = std::pair<T, T>;
        vector<vector<cached_edge<T, U>>> _residual_network;
        std::unique_ptr<vertex[]> _vertices;
        data_structures::circular_queue<T> _q;
        data_structures::queue<pair> _distance_q;
        T _source, _sink, _relabel_progress { 0 }, _relabel_threshold;

        //statistics
        uint64_t _push_cnt { 0 }, _relabel_cnt { 0 }, _global_relabel_cnt { 0 };
    public:
        max_flow_instance ( vector<vector<cached_edge<T, U>>> graph, T source, T sink )
                :
                _residual_network ( std::move ( graph ) ),
                _vertices ( std::make_unique<vertex[]> ( _residual_network . size () ) ),
                _q ( data_structures::circular_queue<T> { _residual_network . size () } ),
                _distance_q ( data_structures::queue<pair> { _residual_network . size () } ),
                _source ( source ), _sink ( sink )
        {
            init ();
        }

        U find_max_flow ( )
        {
            find_max_flow_inner ();
            #ifdef DEBUG
            std::cout << "pushes:\t\t" << _push_cnt << std::endl;
            std::cout << "relabels:\t" << _relabel_cnt << std::endl;
            std::cout << "global updates:\t" << _global_relabel_cnt << std::endl;
            #endif
            return _vertices[_sink] . excess;
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


        void init ( )
        {
            for ( auto & edge : _residual_network[_source] )
            {
                if ( edge . r_capacity > 0 )
                {
                    _vertices[edge . dst_vertex] . excess = edge . r_capacity;
                    edge . reverse_r_capacity += edge . r_capacity;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . r_capacity += edge . r_capacity;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . reverse_r_capacity -= edge . r_capacity;
                    edge . r_capacity = 0;
                    ++_push_cnt;
                }
            }

            T m = 0;
            for ( auto & vec : _residual_network )
                m += vec . size ();
            _relabel_threshold = _residual_network . size () * ALPHA + m / 2;
        }

        void find_max_flow_inner ( )
        {
            global_relabel ();
            for ( ;; )
            {
                if ( _q . empty () )
                    return;

                auto vertex = _q . pop ();
                auto label = _vertices[vertex] . label;
                discharge ( vertex, label );

                if ( _relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold )
                {
                    _relabel_progress = 0;
                    global_relabel ();
                }
            }

        }


        void discharge ( const T vertex, T label )
        {
            while ( label < _residual_network . size () )
            {
                if ( push ( vertex, label ) )
                    return;
                label = relabel ( vertex );
            }
        }


        bool push ( const T vertex, const T label )
        {
            const auto target_label = label - 1;
            for ( auto & edge : _residual_network[vertex] )
            {
                if ( edge . r_capacity > 0 && _vertices[edge . dst_vertex] . label == target_label )
                {
                    ++_push_cnt;
                    auto flow = std::min ( _vertices[vertex] . excess, edge . r_capacity );
                    if ( _vertices[edge . dst_vertex] . excess == 0 && edge . dst_vertex != _sink )
                        _q . push ( edge . dst_vertex );
                    _vertices[vertex] . excess -= flow;
                    _vertices[edge . dst_vertex] . excess += flow;
                    edge . r_capacity -= flow;
                    edge . reverse_r_capacity += flow;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . reverse_r_capacity -= flow;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . r_capacity += flow;
                    if ( _vertices[vertex] . excess == 0 )
                        return true;
                }
            }
            return false;
        }


        T relabel ( const T vertex )
        {
            ++_relabel_cnt;
            _relabel_progress += BETA;
            _vertices[vertex] . label = calculate_new_label ( vertex );
            return _vertices[vertex] . label;
        }


        T calculate_new_label ( const T vertex )
        {
            T increase_to = _residual_network . size () - 1;
            for ( auto & edge : _residual_network[vertex] )
            {
                if ( edge . r_capacity == 0 )
                    continue;
                increase_to = std::min ( increase_to, _vertices[edge . dst_vertex] . label );
            }
            _relabel_progress += _residual_network[vertex] . size ();
            return increase_to + 1;
        }


        void global_relabel ( )
        {
            ++_global_relabel_cnt;
            auto not_reached = _residual_network . size ();

            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
                _vertices[i] . label = not_reached;

            _q . reset ();
            _distance_q . reset ();
            _distance_q . push ( std::make_pair ( _sink, 0 ) );
            _vertices[_sink] . label = 0;

            while ( !_distance_q . empty () )
            {
                auto current_elem = _distance_q . pop ();
                auto current_vertex = current_elem . first;
                auto current_distance = current_elem . second;
                for ( auto & edge : _residual_network[current_vertex] )
                {
                    if ( edge . reverse_r_capacity > 0 && _vertices[edge . dst_vertex] . label == not_reached )
                    {
                        _vertices[edge . dst_vertex] . label = current_distance + 1;
                        _distance_q . push ( std::make_pair ( edge . dst_vertex, current_distance + 1 ) );
                        if ( _vertices[edge . dst_vertex] . excess > 0 )
                            _q . push ( edge . dst_vertex );
                    }
                }
            }
        }
    };
}


#endif //MAXFLOW_PUSH_RELABEL_FIFO_H
