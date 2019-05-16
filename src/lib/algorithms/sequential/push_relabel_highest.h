//
// Created by Jan Groschaft on 19.11.18.
//

/*
 * Push-relabel, highest active vertex selection.
 */

#ifndef MAXFLOW_PUSH_RELABEL_HIGHEST_H
#define MAXFLOW_PUSH_RELABEL_HIGHEST_H


#include "../../common_types.h"
#include "../../data_structures/queue.h"
#include "../../data_structures/linked_list.h"
#include <memory>
#include <queue>
#include <cassert>
#include <chrono>
#include <cstring>

namespace push_relabel_highest
{
    template <template <class> typename vector, typename T, typename U>
    class max_flow_instance
    {
        struct vertex
        {
            vertex * next;
            vertex * prev;
            U excess { 0 };
            T label;
        };

        struct label_info
        {
            data_structures::linked_list<vertex> active_vertices;
            data_structures::linked_list<vertex> inactive_vertices;
        };

        using pair = std::pair<T, T>;
        vector<vector<cached_edge<T, U>>> _residual_network;
        std::unique_ptr<label_info[]> _labels;
        std::unique_ptr<vertex[]> _vertices;
        data_structures::queue<pair> _distance_q;
        T _source, _sink, _highest_active, _highest_vertex, _relabel_progress, _relabel_threshold;

        //statistics
        uint64_t _push_cnt { 0 }, _relabel_cnt { 0 }, _gap_cnt { 0 }, _gap_nodes { 0 }, _global_relabel_cnt { 0 };
    public:
        max_flow_instance ( vector<vector<cached_edge<T, U>>> graph, T source, T sink )
                :
                _residual_network ( std::move ( graph ) ),
                _labels ( std::make_unique<label_info[]> ( _residual_network . size () + 1 ) ),
                _vertices ( std::make_unique<vertex[]> ( _residual_network . size () ) ),
                _distance_q ( data_structures::queue<pair> { _residual_network . size () } ),
                _source ( source ), _sink ( sink ), _relabel_progress ( 0 )
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
            std::cout << "gaps:\t\t" << _gap_cnt << std::endl;
            std::cout << "gap nodes:\t" << _gap_nodes << std::endl;
            #endif
            return _vertices[_sink] . excess;
        }

        void preflow_to_flow ( )
        {
            std::swap ( _source, _sink );
            _highest_vertex = _residual_network . size ();
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
                _vertices[edge . dst_vertex] . excess = edge . r_capacity;
                edge . reverse_r_capacity += edge . r_capacity;
                _residual_network[edge . dst_vertex][edge . reverse_edge_index] . r_capacity += edge . r_capacity;
                _residual_network[edge . dst_vertex][edge . reverse_edge_index] . reverse_r_capacity -= edge . r_capacity;
                edge . r_capacity = 0;
            }
            _push_cnt += _residual_network[_source] . size ();


            T m = 0;
            for ( std::size_t i = 0; i < _residual_network . size (); ++i )
                m += _residual_network[i] . size ();
            _relabel_threshold = _residual_network . size () * ALPHA + m / 2;
            _highest_vertex = 1;
        }


        void find_max_flow_inner ( )
        {
            global_relabel ();

            for ( ;; )
            {
                auto node = get_active_vertex ();
                auto label = _highest_active;
                if ( node == nullptr )
                    return;

                discharge ( node, label );

                if ( _relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold )
                {
                    _relabel_progress = 0;
                    global_relabel ();
                }
            }
        }

        inline auto get_active_vertex ( )
        {
            for ( T i = 0; i <= _highest_active; ++i )
            {
                if ( _labels[_highest_active - i] . active_vertices . empty () )
                    continue;
                auto * node = _labels[_highest_active - i] . active_vertices . pop ();
                _highest_active -= i;
                return node;
            }
            return static_cast<vertex *> (nullptr);
        }


        inline T get_vertex_idx ( vertex * n )
        {
            return std::distance ( _vertices . get (), n );
        }


        inline void discharge ( vertex * n, T label )
        {
            T vertex = get_vertex_idx ( n );
            for ( ;; )
            {
                if ( push ( vertex, label ) )
                {
                    _labels[label] . inactive_vertices . push ( n );
                    return;
                }
                auto new_label = relabel ( vertex, label );
                label = new_label;
                if ( label == _residual_network . size () )
                    return;
            }
        }


        inline bool push ( const T vertex, const T label )
        {
            const auto target_label = label - 1;
            for ( auto & edge : _residual_network[vertex] )
            {
                if ( edge . r_capacity > 0 && _vertices[edge . dst_vertex] . label == target_label )
                {
                    ++_push_cnt;
                    auto flow = std::min ( _vertices[vertex] . excess, edge . r_capacity );
                    if ( _vertices[edge . dst_vertex] . excess == 0 && edge . dst_vertex != _sink )
                    {
                        auto * node = &_vertices[edge . dst_vertex];
                        _labels[target_label] . inactive_vertices . remove ( node );
                        _labels[target_label] . active_vertices . push ( node );
                    }
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


        inline T relabel ( const T vertex, const T current_label )
        {
            ++_relabel_cnt;
            _relabel_progress += BETA;
            const auto new_label = calculate_new_label ( vertex );
            _vertices[vertex] . label = new_label;

            if ( new_label != _residual_network . size () )
            {
                _highest_vertex = std::max ( _highest_vertex, new_label );
                _highest_active = new_label - 1;
            }

            if ( _labels[current_label] . active_vertices . empty () &&
                 _labels[current_label] . inactive_vertices . empty () )
            {
                gap_relabel ( current_label );
                _vertices[vertex] . label = _residual_network . size ();
            }

            return _vertices[vertex] . label;
        }


        inline T calculate_new_label ( const T vertex )
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
            _distance_q . reset ();
            _distance_q . push ( std::make_pair ( _sink, 0 ) );
            _vertices[_sink] . label = 0;

            for ( std::size_t i = 0; i <= _highest_vertex; ++i )
            {
                _labels[i] . active_vertices . clear ();
                _labels[i] . inactive_vertices . clear ();
            }

            _highest_active = _highest_vertex = 0;
            while ( !_distance_q . empty () )
            {
                auto current_elem = _distance_q . pop ();
                auto current_vertex = current_elem . first;
                auto current_distance = current_elem . second;
                _highest_vertex = std::max ( _highest_vertex, current_distance );
                for ( auto & edge : _residual_network[current_vertex] )
                {
                    if ( edge . reverse_r_capacity > 0 && _vertices[edge . dst_vertex] . label == not_reached )
                    {
                        _vertices[edge . dst_vertex] . label = current_distance + 1;
                        _distance_q . push ( std::make_pair ( edge . dst_vertex, current_distance + 1 ) );
                        auto * node = &_vertices[edge . dst_vertex];
                        if ( _vertices[edge . dst_vertex] . excess > 0 )
                        {
                            _highest_active = std::max ( _highest_active, _vertices[edge . dst_vertex] . label );
                            _labels[current_distance + 1] . active_vertices . push ( node );
                        } else
                            _labels[current_distance + 1] . inactive_vertices . push ( node );
                    }
                }
            }
        }

        void gap_relabel ( const T gap_height )
        {
            ++_gap_cnt;
            for ( auto current_height = gap_height + 1; current_height <= _highest_vertex; ++current_height )
            {
                while ( !_labels[current_height] . active_vertices . empty () )
                {
                    ++_gap_nodes;
                    auto * ptr = _labels[current_height] . active_vertices . pop ();
                    auto vertex_idx = get_vertex_idx ( ptr );
                    _vertices[vertex_idx] . label = _residual_network . size ();
                }
                while ( !_labels[current_height] . inactive_vertices . empty () )
                {
                    ++_gap_nodes;
                    auto * ptr = _labels[current_height] . inactive_vertices . pop ();
                    auto vertex_idx = get_vertex_idx ( ptr );
                    _vertices[vertex_idx] . label = _residual_network . size ();
                }
            }
            _highest_vertex = _highest_active = gap_height - 1;
        }
    };
}
#endif //MAXFLOW_PUSH_RELABEL_HIGHEST_H