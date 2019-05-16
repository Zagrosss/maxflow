//
// Created by Jan Groschaft on 12/2/18.
//

/*
 * Implementation closely follows pseudocode descriped in the original paper:
 * K. Ahuja, Ravindra and Orlin, James, A Fast and Simple Algorithm for the Maximum Flow Problem, Operations Research, 1989
 */

#ifndef MAXFLOW_AHUJA_ORLIN_H
#define MAXFLOW_AHUJA_ORLIN_H


#include "../../common_types.h"
#include "../../data_structures/queue.h"
#include "../../data_structures/linked_list.h"
#include <memory>
#include <chrono>
#include <cmath>

namespace ahuja_orlin
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


    private:
        using pair = std::pair<T, T>;
        vector<vector<cached_edge<T, U>>> _residual_network;
        std::unique_ptr<label_info[]> _labels;
        std::unique_ptr<vertex[]> _vertices;
        data_structures::queue<pair> _distance_q;
        T _source, _sink, _lowest_active { 0 }, _highest_active { 0 }, _highest_vertex { 0 }, _relabel_progress {
                0 }, _relabel_threshold;
        U _max_cap;

        //statistics
        uint64_t _push_cnt { 0 }, _relabel_cnt { 0 }, _gap_cnt { 0 }, _gap_nodes { 0 }, _global_relabel_cnt { 0 };
    public:
        max_flow_instance ( vector<vector<cached_edge<T, U>>> graph, T source, T sink )
                :
                _residual_network ( std::move ( graph ) ),
                _labels ( std::make_unique<label_info[]> ( _residual_network . size () + 1 ) ),
                _vertices ( std::make_unique<vertex[]> ( _residual_network . size () ) ),
                _distance_q ( data_structures::queue<pair> { _residual_network . size () } ),
                _source ( source ), _sink ( sink )
        {
            init ();
        }

        U find_max_flow ( ) noexcept
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

        void init ( ) noexcept
        {
            _max_cap = 0;
            for ( auto & edge : _residual_network[_source] )
            {
                _max_cap = std::max ( _max_cap, edge . r_capacity );
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
            auto K = static_cast<U> ( std::ceil ( std::log2 ( _max_cap ) ));
            global_relabel ( _max_cap );
            for ( U k = 0; k <= K; ++k )
            {
                auto delta = static_cast<U> ( std::pow ( 2, K - k ));

                _lowest_active = _residual_network . size ();
                _highest_active = 1;

                for ( std::size_t i = 0; i <= _highest_vertex; ++i )
                {
                    _labels[i] . active_vertices . clear ();
                    _labels[i] . inactive_vertices . clear ();
                }

                for ( std::size_t i = 0; i < _residual_network . size (); ++i )
                {
                    if ( _vertices[i] . excess > delta / 2 && i != _source &&
                         i != _sink && _vertices[i] . label < _residual_network . size () )
                    {
                        _lowest_active = std::min ( _lowest_active, _vertices[i] . label );
                        _highest_active = std::max ( _highest_active, _vertices[i] . label );
                        _labels[_vertices[i] . label] . active_vertices . push ( &_vertices[i] );
                    } else
                        _labels[_vertices[i] . label] . inactive_vertices . push ( &_vertices[i] );
                }

                while ( _lowest_active <= _highest_active )
                {
                    if ( _labels[_lowest_active] . active_vertices . empty () )
                    {
                        ++_lowest_active;
                        continue;
                    }

                    auto vertex = get_vertex_idx ( _labels[_lowest_active] . active_vertices . front () );
                    process ( vertex, delta );

                    if ( _relabel_progress * GLOBAL_RELABEL_FREQ >= _relabel_threshold )
                    {
                        _relabel_progress = 0;
                        global_relabel ( delta );
                    }
                }
            }
        }


        T get_vertex_idx ( vertex * n ) noexcept
        {
            return std::distance ( _vertices . get (), n );
        }


        inline void process ( const T vertex, const U delta ) noexcept
        {
            const auto label = _vertices[vertex] . label;
            if ( push ( vertex, label, delta ) )
                return;
            relabel ( vertex, label );
        }


        inline bool push ( const T vertex, const T label, const U delta ) noexcept
        {
            for ( auto & edge : _residual_network[vertex] )
            {
                if ( edge . r_capacity > 0 && label == _vertices[edge . dst_vertex] . label + 1 )
                {
                    ++_push_cnt;
                    auto flow = std::min ( _vertices[vertex] . excess, edge . r_capacity );
                    if ( edge . dst_vertex != _sink )
                        flow = std::min ( flow, delta - _vertices[edge . dst_vertex] . excess );

                    _vertices[vertex] . excess -= flow;
                    _vertices[edge . dst_vertex] . excess += flow;
                    edge . r_capacity -= flow;
                    edge . reverse_r_capacity += flow;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . reverse_r_capacity -= flow;
                    _residual_network[edge . dst_vertex][edge . reverse_edge_index] . r_capacity += flow;

                    bool ret = false;
                    if ( _vertices[vertex] . excess <= delta / 2 )
                    {
                        _labels[label] . active_vertices . remove ( &_vertices[vertex] );
                        _labels[label] . inactive_vertices . push ( &_vertices[vertex] );
                        ret = true;
                    }

                    if ( _vertices[edge . dst_vertex] . excess > delta / 2 && edge . dst_vertex != _source &&
                         edge . dst_vertex != _sink )
                    {
                        _labels[label - 1] . inactive_vertices . remove ( &_vertices[edge . dst_vertex] );
                        _labels[label - 1] . active_vertices . push ( &_vertices[edge . dst_vertex] );
                        --_lowest_active;
                        ret = true;
                    }
                    if ( ret )
                        return true;
                }
            }
            return false;
        }


        inline void relabel ( const T vertex, const T current_label ) noexcept
        {
            ++_relabel_cnt;
            _relabel_progress += BETA;
            auto new_label = calculate_new_label ( vertex );
            _labels[current_label] . active_vertices . remove ( &_vertices[vertex] );
            _vertices[vertex] . label = new_label;

            if ( new_label != _residual_network . size () )
            {
                _highest_vertex = std::max ( _highest_vertex, new_label );
                _highest_active = std::max ( _highest_active, new_label );
                _labels[new_label] . active_vertices . push ( &_vertices[vertex] );
            }

            if ( _labels[current_label] . active_vertices . empty () &&
                 _labels[current_label] . inactive_vertices . empty () )
            {
                gap_relabel ( current_label );
                _vertices[vertex] . label = _residual_network . size ();
            }
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
            _relabel_progress += _residual_network[vertex] . size ();
            return increase_to + 1;
        }

        void global_relabel ( const U delta ) noexcept
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

            _lowest_active = _residual_network . size ();
            _highest_active = _highest_vertex = 1;

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
                        if ( edge . dst_vertex != _source )
                        {
                            auto * node = &_vertices[edge . dst_vertex];
                            if ( _vertices[edge . dst_vertex] . excess > delta / 2 )
                            {
                                _lowest_active = std::min ( _lowest_active, _vertices[edge . dst_vertex] . label );
                                _highest_active = std::max ( _highest_active, _vertices[edge . dst_vertex] . label );
                                _labels[_vertices[edge . dst_vertex] . label] . active_vertices . push ( node );
                            } else
                                _labels[_vertices[edge . dst_vertex] . label] . inactive_vertices . push ( node );
                        }
                    }
                }
            }
        }

        void gap_relabel ( const T gap_height ) noexcept
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

#endif //MAXFLOW_AHUJA_ORLIN_H
