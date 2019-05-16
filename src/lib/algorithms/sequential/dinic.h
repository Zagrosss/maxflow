//
// Created by Jan Groschaft on 21.10.18.
//

#ifndef DINIC_H
#define DINIC_H


#include <cstdint>
#include <limits>
#include <cstring>
#include <memory>
#include "../../data_structures/queue.h"
#include "../../common_types.h"

namespace dinic
{
    template <template <class> typename vector, typename T, typename U>
    class max_flow_instance
    {
        vector<vector<basic_edge<T, U>>> _residual_network;
        std::unique_ptr<T[]> _levels, _progress;
        data_structures::queue<T> _q;
        T _source, _sink;
    public:
        max_flow_instance ( vector<vector<basic_edge<T, U>>> graph, T source, T sink ) :
                _residual_network ( std::move ( graph ) ),
                _levels ( std::make_unique<T[]> ( _residual_network . size () ) ),
                _progress ( std::make_unique<T[]> ( _residual_network . size () ) ),
                _q ( data_structures::queue<T> ( _residual_network . size () ) ),
                _source ( source ), _sink ( sink )
        {

        }

        U find_max_flow ( )
        {
            U max_flow = 0;
            while ( build_level_graph () )
            {
                std::memset ( _progress . get (), 0, _residual_network . size () * sizeof ( T ) );
                while ( auto flow = augment_paths ( _source, std::numeric_limits<U>::max () ) )
                    max_flow += flow;
            }
            return max_flow;
        }

        auto steal_network ( )
        {
            return std::move ( _residual_network );
        }

    private:
        bool build_level_graph ( )
        {
            auto not_reached = std::numeric_limits<T>::max ();
            std::fill_n ( _levels . get (), _residual_network . size (), not_reached );
            _q . reset ();
            _q . push ( _source );
            _levels[_source] = 0;

            while ( !_q . empty () )
            {
                auto current = _q . pop ();
                if ( current == _sink ) return true;
                for ( auto edge : _residual_network[current] )
                {
                    if ( edge . r_capacity > 0 && _levels[edge . dst_vertex] == not_reached )
                    {
                        _levels[edge . dst_vertex] = _levels[current] + 1;
                        _q . push ( edge . dst_vertex );
                    }
                }
            }
            return false;
        }


        U augment_paths ( T current, U flow )
        {
            if ( current == _sink ) return flow;
            for ( ; _progress[current] < _residual_network[current] . size (); ++_progress[current] )
            {
                auto & edge = _residual_network[current][_progress[current]];
                if ( edge . r_capacity > 0 && _levels[edge . dst_vertex] == _levels[current] + 1 )
                {
                    auto result_flow = augment_paths ( edge . dst_vertex, std::min ( flow, edge . r_capacity ) );
                    if ( result_flow )
                    {
                        edge . r_capacity -= result_flow;
                        _residual_network[edge . dst_vertex][edge . reverse_edge_index] . r_capacity += result_flow;
                        return result_flow;
                    }
                }
            }
            return 0;
        }

    };
}
#endif //DINIC_H
