//
// Created by Jan Groschaft on 26.3.18.
//

#ifndef EDMONDS_KARP_H
#define EDMONDS_KARP_H

#include <queue>
#include <limits>
#include <cstring>
#include <memory>
#include "../../data_structures/queue.h"
#include "../../common_types.h"

namespace edmonds_karp
{
    template <template <class> typename vector, typename T, typename U>
    class max_flow_instance
    {
        vector<vector<basic_edge<T, U>>> _residual_network;
        std::unique_ptr<basic_edge<T, U> * []> _parents;
        std::unique_ptr<bool[]> _visited;
        data_structures::queue<T> _q;
        T _source, _sink;

    public:
        max_flow_instance ( vector<vector<basic_edge<T, U>>> graph, T source, T sink ) :
                _residual_network ( std::move ( graph ) ),
                _parents ( std::make_unique<basic_edge<T, U> * []> ( _residual_network . size () ) ),
                _visited ( std::make_unique<bool[]> ( _residual_network . size () ) ),
                _q ( data_structures::queue<T> ( _residual_network . size () ) ),
                _source ( source ), _sink ( sink )
        {

        }

        U find_max_flow ( )
        {
            U max_flow = 0;
            while ( auto flow = find_augmenting_path () )
                max_flow += flow;
            return max_flow;
        }

        auto steal_network ( )
        {
            return std::move ( _residual_network );
        }

    private:
        U get_minimal_capacity_and_adjust_flow ( )
        {
            //get flow value
            auto current = _sink;
            auto flow = std::numeric_limits<U>::max ();
            while ( current != _source )
            {
                flow = std::min ( _parents[current] -> r_capacity, flow );
                current = _residual_network[_parents[current] -> dst_vertex][_parents[current] -> reverse_edge_index] . dst_vertex;
            }
            //update capacities
            current = _sink;
            while ( current != _source )
            {
                _parents[current] -> r_capacity -= flow;
                _residual_network[_parents[current] -> dst_vertex][_parents[current] -> reverse_edge_index] . r_capacity += flow;
                current = _residual_network[_parents[current] -> dst_vertex][_parents[current] -> reverse_edge_index] . dst_vertex;
            }
            return flow;
        }

        U find_augmenting_path ( )
        {
            std::memset ( _visited . get (), 0, _residual_network . size () );
            _q . reset ();
            _q . push ( _source );
            _visited[_source] = true;

            while ( !_q . empty () )
            {
                auto current = _q . pop ();

                if ( current == _sink )
                    return get_minimal_capacity_and_adjust_flow ();

                for ( auto & edge : _residual_network[current] )
                {
                    if ( edge . r_capacity > 0 && !_visited[edge . dst_vertex] )
                    {
                        _visited[edge . dst_vertex] = true;
                        _parents[edge . dst_vertex] = &edge;
                        _q . push ( edge . dst_vertex );
                    }
                }
            }
            return 0;
        }
    };
}

#endif //EDMONDS_KARP_H
