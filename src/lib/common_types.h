//
// Created by Jan Groschaft on 23.10.18.
//

#ifndef COMMON_TYPES_H
#define COMMON_TYPES_H

#include <cstdint>
#include <utility>
#include <limits>

template <typename T = uint32_t, typename U = uint32_t>
struct basic_edge
{
    basic_edge ( ) = default;

    basic_edge ( T dst_vertex, U capacity, T reverse_edge_index )
            : dst_vertex ( dst_vertex ), reverse_edge_index ( reverse_edge_index ), r_capacity ( capacity )
    { }

    T dst_vertex;
    T reverse_edge_index;
    U r_capacity;
};


template <typename T = uint32_t, typename U = uint32_t>
struct cached_edge
{
    cached_edge ( ) = default;

    cached_edge ( T dst_vertex, U capacity, T reverse_edge_index )
            : dst_vertex ( dst_vertex ), reverse_edge_index ( reverse_edge_index ), r_capacity ( capacity )
    { }

    T dst_vertex;
    T reverse_edge_index;
    U r_capacity;
    U reverse_r_capacity;
};


#endif //COMMON_TYPES_H
