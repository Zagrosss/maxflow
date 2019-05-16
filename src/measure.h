//
// Created by Jan Groschaft on 28.10.18.
//

#ifndef MAXFLOW_MEASURE_H
#define MAXFLOW_MEASURE_H

#include <vector>
#include <iostream>
#include "common_types.h"
#include <chrono>
#include <omp.h>
#include <cassert>
#include <optional>

template <typename U>
struct measurement_result
{
    U max_flow;
    std::chrono::milliseconds time_read;
    std::chrono::milliseconds time_init;
    std::chrono::milliseconds time_solve;
};

template <typename S, typename U>
auto measure_single ( S & solver, measurement_result<U> & result )
{
    auto start = std::chrono::high_resolution_clock::now ();
    auto res = solver . find_max_flow ();
    auto end = std::chrono::high_resolution_clock::now ();
    result . max_flow = res;
    result . time_solve = std::chrono::duration_cast<std::chrono::milliseconds> ( end - start );
}


template <template <template <typename> typename, typename, typename> typename S,
        template <typename> typename vector, typename T, typename U,
        template <typename, typename> typename EDGE>
auto measure_parallel ( vector<vector<EDGE<T, U>>> graph, T source, T sink, std::size_t thread_count )
{
    measurement_result<U> measurement;
    auto start = std::chrono::high_resolution_clock::now ();
    S solver ( std::move ( graph ), source, sink, thread_count );
    auto end = std::chrono::high_resolution_clock::now ();
    measurement . time_init = std::chrono::duration_cast<std::chrono::milliseconds> ( end - start );
    measure_single ( solver, measurement );
    return measurement;
}

template <template <template <typename> typename, typename, typename> typename S,
        template <typename> typename vector, typename T, typename U,
        template <typename, typename> typename EDGE>
auto measure_sequential ( vector<vector<EDGE<T, U>>> graph, T source, T sink )
{
    measurement_result<U> measurement;
    auto start = std::chrono::high_resolution_clock::now ();
    S solver ( std::move ( graph ), source, sink );
    auto end = std::chrono::high_resolution_clock::now ();
    measurement . time_init = std::chrono::duration_cast<std::chrono::milliseconds> ( end - start );
    measure_single ( solver, measurement );
    return measurement;
}


#endif //MAXFLOW_MEASURE_H
