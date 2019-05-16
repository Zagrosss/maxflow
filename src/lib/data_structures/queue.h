//
// Created by Jan Groschaft on 21.10.18.
//

/*
 * Simple non-growing, non-circular queue. Can only enqueue given number of elements, after that it has to be reset.
 */

#ifndef FORD_FULKERSON_QUEUE_H
#define FORD_FULKERSON_QUEUE_H

#include <cstdint>
#include <memory>

namespace data_structures
{

    template <typename T>
    class queue
    {
        std::unique_ptr<T[]> _data { nullptr };
        std::size_t _back { 0 };
        std::size_t _front { 0 };
    public:
        queue ( ) = default;

        explicit queue ( std::size_t size ) : _data ( std::make_unique<T[]> ( size ) )
        { }

        queue ( queue<T> && other ) noexcept = default;

        queue<T> & operator = ( queue<T> && other ) noexcept = default;

        queue ( const queue<T> & other ) = delete;

        queue<T> & operator = ( const queue<T> & other ) = delete;

        void push ( T val ) noexcept
        { _data[_back++] = val; }

        T pop ( ) noexcept
        { return _data[_front++]; }

        bool empty ( ) const noexcept
        { return _front == _back; }

        void reset ( ) noexcept
        { _front = _back = 0; }
    };
}

#endif //FORD_FULKERSON_QUEUE_H
