//
// Created by Jan Groschaft on 2/22/19.
//

/*
 * Thread safe array, data is retrieved using the swap method.
 */

#ifndef MAXFLOW_THREAD_SAFE_BUFFER_H
#define MAXFLOW_THREAD_SAFE_BUFFER_H

#include <memory>
#include <atomic>
#include <cstring>

namespace data_structures
{
    template <typename T>
    class thread_safe_buffer
    {
        std::unique_ptr<T[]> _buffer;
        std::atomic<std::size_t> _pos;
    public:
        explicit thread_safe_buffer ( std::size_t size ) : _buffer ( std::make_unique<T[]> ( size ) ), _pos ( 0 )
        {
        }

        void clear ( ) const noexcept
        {
            _pos = 0;
        }

        bool empty ( ) const noexcept
        {
            return _pos == 0;
        }

        void append ( const T * const data, std::size_t size ) noexcept
        {
            auto begin = _pos . fetch_add ( size, std::memory_order_relaxed );
            std::memcpy ( _buffer . get () + begin, data, size * sizeof ( T ) );
        }

        void push_back ( const T & val ) noexcept
        {
            auto current = _pos . fetch_add ( 1, std::memory_order_relaxed );
            _buffer[current] = val;
        }

        std::size_t size ( ) const noexcept
        {
            return _pos . load ( std::memory_order_relaxed );
        }

        void swap_data ( std::unique_ptr<T[]> & other ) noexcept
        {
            using std::swap;
            _buffer . swap ( other );
            _pos . store ( 0, std::memory_order_relaxed );
        }
    };
}

#endif //MAXFLOW_THREAD_SAFE_BUFFER_H
