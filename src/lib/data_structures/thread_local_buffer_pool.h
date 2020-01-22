//
// Created by Jan Groschaft on 2/22/19.
//

/*
 * Buffer which provides concurrent access to its underlying container, without the need of synchronization for
 * every insert. Each thread stores data to its own, smaller buffer, which is later copied to the shared thread_safe_buffer.
 */

#ifndef MAXFLOW_THREAD_LOCAL_BUFFER_POOL_H
#define MAXFLOW_THREAD_LOCAL_BUFFER_POOL_H

#include <vector>
#include "thread_safe_buffer.h"

#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif

namespace data_structures
{
    template <typename T>
    class thread_local_buffer_pool
    {
        struct alignas (CACHE_LINE_SIZE) aligned_vector
        {
            std::vector<T> data;
        };

        thread_safe_buffer<T> _buffer;
        std::unique_ptr<aligned_vector[]> _pool;
        const std::size_t _thread_count;
    public:
        thread_local_buffer_pool ( std::size_t thread_count, std::size_t total_max_size ) : _buffer ( total_max_size ),
                                                                                            _pool ( std::make_unique<aligned_vector[]> ( thread_count ) ),
                                                                                            _thread_count ( thread_count )
        {
            for ( std::size_t i = 0; i < _thread_count; ++i )
                _pool[i] . data . reserve ( std::max ( 1ul << 4, total_max_size / _thread_count ) );
        }

        void push_back ( const T & value, std::size_t thread_id ) noexcept
        {
            auto & target_buffer = _pool[thread_id] . data;
            if ( target_buffer . size () == target_buffer . capacity () )
            {
                _buffer . append ( target_buffer . data (), target_buffer . size () );
                target_buffer . clear ();
            }
            _pool[thread_id] . data . push_back ( value );
        }

        bool empty ( ) const noexcept
        {
            if ( !_buffer . empty () )
                return false;
            for ( std::size_t i = 0; i < _thread_count; ++i )
                if ( !_pool[i] . data . empty () )
                    return false;
            return true;
        }

        auto swap_data ( std::unique_ptr<T[]> & other ) noexcept
        {
            #pragma omp parallel for schedule(static)
            for ( std::size_t i = 0; i < _thread_count; ++i )
            {
                auto & target_buffer = _pool[i] . data;
                if ( target_buffer . size () > 0 )
                    _buffer . append ( target_buffer . data (), target_buffer . size () );
                target_buffer . clear ();
            }
            auto prev_size = _buffer . size ();
            _buffer . swap_data ( other );
            return prev_size;
        }
    };

}


#endif //MAXFLOW_THREAD_LOCAL_BUFFER_POOL_H
