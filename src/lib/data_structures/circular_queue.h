//
// Created by Jan Groschaft on 26.3.19.
//

/*
 * Circular buffer, fixed size.
 */

#ifndef MAXFLOW_CIRCULAR_QUEUE_H
#define MAXFLOW_CIRCULAR_QUEUE_H

#include <memory>

namespace data_structures
{
    template <typename T>
    class circular_queue
    {
        std::unique_ptr<T[]> _data { nullptr };
        std::size_t _front { 0 };
        std::size_t _back { 0 };
        const std::size_t _size { 0 };
        bool _flag { false };
    public:
        circular_queue ( ) = default;

        explicit circular_queue ( std::size_t size ) : _data ( std::make_unique<T[]> ( size ) ), _size ( size )
        { }

        circular_queue ( circular_queue<T> && other ) noexcept = default;

        circular_queue<T> & operator = ( circular_queue<T> && other ) noexcept = default;

        circular_queue ( const circular_queue<T> & other ) = delete;

        circular_queue<T> & operator = ( const circular_queue<T> & other ) = delete;

        void push ( T val ) noexcept
        {
            _data[_back++] = val;
            if ( _back == _size )
                _back = 0;
            _flag = true;
        }

        T pop ( ) noexcept
        {
            auto ret = _data[_front++];
            if ( _front == _size )
                _front = 0;
            _flag = false;
            return ret;
        }

        void reset ( ) noexcept
        { _front = _back = 0; }

        bool empty ( ) const noexcept
        { return _front == _back && !_flag; }
    };
}


#endif //MAXFLOW_CIRCULAR_QUEUE_H
