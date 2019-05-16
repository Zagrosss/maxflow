//
// Created by Jan Groschaft on 20.11.18.
//

/*
 * Simple linked list using preallocated nodes.
 */

#ifndef MAXFLOW_LINKED_LIST_H
#define MAXFLOW_LINKED_LIST_H


#include <utility>
#include <iostream>
#include <cassert>

namespace data_structures
{
    template <typename node>
    class linked_list
    {
        node _head = node { }, _tail = node { };
        std::size_t _size { 0 };
    public:
        linked_list ( )
        {
            _head . next = &_tail;
            _tail . prev = &_head;
        }

        linked_list ( const linked_list & other ) = delete;

        linked_list ( linked_list && other ) = delete;

        linked_list & operator = ( const linked_list & other ) = delete;

        linked_list & operator = ( linked_list && other ) = delete;

        node * pop ( ) noexcept
        {
            auto * ret = _head . next;
            _head . next = _head . next -> next;
            _head . next -> prev = &_head;
            --_size;
            return ret;
        }

        void push ( node * n ) noexcept
        {
            n -> next = &_tail;
            n -> prev = _tail . prev;
            _tail . prev -> next = n;
            _tail . prev = n;
            ++_size;
        }

        void push_front ( node * n ) noexcept
        {
            n -> next = _head . next;
            n -> prev = &_head;
            _head . next -> prev = n;
            _head . next = n;
            ++_size;
        }

        void remove ( node * n ) noexcept
        {
            n -> prev -> next = n -> next;
            n -> next -> prev = n -> prev;
            --_size;
        }

        node * front ( ) const noexcept
        {
            return _head . next;
        }

        node * back ( ) const noexcept
        {
            return _tail . prev;
        }

        bool empty ( ) const noexcept
        {
            return _size == 0;
        }

        void clear ( ) noexcept
        {
            _head . next = &_tail;
            _tail . prev = &_head;
            _size = 0;
        }

        std::size_t size ( ) const noexcept
        {
            return _size;
        }

        void append_list ( linked_list & other ) noexcept
        {
            if ( other . empty () )
                return;
            auto other_head = other . front ();
            auto other_tail = other . back ();
            this -> back () -> next = other_head;
            other_head -> prev = this -> back ();
            _tail . prev = other_tail;
            other_tail -> next = &_tail;
            _size += other . size ();
            other . clear ();
        }
    };
}


#endif //MAXFLOW_LINKED_LIST_H
