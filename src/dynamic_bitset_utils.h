#pragma once

#include "BasicTypedefs.h"

#include <boost/dynamic_bitset.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_free.hpp>

#define foreach_bit( iterator_type, iterator_name, bitset ) \
    for ( iterator_type iterator_name = (bitset).find_first(); iterator_name != boost::dynamic_bitset<>::npos; iterator_name = (bitset).find_next( iterator_name ) )

BOOST_SERIALIZATION_SPLIT_FREE( boost::dynamic_bitset<> )

namespace boost {
    namespace serialization {
        template <typename Archive> 
        void save(Archive& ar, const boost::dynamic_bitset<>& bitset, const unsigned int) {
            std::string bits;
            boost::to_string( bitset, bits );
            ar << BOOST_SERIALIZATION_NVP( bits );
        }

        template <typename Archive> 
        void load(Archive& ar, boost::dynamic_bitset<>& bitset, const unsigned int) {
            std::string bits;
            ar >> BOOST_SERIALIZATION_NVP( bits );
            bitset = boost::dynamic_bitset<>( bits );
        } 
    }
}

/**
 *  Simplest const_iterator for dynamic_bitset
 */
class dynamic_bitset_const_iterator {
private:
    typedef boost::dynamic_bitset<> bitset;
    typedef size_t index_type;

    const bitset*   _iterable;
    index_type      _ix;

    dynamic_bitset_const_iterator( const bitset& iterable, index_type ix )
    : _iterable( &iterable ), _ix( ix )
    {};

public:
    dynamic_bitset_const_iterator operator++() {
        if ( _ix != bitset::npos ) {
            _ix = _iterable->find_next( _ix );
        }
        return ( *this );
    }

    dynamic_bitset_const_iterator operator++(int) const {
        dynamic_bitset_const_iterator  copy( *this );
        return ( ++copy );
    }

    const index_type* const operator->() const {
        return ( &_ix );
    }

    index_type operator*() const {
        return ( _ix );
    }

    bool operator==( const dynamic_bitset_const_iterator& that ) const {
        if ( _iterable != that._iterable ) throw std::invalid_argument( "Incompatible iterator" );
        return ( _ix == that._ix );
    }

    bool operator!=( const dynamic_bitset_const_iterator& that ) const {
        if ( _iterable != that._iterable ) throw std::invalid_argument( "Incompatible iterator" );
        return ( _ix != that._ix );
    }

    static dynamic_bitset_const_iterator Begin( const bitset& iterable )
    {
        return ( dynamic_bitset_const_iterator( iterable, iterable.find_first() ) );
    }

    static dynamic_bitset_const_iterator End( const bitset& iterable )
    {
        return ( dynamic_bitset_const_iterator( iterable, bitset::npos ) );
    }
};

/**
 *  View (section) of the bitset.
 */
class dynamic_bitset_view {
private:
    const boost::dynamic_bitset<>&    _bitset;
    size_t  _offset;
    size_t  _stride;
    size_t  _size;

public:
    dynamic_bitset_view(
        const boost::dynamic_bitset<>&    bitset,
        size_t  offset,
        size_t  stride,
        size_t  size
    ) : _bitset( bitset )
    , _offset( offset ), _stride( stride )
    , _size( size )
    {}

    dynamic_bitset_view( const boost::dynamic_bitset<>&  bitset )
    : _bitset( bitset ), _offset( 0 ), _stride( 1 ), _size( bitset.size() )
    {}

    size_t size() const {
        return ( _size );
    }

    bool test( size_t index ) const {
        BOOST_ASSERT( index < _size );
        return ( _bitset.test( _offset + index * _stride ) );
    }
};
