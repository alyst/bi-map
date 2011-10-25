#pragma once

#include "BasicTypedefs.h"

#include <vector>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/nvp.hpp>

template<typename T>
class symmetric_array2d {
private:
    typedef typename std::vector<T> container_type;
    size_t          _size;
    container_type  _data;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "size", _size );
        ar & boost::serialization::make_nvp( "data", _data );
    }

    static size_t Index( size_t size, size_t i, size_t j ) {
        BOOST_ASSERT( ( i < size ) && ( j < size ) );
        size_t min_ix = std::min( i, j );
        return ( ( 2 * size - ( min_ix - 1 ) ) * min_ix / 2 + ( std::max( i, j ) - min_ix ) );
    }
    size_t index( size_t i, size_t j ) const {
        return ( Index( _size, i, j ) );
    }

public:
    symmetric_array2d( size_t size = 0, const T& val = T() )
    : _size( size ), _data( size * ( size + 1 ) / 2, val )
    {}

    void reset( const T& val = T() )
    {
        _data.assign( _data.size(), val );
    }

    void resize( size_t size, const T& val = T() )
    {
        _size = size;
        _data.assign( _size * ( _size + 1 ) / 2, val );
    }

    bool empty() const {
        return ( _size == 0 );
    }

    size_t size() const {
        return ( _size );
    }

    const T& operator()( size_t i, size_t j ) const {
        return ( _data[ index( i, j ) ] );
    }
    T& operator()( size_t i, size_t j ) {
        return ( _data[ index( i, j ) ] );
    }
    const T* data( size_t i = 0, size_t j = 0 ) const {
        return ( &_data[ index( i, j ) ] );
    }
    T* data( size_t i = 0, size_t j = 0) {
        return ( &_data[ index( i, j ) ] );
    }

    void insert( size_t ix, const T& val = T() )
    {
        BOOST_ASSERT( _size >= ix );
        container_type oldData( _data );
        _size++;
        _data.resize( ( _size + 1 ) * _size / 2, val );
        if ( !oldData.empty() ) {
            size_t srcIx = 0;
            size_t destIx = 0;
            for ( size_t i = 0; i < _size; i++ ) {
                for ( size_t j = i; j < _size; j++ ) {
                    _data[ destIx++ ] = ( i == ix || j == ix ) 
                                      ? val : oldData[ srcIx++ ];
                }
            }
        }
    }

    void remove( size_t ix, const T& val = T() )
    {
        BOOST_ASSERT( ix < size() );
        _size--;
        if ( _size > 0 ) {
            container_type oldData( _data );
            _data.resize( ( _size + 1 ) * _size / 2, val );
            size_t srcIx = 0;
            size_t destIx = 0;
            for ( size_t i = 0; i < _size + 1; i++ ) {
                for ( size_t j = i; j < _size + 1; j++ ) {
                    if ( i != ix && j != ix ) {
                        _data[ destIx++ ] = oldData[ srcIx ];
                    }
                    srcIx++;
                }
            }
        }
    }
    void fill( const T& val )
    {
        for ( size_t i = 0; i < _data.size(); i++ ) {
            _data[ i ] = val;
        }
    }
    void fill( size_t ix, const T& val )
    {
        BOOST_ASSERT( ix < _size );
        if ( !_data.empty() ) {
            for ( size_t i = 0; i < _size; i++ ) {
                _data[ index( i, ix ) ] = val;
            }
        }
    }
    symmetric_array2d<T>& operator=( const symmetric_array2d<T>& that )
    {
        _data = that._data;
        _size = that._size;
        return ( *this );
    }
};
