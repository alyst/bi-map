#pragma once

#include "BasicTypedefs.h"

#include <stdexcept>
#include <vector>
#include <boost/serialization/access.hpp>
#include <boost/serialization/nvp.hpp>

/**
 *  2D array.
 */
template<typename T>
class array2d {
private:
    typedef typename std::vector<T> container_type;
    size_t          _size1;
    size_t          _size2;
    container_type  _data;

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "size1", _size1 );
        ar & boost::serialization::make_nvp( "size2", _size2 );
        ar & boost::serialization::make_nvp( "data", _data );
    }

public:
    array2d( size_t size1 = 0, size_t size2 = 0, const T& val = T() )
    : _size1( size1 ), _size2( size2 )
    , _data( size1 * size2, val )
    {}

    void reset( size_t size1 = 0, size_t size2 = 0, const T& val = T() )
    {
        _size1 = size1;
        _size2 = size2;
        _data.assign( _size1 * _size2, val );
    }

    bool empty() const {
        return ( _size1 == 0 || _size2 == 0 );
    }

    size_t size1() const {
        return ( _size1 );
    }
    size_t size2() const {
        return ( _size2 );
    }

    const T& operator()( size_t i, size_t j ) const {
        BOOST_ASSERT( ( i < _size1 ) && ( j < _size2 ) );
        return ( _data[ i * _size2 + j ] );
    }
    T& operator()( size_t i, size_t j ) {
        BOOST_ASSERT( ( i < _size1 ) && ( j < _size2 ) );
        return ( _data[ i * _size2 + j ] );
    }
    const T* data( size_t i = 0, size_t j = 0 ) const {
        BOOST_ASSERT( ( i < _size1 ) && ( j < _size2 ) );
        return ( &_data[ i * _size2 + j ] );
    }
    T* data( size_t i = 0, size_t j = 0) {
        BOOST_ASSERT( ( i < _size1 ) && ( j < _size2 ) );
        return ( &_data[ i * _size2 + j ] );
    }

    void insert1( size_t ix1, const T& val = T() )
    {
        BOOST_ASSERT( _size1 >= ix1 );
        _data.resize( (_size1+1) * _size2, val );
        if ( !_data.empty() ) {
            const size_t offsDst = ( ix1 + 1 ) * _size2;
            const size_t offsSrc = ix1 * _size2;
            const size_t len = ( _size1 - ix1 ) * _size2;
            for ( int i = (int)len-1; i >= 0; i-- ) {
                _data[ offsDst + i ] = _data[ offsSrc + i ];
            }
            for ( size_t i = 0; i < _size2; i++ ) {
                _data[ offsSrc + i ] = val;
            }
        }
        _size1++;
    }

    void insert2( size_t ix2, const T& val = T() )
    {
        BOOST_ASSERT( size2() >= ix2 );
        _data.resize( _size1 * (_size2+1), val );
        if ( !_data.empty() ) {
            for ( int i = (int)_size1-1; i >= 0; i-- ) {
                size_t offsDst = i * ( _size2 + 1 ) + 1;
                const size_t offsSrc = i * _size2;
                for ( int j = (int)_size2-1; j >= (int)ix2; j-- ) {
                    _data[ offsDst + j ] = _data[ offsSrc + j ];
                }
                offsDst--;
                _data[ offsDst + ix2 ] = val;
                if ( i > 0 ) {
                    for ( int j = (int)ix2 - 1; j >= 0; j-- ) {
                        _data[ offsDst + j ] = _data[ offsSrc + j ];
                    }
                }
            }
        }
        _size2++;
    }

    void remove1( size_t ix1, const T& val = T() )
    {
        BOOST_ASSERT( ix1 < size1() );
        if ( !_data.empty() ) {
            const size_t offsDst = ix1 * _size2;
            const size_t offsSrc = ( ix1 + 1 ) * _size2;
            const size_t len = ( _size1 - ix1 - 1 ) * _size2;
            for ( size_t i = 0; i < len; i++ ) {
                _data[ offsDst + i ] = _data[ offsSrc + i ];
            }
            _data.resize( ( _size1 - 1 ) * _size2 );
        }
        _size1--;
    }
    void remove2( size_t ix2, const T& val = T() )
    {
        BOOST_ASSERT( ix2 < _size2 );
        if ( !_data.empty() ) {
            for ( size_t i = 0; i < _size1; i++ ) {
                size_t offsDst = i * ( _size2 - 1 );
                const size_t offsSrc = i * _size2;
                if ( i > 0 ) {
                    for ( size_t j = 0; j < ix2; j++ ) {
                        _data[ offsDst + j ] = _data[ offsSrc + j ];
                    }
                }
                offsDst--;
                for ( size_t j = ix2 + 1; j < _size2; j++ ) {
                    _data[ offsDst + j ] = _data[ offsSrc + j ];
                }
            }
            _data.resize( _size1 * (_size2-1) );
        }
        _size2--;
    }
    void fill( const T& val )
    {
        for ( size_t i = 0; i < _data.size(); i++ ) {
            _data[ i ] = val;
        }
    }
    void fill1( size_t ix1, const T& val )
    {
        BOOST_ASSERT( ix1 < _size1 );
        if ( !_data.empty() ) {
            const size_t lastIx = ( ix1 + 1 ) * _size2;
            for ( size_t i = ix1 * _size2; i < lastIx; i++ ) {
                _data[ i ] = val;
            }
        }
    }
    void fill2( size_t ix2, const T& val )
    {
        BOOST_ASSERT( ix2 < _size2 );
        for ( size_t i = ix2; i < _size1 * _size2; i+= _size2 ) {
            _data[ i ] = val;
        }
    }

    /**
     *  Section of array2d in selected dimension.
     */
    class Section {
    private:
        typedef array2d<T> array2d_type;
        friend class array2d<T>;

        array2d_type&   _array;
        size_t          _dim;   // 0-based dimension of the section
        size_t          _ix;    // fixed index of the sectioned dimension

        Section( array2d_type& array, size_t dim, size_t ix )
        : _array( array ), _dim( dim ), _ix( ix )
        {
            BOOST_ASSERT( ix < ( _dim == 0 ? _array.size1() : _array.size2() ) );
        }

    public:
        size_t size() const {
            return ( _dim == 0 ? _array.size2() : _array.size1() );
        }
        const T& operator[]( size_t i ) const {
            BOOST_ASSERT( i < size() );
            return ( _dim == 0 ? _array( _ix, i ) : _array( i, _ix ) );
        }
        T& operator[]( size_t i ) {
            BOOST_ASSERT( i < size() );
            return ( _dim == 0 ? _array( _ix, i ) : _array( i, _ix ) );
        }
        /**
         *  Copy the contents of one section to another.
         */
        Section& operator=( const Section& that )
        {
            // @todo: speed up
            if ( size() != that.size() ) throw std::runtime_error( "Incompatible section sizes" );
            for ( size_t i = 0; i < size(); i++ ) {
                operator[]( i ) = that[ i ];
            }
            return ( *this );
        }
        /**
         *  Copy the contents of vector to section.
         */
        template<typename Vector>
        Section& operator=( const Vector& v )
        {
            // @todo: speed up
            if ( size() != v.size() ) throw std::runtime_error( "Incompatible sizes of section and container" );
            for ( size_t i = 0; i < size(); i++ ) {
                operator[]( i ) = v[ i ];
            }
            return ( *this );
        }

        operator std::vector<T>() const
        {
            std::vector<T> res( size() );
            for ( size_t i = 0; i < size(); i++ ) {
                res[i] = operator[]( i );
            }
            return ( res );
        }
    };

    typedef Section section_type;

    section_type section1( size_t ix ) {
        return ( Section( *this, 0, ix ) );
    }
    const section_type section1( size_t ix ) const {
        return ( Section( const_cast< array2d< T >& >( *this ), 0, ix ) );
    }
    section_type section2( size_t ix ) {
        return ( Section( *this, 1, ix ) );
    }
    const section_type section2( size_t ix ) const {
        return ( Section( const_cast< array2d< T >& >( *this ), 1, ix ) );
    }
};
