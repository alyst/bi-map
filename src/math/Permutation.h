#pragma once

#include "../BasicTypedefs.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>

/**
    C++ wrapper for \a gsl_permutation API.
 */
class Permutation {
private:
    gsl_permutation*    _data;

public:
    Permutation( size_t size )
        : _data( gsl_permutation_calloc( size ) )
    {
    }

    Permutation( const Permutation& a )
        : _data( gsl_permutation_alloc( a.size() ) )
    {
        if ( _data ) {
            gsl_permutation_memcpy( _data, a._data );
        }
    }

    void swap( Permutation& a )
    {
        std::swap( _data, a._data );
    }

    Permutation& operator=( const Permutation& a )
    {
        if ( a.size() == size() ) {
            // don't reallocate permutation, if the same size
            gsl_permutation_memcpy( _data, a._data );
        } else {
            Permutation aCopy( a );
            swap( aCopy );
        }
        return ( *this );
    }

    void init()
    {
        gsl_permutation_init( _data );
    }

    size_t size() const {
        return ( gsl_permutation_size( _data ) );
    }
    
    size_t operator[]( size_t index ) const {
        return ( gsl_permutation_get( _data, index ) );
    }
    
    void shuffle( const gsl_rng* r )
    {
        gsl_ran_shuffle( r, _data->data, gsl_permutation_size( _data ), sizeof( size_t ) );
    }

    ~Permutation()
    {
        gsl_permutation_free( _data );
    }
    
    static Permutation Random( const gsl_rng* r, size_t size )
    {
        Permutation res( size );
        res.shuffle( r );
        return ( res );
    }
};
