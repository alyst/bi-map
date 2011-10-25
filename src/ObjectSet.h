#pragma once

#include "BasicTypedefs.h"

#include <boost/functional/hash.hpp>

typedef std::set<object_index_t>     object_set_t;

inline std::size_t hash_value( const object_set_t& value )
{
    return ( boost::hash_range( value.begin(), value.end() ) );
}

struct ObjectSetDistance {
    typedef size_t distance_type;

    distance_type operator()( const object_set_t& a, const object_set_t& b, distance_type threshold )
    {
        const object_set_t::const_iterator aEnd = a.end();
        const object_set_t::const_iterator bEnd = b.end();

        distance_type  res = 0;
        for ( object_set_t::const_iterator aIt = a.begin(); aIt != aEnd; ++aIt ) {
            if ( b.find( *aIt ) == bEnd ) {
                res++;
                if ( res > threshold ) return ( res );
            }
        }
        for ( object_set_t::const_iterator bIt = b.begin(); bIt != bEnd; ++bIt ) {
            if ( a.find( *bIt ) == aEnd ) {
                res++;
                if ( res > threshold ) return ( res );
            }
        }
        return ( res );
    }
};

