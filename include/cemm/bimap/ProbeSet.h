#pragma once

#include <cemm/containers/dynamic_bitset_hash.h>

namespace cemm { namespace bimap {

typedef boost::dynamic_bitset<>                  probe_bitset_t;

struct ProbeSetDistance {
    typedef size_t distance_type;

    distance_type operator()( const probe_bitset_t& a, const probe_bitset_t& b, distance_type threshold )
    {
        size_t  ab = ( a - b ).count();
        if ( ab > threshold ) return ( ab );
        return ( ab + ( b - a ).count() );
    }
};

} }