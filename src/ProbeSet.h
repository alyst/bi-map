#pragma once

#include "BasicTypedefs.h"

#include <boost/dynamic_bitset.hpp>
#include <boost/functional/hash.hpp>

#include "BasicTypedefs.h"

typedef boost::dynamic_bitset<>                  probe_bitset_t;

inline std::size_t hash_value( const probe_bitset_t& value ) {
    std::vector<long unsigned int> blocks;
    blocks.reserve( value.num_blocks() );
    boost::to_block_range( value, blocks.begin() );
    return ( boost::hash_value( blocks ) );
}

struct ProbeSetDistance {
    typedef size_t distance_type;
    
    distance_type operator()( const probe_bitset_t& a, const probe_bitset_t& b, distance_type threshold )
    {
        size_t  ab = ( a - b ).count();
        if ( ab > threshold ) return ( ab );
        return ( ab + ( b - a ).count() );
    }
};
