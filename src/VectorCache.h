#pragma once

#include "BasicTypedefs.h"

#include <stack>
#include <vector>

/**
 * Cache of vector objects.
 */
template<class T>
class VectorCache {
private:
    typedef std::vector<T> vector_type;
    std::stack<vector_type*> cache;

public:
    ~VectorCache()
    {
        while ( !cache.empty() ) {
            delete cache.top();
            cache.pop();
        }
    }

    /**
     *  Borrow empty vector from the cache.
     */
    vector_type& borrow( size_t size, const T& def = T() )
    {
        if ( cache.empty() ) {
            return ( *(new vector_type( size )) );
        } else {
            vector_type* res = cache.top();
            res->resize( size, def );
            cache.pop();
            return ( *res );
        }
    }

    /**
     *  Borrow vector from the cache,
     *  which is a copy of a.
     */
    vector_type& borrow( const vector_type& a )
    {
        if ( cache.empty() ) {
            return ( *(new vector_type( a )) );
        } else {
            vector_type* res = cache.top();
            *res = a;
            cache.pop();
            return ( *res );
        }
    }

    /**
     *  Return vector to the cache.
     */
    void redeem( vector_type& v )
    {
        BOOST_ASSERT( cache.empty() || ( &v != cache.top() ) );
        cache.push( &v );
    }
};
