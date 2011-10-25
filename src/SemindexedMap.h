#pragma once

#include "BasicTypedefs.h"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/mem_fun.hpp>
#include <boost/multi_index/composite_key.hpp>
#include <boost/multi_index/indexed_by.hpp>

#include "EntityIndexing.h"

template<class Key1, class Key2, class Value, class Key2Distance>
class SemindexedMap;

template<class Key1, class Key2, class Value, class Key2Distance>
class SemindexedMap {
public:
    typedef Key1 key1_type;
    typedef Key2 key2_type;
    typedef Value value_type;
    typedef EntityIndexing<key2_type> key2_indexing_type;
    typedef typename EntityIndexing<key2_type>::entity_pointer_type key2_pointer_type;

private:
    typedef IndexedEntity<key2_type> key2_indexed_type;
    typedef typename key2_indexed_type::serial_type key2_serial_type;
    typedef typename Key2Distance::distance_type distance_type;

    struct SemindexedEntry {
    public:
        typedef Key1 key1_type;
        typedef Key2 key2_type;
        typedef Value value_type;
        typedef EntityIndexing<key2_type> key2_indexing_type;
    private:
        typedef IndexedEntity<key2_type> key2_indexed_type;
        typedef typename key2_indexed_type::serial_type key2_serial_type;
        typedef typename key2_indexed_type::pointer_type key2_pointer_type;

        key1_type           _key1;
        key2_pointer_type   _key2;
        value_type          _value;

        friend class SemindexedMap<Key1, Key2, Value, Key2Distance>;

    public:
        SemindexedEntry( const key1_type& key1, const key2_pointer_type& key2, const value_type& value )
        : _key1( key1 ), _key2( key2 ), _value( value )
        {}

        const key2_serial_type key2serial() const {
            return ( _key2->serial() );
        }
        const value_type& value() const {
            return ( _value );
        }
        void setValue( const value_type& value ) {
            _value = value;
        }
    };

public:
    typedef SemindexedEntry indexed_entity_type;
    
private:
    typedef typename boost::multi_index::member<indexed_entity_type, const key1_type, &indexed_entity_type::_key1> key1_accessor_type;
    typedef typename boost::multi_index::const_mem_fun<indexed_entity_type, const key2_serial_type, &indexed_entity_type::key2serial> key2_accessor_type;

    typedef typename boost::multi_index::multi_index_container<
        indexed_entity_type,
        boost::multi_index::indexed_by<
            // by key1-key2
            boost::multi_index::ordered_unique<boost::multi_index::composite_key<
                indexed_entity_type, 
                key1_accessor_type, key2_accessor_type> >,
            // by key2
            boost::multi_index::hashed_non_unique<key2_accessor_type>
        > > semindex_type;

    typedef typename semindex_type::template nth_index<0>::type compound_index_type;
    typedef typename semindex_type::template nth_index<1>::type key2_index_type;
    typedef typename key2_index_type::const_iterator const_key2_iterator;
    typedef typename boost::tuples::tuple<key1_type, key2_serial_type> compound_key_type;

    key2_indexing_type&     _key2indexing;
    semindex_type           _semindex;

public:
    typedef typename compound_index_type::iterator       compound_iterator;
    typedef typename compound_index_type::const_iterator const_compound_iterator;
    typedef typename std::pair<const_compound_iterator, const_compound_iterator> key1_equal_range_type;
    typedef typename std::pair<compound_iterator, bool>  insert_result;

protected:
    static compound_key_type compoundKey( const key1_type& key1, const key2_serial_type& key2serial ) {
        return boost::make_tuple( key1, key2serial );
    }
    compound_key_type compoundKey( const key1_type& key1, const key2_type& key2 ) const {
        return compoundKey( key1, _key2indexing.index( key2 )->serial() );
    }
    indexed_entity_type entity( const key1_type& key1, const key2_pointer_type& pKey2, const value_type& value ) const {
        return ( indexed_entity_type( key1, pKey2, value ) );
    }
    
public:
    SemindexedMap( key2_indexing_type& key2indexing )
    : _key2indexing( key2indexing )
    {
    }
    
    key2_pointer_type key2index( const key2_type& key2 ) const {
        return ( _key2indexing.index( key2 ) );
    }
    
    const_compound_iterator find( const key1_type& key1, const key2_type& key2 ) const {
        return ( _semindex.find( compoundKey( key1, key2 ) ) );
    }
    const_compound_iterator find( const key1_type& key1, const key2_serial_type& key2serial ) const {
        return ( _semindex.find( compoundKey( key1, key2serial ) ) );
    }
    const_compound_iterator find( const key1_type& key1, const key2_pointer_type& pKey2 ) const {
        return ( _semindex.find( compoundKey( key1, pKey2->serial() ) ) );
    }
    
    const_compound_iterator findClosest( const key1_type& key1, const key2_type& key2, const distance_type& threshold ) const {
        if ( threshold <= 0 ) return ( end() );
        const key1_equal_range_type key1er = _semindex.equal_range( key1 );
        bool hasMatch = false;
        const_compound_iterator matchIt = _semindex.end();
        distance_type minDist = threshold;
        Key2Distance  key2dist;
        for ( const_compound_iterator it = key1er.first; it != key1er.second; ++it ) {
            const distance_type dist = key2dist( *it->_key2, key2, minDist );
            if ( ( dist < minDist ) || ( dist == minDist && !hasMatch ) ) {
                if ( dist == 0 ) return ( it );
                minDist = dist;
                matchIt = it;
                hasMatch = true;
            }
        }
        return ( matchIt );
    }
    const_compound_iterator findClosest( const key1_type& key1, const key2_serial_type& key2serial, const distance_type threshold ) const {
        return ( findClosest( key1, *_key2indexing.entity( key2serial ) ), threshold );
    }
    const_compound_iterator findClosest( const key1_type& key1, const key2_pointer_type& pKey2, const distance_type threshold ) const {
        return ( findClosest( key1, *pKey2 , threshold ) );
    }

    const_compound_iterator end() const {
        return ( _semindex.end() );
    }
    
    bool put( const key1_type& key1, const key2_type& key2, const value_type& value ) {
        key2_pointer_type pKey2 = key2index( key2 );
        return ( put( key1, pKey2, value ) );
    }
    bool put( const key1_type& key1, const key2_serial_type& key2Serial, const value_type& value ) {
        key2_pointer_type pKey2 = _key2indexing.entity( key2Serial );
        return ( put( key1, pKey2, value ) );
    }
    bool put( const key1_type& key1, const key2_pointer_type& pKey2, const value_type& value ) {
        indexed_entity_type entity( key1, pKey2, value );
        insert_result res = _semindex.insert( indexed_entity_type( key1, pKey2, value ) );
        if ( !res.second ) {
            return ( _semindex.replace( res.first, entity ) );
        } else {
            return ( true );
        }
    }
};
