#pragma once

#include "BasicTypedefs.h"

#include <boost/serialization/set.hpp>

#include "EntityIndexing.h"

template<class Element>
struct CollectionEntitySerializer;

template<class Element>
struct IndexedEntityCompare {
    typedef EntityIndexing<Element>                 element_indexing_type;
    typedef typename element_indexing_type::entity_pointer_type element_pointer_type;

    bool operator()( const element_pointer_type& a, const element_pointer_type& b ) const
    {
        return ( a->serial() < b->serial() );
    }
};

template<class Element>
struct EntityIndexingTraits< std::set<typename IndexedEntity<Element>::pointer_type, 
                                      IndexedEntityCompare<Element> > >
{
    typedef CollectionEntitySerializer<Element> entity_serializer_type;
};

template<class Element>
class EntityCollectionIndexing: protected EntityIndexing< std::set<typename IndexedEntity<Element>::pointer_type, 
                                                                   IndexedEntityCompare<Element> > >
{
public:
    typedef EntityIndexing<Element>                 element_indexing_type;
    typedef typename element_indexing_type::entity_pointer_type element_pointer_type;
    typedef std::set<element_pointer_type, IndexedEntityCompare<Element> >          collection_type;

private:
    typedef EntityIndexing<collection_type>         super;
    typedef EntityCollectionIndexing<Element>       this_type;

    element_indexing_type           _elementIndexing;

    friend class boost::serialization::access;
    friend class CollectionEntitySerializer<Element>;

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << boost::serialization::make_nvp( "elementsIndex", _elementIndexing );
        ar << boost::serialization::make_nvp( "setsIndex", boost::serialization::base_object<super>( const_cast<this_type&>( *this ) ) );
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> boost::serialization::make_nvp( "elementsIndex", _elementIndexing );
        ar >> boost::serialization::make_nvp( "setsIndex", boost::serialization::base_object<super>( *this ) );
    }

public:
    using super::const_value_iterator;
    using super::const_serial_iterator;
    typedef typename super::size_type size_type;
    typedef typename super::entity_pointer_type collection_pointer_type;
    typedef typename super::indexed_entity_type indexed_collection_type;
   
    EntityCollectionIndexing( size_type min_size /*= MinIndexSize*/, size_type max_size /*= MaxIndexSize*/ )
    : super( min_size, max_size )
    , _elementIndexing( min_size, max_size )
    {
    }
    ~EntityCollectionIndexing()
    {
    }

    using super::iterator_to;
    using super::index;
    using super::valueMap;
    using super::serialMap;

    template<class Iterator>
    collection_pointer_type index( const Iterator& begin, const Iterator& end )
    {
        collection_type     coll;
        for ( Iterator it = begin; it != end; ++it ) {
            coll.insert( _elementIndexing.index( it->items() ) );
        }
        return index( coll );
    }

    const element_indexing_type& elementIndexing() const {
        return ( _elementIndexing );
    }

    element_indexing_type& elementIndexing() {
        return ( _elementIndexing );
    }

    /**
     *  Removes collections lacking outside references.
     */
    void remove_unreferenced() {
        super::remove_unreferenced();
        _elementIndexing.remove_unreferenced();
    }
};

#if 0
struct BielementsByIndexCompare {
    typedef BielementsIndexing::entity_pointer_type bielemnts_pointer_type;
    
    bool operator()( const bielemnts_pointer_type& a, const bielemnts_pointer_type& b ) const
    {
        return ( a->serial() < b->serial() );
    }
};
#endif

template<class Element>
struct CollectionEntitySerializer {
    typedef EntityCollectionIndexing<Element> entity_indexing_type;
    typedef typename entity_indexing_type::super base_entity_indexing_type;
    typedef typename entity_indexing_type::value_type entity_type;

    entity_indexing_type& indexing;
    unsigned int entity_version;

    CollectionEntitySerializer( const base_entity_indexing_type& indexing, unsigned int entity_version )
    : indexing( const_cast<entity_indexing_type&>( static_cast<const entity_indexing_type&>( indexing ) ) )
    , entity_version( entity_version )
    {
    }

    template<class Archive>
    void save_entity(Archive & ar, const entity_type& v) const
    {
        boost::serialization::collection_size_type elems_count( v.size() );
        ar << BOOST_SERIALIZATION_NVP( elems_count );
        for ( typename entity_type::const_iterator it = v.begin(); it != v.end(); ++it ) {
            default_serial_t elem_serial = (*it)->serial();
            ar << BOOST_SERIALIZATION_NVP( elem_serial );
        }
    }

    template<class Archive>
    entity_type load_entity(Archive & ar) const
    {
        boost::serialization::collection_size_type elems_count;
        ar >> BOOST_SERIALIZATION_NVP( elems_count );
        LOG_DEBUG3("Loaded elements count: " << elems_count );
        entity_type v;
        while ( elems_count-- ) {
            boost::serialization::collection_size_type elem_serial;
            ar >> BOOST_SERIALIZATION_NVP( elem_serial );
            LOG_DEBUG3( "Loading elements serial " << elem_serial << ", " << elems_count << " left" );
            typename entity_indexing_type::element_indexing_type::const_serial_iterator elemIt = indexing.elementIndexing().iterator_to( elem_serial );
            BOOST_ASSERT( elemIt != indexing.elementIndexing().serial_not_found() );
            BOOST_VERIFY( v.insert( *elemIt ).second );
        }
        return ( v ); 
    }
};
