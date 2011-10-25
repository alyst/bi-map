#pragma once

#include "BasicTypedefs.h"

#include <stdexcept>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/indexed_by.hpp>

#include <boost/intrusive_ptr.hpp>

#include <boost/serialization/version.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/base_object.hpp>

typedef unsigned int default_serial_t;

template<class ValueType, typename SerialType = default_serial_t>
class EntityIndexing;

template<class ValueType, typename SerialType = default_serial_t>
struct IndexedEntity {
    friend class boost::serialization::access;

public:
    typedef SerialType serial_type;
    typedef ValueType value_type;
    typedef IndexedEntity<value_type, serial_type> this_type;
    typedef typename boost::intrusive_ptr<this_type> pointer_type;

private:
    value_type  _value;
    serial_type _serial;
    mutable volatile std::size_t _refCount;

    typedef EntityIndexing<ValueType, SerialType> entity_indexing_type;
    friend class EntityIndexing<value_type, serial_type>;

#if defined(_DEBUG)
    const entity_indexing_type* _owner;
#endif

protected:

    IndexedEntity( const value_type& value, serial_type serial
#if defined(_DEBUG)
        , const entity_indexing_type* owner
#endif
        )
        : _value( value ), _serial( serial ), _refCount( 0 )
#if defined(_DEBUG)
        , _owner( owner )
#endif
    {
    }

    /// required by pointer_value_type for maintaining reference count
    friend void intrusive_ptr_add_ref( const IndexedEntity<value_type, serial_type>* pIEnt ) {
        pIEnt->_refCount++;
    }

    /// required by pointer_value_type for maintaining reference count
    friend void intrusive_ptr_release( const IndexedEntity<value_type, serial_type>* pIEnt ) {
        pIEnt->_refCount--;
        if ( pIEnt->_refCount == 0 ) {
            delete pIEnt;
            pIEnt = NULL;
        }
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "serial", _serial );
        ar & boost::serialization::make_nvp( "value", _value );
    }

 public:
    serial_type serial() const { return _serial; }
    const value_type& value() const { return _value; }
    operator const value_type&() const { return value(); }
    size_t refCount() const { return _refCount; }

    /// indexed value without owner
    static this_type orphan(const value_type& value) {
        this_type res( value, 0
#if defined(_DEBUG)
            , NULL 
#endif
            );
        ++res._refCount; // because referenced from current scope and would be autodeleted there
        return res;
    }

#if defined(_DEBUG)
    const entity_indexing_type* owner() const {
        return _owner;
    }
#endif
};

template<class ValueType>
struct IndexedEntityValueHasher
{
    std::size_t operator()( const ValueType& value ) const {
        return ( hash_value( value ) );
    }
};

template<class ValueType, typename SerialType>
struct DefaultEntitySerializer
{
    typedef EntityIndexing<ValueType, SerialType> entity_indexing_type;

    unsigned int entity_version;

    DefaultEntitySerializer( const entity_indexing_type& indexing, unsigned int entity_version )
    : entity_version( entity_version )
    {
    }

    template<class Archive>
    void save_entity(Archive & ar, const ValueType& entity) const
    {
        ar << boost::serialization::make_nvp( "entity", const_cast<ValueType&>( entity ) );
    }

    template<class Archive>
    ValueType load_entity(Archive & ar) const
    {
        ValueType entity;
        ar >> boost::serialization::make_nvp( "entity", entity );
        return ( entity );
    }
};

template<class ValueType, typename SerialType = default_serial_t>
struct EntityIndexingTraits
{
    typedef DefaultEntitySerializer<ValueType, SerialType> entity_serializer_type;
};

/// maintains collection
template<class ValueType, class SerialType>
class EntityIndexing {
public:
    friend class boost::serialization::access;

    typedef std::size_t size_type;
    typedef ValueType value_type;
    typedef SerialType serial_type;
    typedef IndexedEntity<value_type, serial_type> indexed_entity_type;
    typedef typename indexed_entity_type::pointer_type entity_pointer_type;

protected:
    typedef typename EntityIndexingTraits<value_type, serial_type>::entity_serializer_type entity_serializer_type;
    typedef typename boost::multi_index::multi_index_container<
        entity_pointer_type,
        boost::multi_index::indexed_by<
            // by assigned serial number
            boost::multi_index::hashed_unique<boost::multi_index::member<const indexed_entity_type, const SerialType, &indexed_entity_type::_serial> >,
            // by entity contents
            boost::multi_index::hashed_unique<boost::multi_index::member<const indexed_entity_type, const ValueType, &indexed_entity_type::_value>, IndexedEntityValueHasher<ValueType> >
        > > multindex_type;
    typedef EntityIndexing<value_type, serial_type> this_type; 
    typedef typename multindex_type::template nth_index<0>::type serial_index;
    typedef typename multindex_type::template nth_index<1>::type value_index;
    typedef typename value_index::iterator value_iterator;
    typedef typename serial_index::iterator serial_iterator;

public:
    typedef typename value_index::const_iterator const_value_iterator;
    typedef typename serial_index::const_iterator const_serial_iterator;
    
protected:

private:
    multindex_type  _index;
    /// max size after garbage collection
    const size_type _minSize;
    /// max allowed size before garbage collection
    const size_type _maxSize;

    volatile serial_type _nextSerial;

protected:
#if defined(_DEBUG)
    static const size_type MinIndexSize = 10000;
    static const size_type MaxIndexSize = 30000;
#else
    static const size_type MinIndexSize = 100000;
    static const size_type MaxIndexSize = 500000;
#endif

    value_index& valueMap() {
        return _index.template get<1>();
    }

    serial_index& serialMap() {
        return _index.template get<0>();
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        const unsigned int entity_version = boost::serialization::version<value_type>::value;
        ar << BOOST_SERIALIZATION_NVP(entity_version);
        serial_type nextSerial = _nextSerial;
        boost::serialization::collection_size_type count( valueMap().size() );
        ar << BOOST_SERIALIZATION_NVP( count );
        ar << BOOST_SERIALIZATION_NVP( nextSerial );

        entity_serializer_type serializer( *this, entity_version );
        for ( const_serial_iterator it = serialMap().begin(); it != serialMap().end(); ++it ) {
            serial_type serial = (*it)->serial();
            ar << BOOST_SERIALIZATION_NVP( serial );
            serializer.save_entity( ar, (*it)->value() );
        }
    }

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        unsigned int entity_version;
        boost::serialization::collection_size_type count;
        ar >> BOOST_SERIALIZATION_NVP(entity_version);
        ar >> BOOST_SERIALIZATION_NVP(count);
        boost::serialization::collection_size_type nextSerial;
        ar >> BOOST_SERIALIZATION_NVP(nextSerial);
        _nextSerial = nextSerial;

        _index = multindex_type();
        entity_serializer_type serializer( *this, entity_version );
        LOG_DEBUG3( count << " entities to read" );
        while ( count-- ) {
            boost::serialization::collection_size_type serial;
            ar >> BOOST_SERIALIZATION_NVP( serial );
            LOG_DEBUG3( "Read entity serial " << serial << ", " << count << " left" );
#if defined(_DEBUG)
            bool res = 
#endif
            serialMap().insert( new indexed_entity_type( serializer.load_entity( ar ), serial
#if defined(_DEBUG)
                , this
#endif
                ) ).second;
            BOOST_ASSERT( res );
        }
    }

public:
    EntityIndexing( size_type minSize = MinIndexSize, size_type maxSize = MaxIndexSize )
    : _minSize( minSize ), _maxSize( maxSize ), _nextSerial( 0 )
    {}

    ~EntityIndexing()
    {}

    const value_index& valueMap() const {
        return _index.template get<1>();
    }

    const serial_index& serialMap() const {
        return _index.template get<0>();
    }

    const const_value_iterator serialToValue( const const_serial_iterator& serIt ) const {
        return ( _index.template project<1>( serIt ) );
    }

    bool operator==( EntityIndexing& that ) const {
        /// @fixme use super class comparison
        return this == &that;
    }

    const_value_iterator value_not_found() const {
        return valueMap().end();
    }

    const_serial_iterator serial_not_found() const {
        return serialMap().end();
    }

    const_serial_iterator iterator_to( serial_type serial ) const {
        return serialMap().find( serial );
    }

    const_value_iterator iterator_to( const value_type& value ) const {
        return valueMap().find( value );
    }

    /// index submitted rowset or retrieve existing index
    entity_pointer_type index( const value_type& value ) {
        value_iterator eit = valueMap().find( value );
        if ( eit == value_not_found() ) {
            if ( is_full() ) {
                throw std::overflow_error("Entity Index is full");
                //collect_garbage();
            }

            // @TODO: exclusive lock
            value_iterator newit = valueMap().insert( eit, new indexed_entity_type( value, _nextSerial++
#if defined(_DEBUG)
                , this
#endif
                ) );
            if ( newit == eit ) throw std::runtime_error( "Indexing of entity failed" );
            return ( *newit );
        }
        else {
            return ( *eit );
        }
    }

    /// retrieve entity by its serial, NULL if not found
    entity_pointer_type entity( serial_type serial ) const {
        const_serial_iterator eit = serialMap().find( serial );
        return ( eit == serial_not_found() ? entity_pointer_type(NULL) : *eit );
    }

    size_t size() const {
        return ( _index.size() );
    }

    bool is_full() const {
        return ( size() > _maxSize );
    }

    /**
     *  Removes objects that have no outside references.
     */
    void remove_unreferenced() {
        for ( value_iterator vit = valueMap().begin(); vit != valueMap().end(); ) {
            if ( (*vit)->refCount() <= 1 ) {
                vit = valueMap().erase( vit );
            }
            else {
                ++vit;
            }
        }
    }

    // size_type collect_garbage();
};
