#pragma once

#include "BasicTypedefs.h"

#include <sstream>
#include <boost/bimap.hpp>
#include <boost/serialization/traits.hpp>

/**
 *  Wrapper around object to enable static tracking.
 */
template<class T>
struct statically_tracked: public boost::serialization::wrapper_traits<const statically_tracked<T> >
{
private:
    typedef T object_type;

    typedef unsigned int id_type;
    typedef boost::bimap<id_type, const object_type*> bimap_type;
    typedef typename bimap_type::value_type bimap_value_type;
    typedef std::pair<id_type, bool> register_result_type;

    const char* _name;
    const object_type*& _pT;
    static bimap_type tracking_bimap;

    static register_result_type register_object( const char* name, const object_type* obj )
    {
        typename bimap_type::right_const_iterator it = tracking_bimap.right.find( obj );
        if ( it == tracking_bimap.right.end() ) {
            // register new object
            id_type id = tracking_bimap.left.size();
            LOG_DEBUG2( "Registering " << name << " id=" << id << ", ptr=0x" << obj );
            tracking_bimap.insert( bimap_value_type( id, obj ) );
            return ( std::make_pair( id, true ) );
        }
        else {
            LOG_DEBUG2( "Already registered " << name << " id=" << it->second << ", ptr=0x" << obj );
            // already registered, return old id
            return ( std::make_pair( it->second, false ) );
        }
    }

    static register_result_type register_object_with_id( const char* name, const object_type* obj, id_type id )
    {
        typename bimap_type::right_const_iterator it = tracking_bimap.right.find( obj );
        if ( it == tracking_bimap.right.end() ) {
            // register new object
            LOG_DEBUG1( "Registering " << name << " with predefined id=" << id << ", ptr=0x" << obj );
            tracking_bimap.insert( bimap_value_type( id, obj ) );
            return ( std::make_pair( id, true ) );
        }
        else {
            if ( id == it->second ) {
                LOG_DEBUG2( "Already registered " << name << " id=" << id << ", ptr=0x" << obj );
                return ( std::make_pair( id, false ) );
            }
            else {
                THROW_RUNTIME_ERROR( name << " ptr=0x" << obj << " or id=" << id << " already tracked" );
            }
        }
    }

    static const object_type* fetch_object( const char* name, id_type id )
    {
        typename bimap_type::left_const_iterator it = tracking_bimap.left.find( id );
        if ( it != tracking_bimap.left.end() ) {
            LOG_DEBUG1( "Fetching " << name << " id=" << id << ", ptr=0x" << it->second );
            return ( it->second );
        }
        else {
            LOG_DEBUG2( name << " with id=" << id << " not registered" );
            return ( NULL );
        }
    }

public:
    explicit statically_tracked(const char * name, const object_type*& pT)
    : _name( name ), _pT( pT )
    {}

    const char * name() const {
        return _name;
    }
    T & value() const {
        return *_pT;
    }

    const T & const_value() const {
        return *_pT;
    }

    template<class Archivex>
    void save( Archivex & ar, const unsigned int /* file_version */ ) const
    {
        register_result_type res = register_object( name(), _pT );
        ar << boost::serialization::make_nvp( "id", res.first );
        if ( res.second ) {
            LOG_DEBUG1( name() << " id=" << res.first << ": not registered, full serialization" );
            ar << boost::serialization::make_nvp( "value", _pT );
        }
        else {
            LOG_DEBUG1( name() << " id=" << res.first << " already registered, no serialization done" );
        }
    }

    template<class Archivex>
    void load( Archivex & ar, const unsigned int /* file_version */ )
    {
        id_type id;
        ar >> boost::serialization::make_nvp( "id", id );
        _pT = fetch_object( name(), id );
        if ( _pT == NULL ) {
            LOG_DEBUG1( name() << " id=" << id << " not registered, full deserialization" );
            ar >> boost::serialization::make_nvp( "value", _pT );
            register_object_with_id( name(), _pT, id );
        }
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

/**
 *  Declare bimap as external.
 *  For particular classes, the variable should be explicitly defined.
 */
template<class T>
typename statically_tracked<T>::bimap_type statically_tracked<T>::tracking_bimap;

/**
 *  Macro to enable static object serialization tracking.
 *  Use it in source file.
 */
#define ENABLE_STATIC_TRACKING( obj_class ) \
template<> \
statically_tracked<obj_class>::bimap_type statically_tracked<obj_class>::tracking_bimap;
