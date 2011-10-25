#pragma once

#include "BasicTypedefs.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/enable_shared_from_this.hpp>

/**
    Wraps an object and clones it if modification is required.
    To clone object, the descendant class should call unshare() method.
 */
template<class Object>
class LazyClone {
protected:
    typedef Object object_type;
    typedef boost::shared_ptr<object_type> shared_object_ptr_type;
    typedef object_type const* const_object_ptr_type;

    bool _isRef;
    const_object_ptr_type pExternal;
    shared_object_ptr_type pShared;

    virtual void wrappedReplaced( shared_object_ptr_type& pShared )
    {
    };

    void unshare()
    {
        if ( _isRef ) {
            BOOST_ASSERT( pShared || pExternal );
            if ( pShared ) {
                if ( !pShared.unique() ) { // skip duplication if object is already unique
                    pShared = boost::make_shared<object_type>( *pShared );
                    LOG_DEBUG2( "Lazy clone unshared: new ext=" << (long int)( pExternal ) << " shared=" << (long int)(pShared.get()) );
                }
            }
            else if ( pExternal ) {
                pShared = boost::make_shared<object_type>( *pExternal );
                pExternal = NULL;
                LOG_DEBUG2( "Lazy clone unshared from external: new ext=" << (long int)( pExternal ) << " shared=" << (long int)(pShared.get()) );
            }
            else {
                throw std::runtime_error( "No object referenced" );
            }
            wrappedReplaced( pShared );
            _isRef = false;
        }
    }

    const object_type& wrapped() const {
        BOOST_ASSERT( pExternal || pShared );
        return ( pExternal ? *pExternal : *pShared );
    }

    /**
        @warning
        Wrapped object is implicitly cloned upon this call.
     */
    object_type& wrapped() {
        unshare();
        BOOST_ASSERT( pShared );
        return ( *pShared );
    }

    class LockedPointer {
    private:
        const object_type* pExternal;
        shared_object_ptr_type pShared;
    public:
        LockedPointer( const object_type* pExternal, const shared_object_ptr_type& pShared )
        : pExternal( pExternal ), pShared( pShared )
        {
            BOOST_ASSERT( ( pShared && !pExternal ) || ( !pShared && pExternal ) );
        }

        operator const const_object_ptr_type() const {
            return ( pExternal != NULL ? pExternal : pShared.get() );
        }
        const object_type& operator*() const {
            return ( pExternal != NULL ? *pExternal : *pShared.get() );
        }
        const const_object_ptr_type operator->() const {
            return ( pExternal != NULL ? pExternal : pShared.get() );
        }

        bool operator==( const LockedPointer& that ) const {
            return ( (const_object_ptr_type)(*this) == (const_object_ptr_type)that );
        }
    };

    typedef LockedPointer locked_ptr_type;

    /**
        Locks object and returns pointer to prevent it's destruction
        or modification of the pointer
     */
    const LockedPointer lock() const {
        return ( LockedPointer( pExternal, pShared ) );
    }

public:
    LazyClone( const shared_object_ptr_type& pObj, bool isRef )
    : _isRef( isRef )
    , pExternal( NULL )
    , pShared( pObj )
    {
        BOOST_ASSERT( pObj );
        LOG_DEBUG2( "Lazy clone ref=" << _isRef << " ext=" << (long int)( pExternal ) << " shared=" << (long int)(pShared.get()) );
    }

    LazyClone( const object_type& obj )
    : _isRef( true )
    , pExternal( &obj )
    {
        LOG_DEBUG2( "Lazy clone ref=" << _isRef << " ext=" << (long int)( pExternal ) << " shared=" << (long int)(pShared.get()) );
    }

    LazyClone( const LazyClone<object_type>& that )
    : _isRef( true ) // note: in 'that' it might be not shared
    , pExternal( that.pExternal )
    , pShared( that.pShared )
    {
        LOG_DEBUG2( "Lazy clone ref=" << _isRef << " ext=" << (long int)( pExternal ) << " shared=" << (long int)(pShared.get()) );
    }

    virtual ~LazyClone()
    {}

    /**
        Checks if instance references object created somewhere else.
        @remark
        Also returns true, if this instance holds the unique pointer
        to object created somewhere else.
     */
    bool isRef() const {
        // tip: retu
        return ( _isRef && ( !pShared || !pShared.unique() ) );
    }

    operator const object_type&() const
    {
        return ( wrapped() );
    }

    const shared_object_ptr_type& shared_ptr() const {
        return ( pShared );
    }

    shared_object_ptr_type& shared_ptr() {
        return ( pShared );
    }
};
