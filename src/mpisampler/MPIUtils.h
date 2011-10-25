#pragma once

#include "../BasicTypedefs.h"

#include <boost/dynamic_bitset.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/timer.hpp>
#include <boost/utility/enable_if.hpp>

#include "../dynamic_bitset_utils.h"

//#define MPI_TIMEOUT_CONTROL

/**
 *  Macros for setting traits of types passed by MPI,
 */

/**
 *  Declares datatype as bitwise-serializable, 
 *  no pointer-tracked and no versioned.
 */
#define BOOST_MPI_SERIALIZATION_TRAITS( datatype ) \
    BOOST_CLASS_IMPLEMENTATION( datatype, object_serializable ) \
    BOOST_CLASS_TRACKING( datatype, track_never ) \
    BOOST_IS_BITWISE_SERIALIZABLE( datatype )

/**
 *  Defines templ_datatype<templ_param> as MPI datatype,
 *  if templ_param is MPI datatype.
 */
#define BOOST_IS_MPI_DATATYPE_TEMPLATE( templ_datatype, templ_param ) \
namespace boost {                                                     \
namespace mpi {                                                       \
template<class templ_param>                                           \
struct is_mpi_datatype< templ_datatype<templ_param> > : is_mpi_datatype<templ_param> {}; \
}}

#define BOOST_IS_BITWISE_SERIALIZABLE_TEMPLATE( templ_datatype, templ_param ) \
namespace boost {                                                     \
namespace serialization {                                             \
template<class templ_param>                                           \
struct is_bitwise_serializable< templ_datatype<templ_param> > : mpl::true_ {}; \
}}

#define BOOST_CLASS_TRACKING_TEMPLATE( templ_datatype, templ_param, level ) \
namespace boost {                                                     \
namespace serialization {                                             \
template<class templ_param>                                           \
struct tracking_level< templ_datatype<templ_param> > { \
    typedef mpl::integral_c_tag tag;         \
    typedef mpl::int_< level > type;         \
    BOOST_STATIC_CONSTANT(                   \
        int,                                 \
        value = tracking_level::type::value  \
    );                                       \
    /* tracking for a class  */              \
    BOOST_STATIC_ASSERT((                    \
        mpl::greater<                        \
            /* that is a prmitive */         \
            implementation_level< templ_datatype<templ_param> >,   \
            mpl::int_<primitive_type>        \
        >::value                             \
    ));                                      \
};                                           \
}}

#define BOOST_CLASS_IMPLEMENTATION_TEMPLATE( templ_datatype, templ_param, level ) \
namespace boost {                                                     \
namespace serialization {                                             \
template<class templ_param>                                           \
struct implementation_level_impl< templ_datatype<templ_param> > {     \
    typedef mpl::int_< level > type;                 \
    BOOST_STATIC_CONSTANT(                           \
        int,                                         \
        value = implementation_level_impl::type::value    \
    );                                               \
}; }}

#define BOOST_MPI_SERIALIZATION_TRAITS_TEMPLATE( templ_datatype, templ_param ) \
    BOOST_CLASS_IMPLEMENTATION_TEMPLATE( templ_datatype, templ_param, object_serializable ) \
    BOOST_IS_BITWISE_SERIALIZABLE_TEMPLATE( templ_datatype, templ_param ) \
    BOOST_CLASS_TRACKING_TEMPLATE( templ_datatype, templ_param, track_never )

template<class T, bool Packed>
struct MPIISendStorageTraits {
    typedef T storage_type;

    static boost::mpi::request isend( boost::mpi::communicator& comm, int dest, int tag, storage_type& storage )
    {
        return ( comm.isend( dest, tag, storage ) );
    }
};

template<class T>
struct MPIISendStorageTraits<T, true> {
    typedef boost::mpi::packed_oarchive storage_type;

    static boost::mpi::request isend( boost::mpi::communicator& comm, int dest, int tag, storage_type& storage )
    {
        return ( comm.isend( dest, tag, storage ) );
    }
};

template<class T>
struct MPIISendStorageTraits<std::vector<T>, false> {
    typedef std::vector<T> storage_type;

    static boost::mpi::request isend( boost::mpi::communicator& comm, int dest, int tag, storage_type& storage )
    {
        return ( comm.isend( dest, tag, storage.data(), storage.size() ) );
    }
};

/**
 *  Heap-allocated non-blocking send request data.
 *  It's essentially non-copyable, since MPI requests need
 *  transiently access the memory of response.
 */
template<class T, bool Packed>
class MPIISendStorageBase: public boost::noncopyable {
public:
    MPIISendStorageBase(
        boost::mpi::communicator&   comm,
        const T&                    data
    );
};

template<class T>
class MPIISendStorageBase<T, false>: public boost::noncopyable {
protected:
    typedef MPIISendStorageTraits<T, false> storage_traits;
    typedef typename storage_traits::storage_type storage_type;

    storage_type        _data;

public:
    MPIISendStorageBase(
        boost::mpi::communicator&   comm,
        const T&                    data
    ) : _data( data )
    {
    }
};

template<class T>
class MPIISendStorageBase<T, true>: public boost::noncopyable {
protected:
    typedef MPIISendStorageTraits<T, true> storage_traits;
    typedef typename storage_traits::storage_type storage_type;

    storage_type        _data;

public:
    MPIISendStorageBase(
        boost::mpi::communicator&   comm,
        const T&                    data
    ) : _data( comm )
    {
        _data << data;
    }
};

/**
 *  Heap-allocated non-blocking send request data.
 *  It's essentially non-copyable, since MPI requests need
 *  transiently access the memory of response.
 */
template<class T, bool Packed>
class MPIISendStorage: protected MPIISendStorageBase<T, Packed> {
private:
    typedef MPIISendStorageBase<T, Packed> base_type;
    typedef typename base_type::storage_traits storage_traits;

#if defined(MPI_TIMEOUT_CONTROL)
    boost::mpi::timer               _timer;
    double                          _timeout;
    int                             _tag;
    int                             _dest;
#endif
    using base_type::_data;

    bool                            _complete;
    boost::mpi::request             _isend;

public:
    MPIISendStorage(
        boost::mpi::communicator&   comm,
        int                         dest,
        int                         tag,
        const T&                    data
    #if defined(MPI_TIMEOUT_CONTROL)
        , double                    timeout = 0
    #endif
    ) : base_type( comm, data )
#if defined( MPI_TIMEOUT_CONTROL )
    , _timeout( timeout )
    , _tag( tag ), _dest( dest )
#endif
    , _complete( false )
    , _isend( storage_traits::isend( comm, dest, tag, _data ) )
    {}

    ~MPIISendStorage()
    {
        LOG_DEBUG2( "MPISendStorage::dtor()" );
        if ( !complete() ) {
            cancel();
        }
    }

#if defined( MPI_TIMEOUT_CONTROL )
    double elapsed() const {
        return ( _timer.elapsed() );
    }
#endif

    void cancel()
    {
        _isend.cancel();
        LOG_DEBUG2( "MPISendStorage::cancel() -- isend() request cancelled" );
    }

    bool complete() {
        LOG_DEBUG2( "#" << _comm.rank() << ": testing isend status" );
        if ( _complete ) return ( true );

        boost::optional<boost::mpi::status> status = _isend.test();
        if ( status ) {
            LOG_DEBUG2( "#" << _comm.rank() << ": isend to #" << status->source() << " completed" );
            _complete = true;
        }
#if defined(MPI_TIMEOUT_CONTROL)
        else if ( _timeout > 0 && elapsed() > _timeout ) {
            LOG_WARN( "#" << _comm.rank() << ": MPIISendStorage<" << _tag << "> timeout " 
                            << boost::format("%.0f") % elapsed()
                            << "s, still waiting for " << _dest << "..." );
        }
#endif
        return ( _complete );
    }
};

/**
 *  Heap-allocated non-blocking broadcast request data.
 *  It's essentially non-copyable, since MPI requests need
 *  transiently access the memory of response.
 */
template<class T, bool Packed>
class MPIIBroadcastStorage: protected MPIISendStorageBase<T, Packed> {
public:
    typedef boost::dynamic_bitset<> ranks_mask_type;

protected:
    boost::mpi::communicator&       _comm;
#if defined(MPI_TIMEOUT_CONTROL)
    boost::mpi::timer               _timer;
    double                          _timeout;
    double                          _lastReported;
#endif

private:
    typedef MPIISendStorageBase<T, Packed> base_type;
    typedef typename base_type::storage_traits storage_traits;

    using base_type::_data;

    ranks_mask_type                     _waitingUnits;
    std::vector<boost::mpi::request>    _isends;   /// MPI isend request statuses

public:
    MPIIBroadcastStorage(
        boost::mpi::communicator&   comm,
        int                         tag,
        const T&                    data,
        const ranks_mask_type&      recipients = ranks_mask_type()
#if defined(MPI_TIMEOUT_CONTROL)
        , double                      timeout = 0
#endif
    ) : base_type( comm, data ), _comm( comm )
#if defined( MPI_TIMEOUT_CONTROL )
    , _timeout( timeout )
    , _lastReported( 0 )
#endif
    , _waitingUnits( _comm.size() )
    , _isends( _comm.size() )
    {
        if ( recipients.size() > 0 ) {
            _waitingUnits = recipients;
        } else {
            _waitingUnits.set();
        }
        _waitingUnits.set( _comm.rank(), false ); // exclude itself
        foreach_bit( size_t, i, _waitingUnits ) {
            _isends[ i ] = storage_traits::isend( _comm, i, tag, _data );
        }
#if defined( MPI_TIMEOUT_CONTROL )
        _timer.restart();
#endif
    }

    ~MPIIBroadcastStorage()
    {
        LOG_DEBUG2( "MPIBroadcastStorage::dtor()" );
        if ( !complete() ) {
            cancel();
        }
    }

#if defined( MPI_TIMEOUT_CONTROL )
    double elapsed() const {
        return ( _timer.elapsed() );
    }
#endif

    void cancel()
    {
        foreach_bit( size_t, i, _waitingUnits ) {
            _isends[i].cancel();
        }
        LOG_DEBUG2( "MPIBroadcastStorage::cancel() -- all isend() requests cancelled" );
    }

    const ranks_mask_type& waitingUnits() const {
        return ( _waitingUnits );
    }

    bool complete() {
        LOG_DEBUG2( "#" << _comm.rank() << ": testing ibroadcast status" );
        if ( _waitingUnits.none() ) return ( true );

        foreach_bit( size_t, i, _waitingUnits ) {
            boost::optional<boost::mpi::status> status = _isends[ i ].test();
            if ( status ) {
                LOG_DEBUG2( "#" << _comm.rank() << ": isend to #" << status->source() << " completed" );
                _waitingUnits.set( i, false );
            }
        }
        if ( _waitingUnits.none() ) {
            LOG_DEBUG1( "#" << _comm.rank() << ": ibroadcast complete" );
            return ( true );
        } else {
#if defined(MPI_TIMEOUT_CONTROL)
            if ( _timeout > 0 && elapsed() > _lastReported + _timeout  ) {
                std::ostringstream strstr;
                foreach_bit( size_t, i, _waitingUnits ) {
                    strstr << ' ' << i;
                }
                LOG_WARN( "#" << _comm.rank() << ": MPIBroadcastStorage timeout " 
                              << boost::format("%.0f") % elapsed() << "s, still waiting for " << strstr.str() << "..." );
                _lastReported = elapsed();
            }
#endif
            return ( false );
        }
    }
};

boost::mpi::graph_communicator create_star_communicator( boost::mpi::communicator& parent, int center );
