#pragma once

#include "../BasicTypedefs.h"

#include <numeric>
#include <map>
#include <boost/unordered_set.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/timer.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "MPIUtils.h"
#include "../eesampler/TurbineCascadeUnit.h"

#define MPI_LANDSCAPE_BROADCAST_TIMEOUT 10
#define MPI_LANDSCAPE_RESPONSE_TIMEOUT 10

namespace mpi = boost::mpi;

enum MPISamplerTag {
    MTAG_TurbineStatuses = 1,

    MTAG_EEJumpRequest = 2,
    MTAG_EEJumpResponse = 3,

    MTAG_EnergyLandscapeRequest = 4,
    MTAG_EnergyLandscapeResponse = 5,

    MTAG_Sample = 6,
    MTAG_AdjustEnergyLadder = 7,

    MTAG_GeneratedParticles = 8,

    MTAG_Shutdown = 9
};

/**
 *  Status of particle turbine.
 */
struct TurbineStatus {
    bool            cascadeBurnIn;  /** cascade burning-in flag */ 
    size_t          iteration;      /** # of iterations done */
    size_t          prisonersCount; /** # of prisoners in the jail */

    TurbineStatus()
    : cascadeBurnIn( true )
    , iteration( 0 ), prisonersCount( 0 )
    {}
    template<class DynamicParticleFactory>
    TurbineStatus( const ParticleTurbineBase<DynamicParticleFactory>& turbine, bool cascadeBurnIn )
    : cascadeBurnIn( cascadeBurnIn )
    , iteration( turbine.iteration() ), prisonersCount( turbine.prisonersCount() )
    {}

    template<class Archive>
    void serialize( Archive & ar, const unsigned int version )
    {
        ar & BOOST_SERIALIZATION_NVP( cascadeBurnIn );
        ar & BOOST_SERIALIZATION_NVP( iteration );
        ar & BOOST_SERIALIZATION_NVP( prisonersCount );
    }
};

BOOST_MPI_SERIALIZATION_TRAITS( TurbineStatus )
BOOST_IS_MPI_DATATYPE( TurbineStatus )

/**
 *  Message of turbine status.
 */
struct TurbineStatusMessage {
    turbine_ix_t    turbineIx;
    TurbineStatus   status;

    TurbineStatusMessage()
    : turbineIx( TURBINE_NA )
    {}
    TurbineStatusMessage( turbine_ix_t turbineIx, const TurbineStatus& status )
    : turbineIx( turbineIx ), status( status )
    {}

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( turbineIx );
        ar & BOOST_SERIALIZATION_NVP( status );
    }
};

BOOST_MPI_SERIALIZATION_TRAITS( TurbineStatusMessage )
BOOST_IS_MPI_DATATYPE( TurbineStatusMessage )

/**
 *  Body of equi-energy jump requesting message (MTAG_EEJumpRequest).
 *  Sent to hot turbine by connected cold turbine.
 */
template<typename ParticleEnergyEval>
struct EEJumpRequestMessage {
    turbine_ix_t    acceptorIx;     /** turbine requesting EE jump */
    turbine_ix_t    donorIx;        /** turbine to process EE jump and respond */
    ParticleEnergyEval      energyEval;     /** evaluator of particle energies */
    EEJumpParams    params;         /** EE jump parameters */

    EEJumpRequestMessage()
    : acceptorIx( TURBINE_NA ), donorIx( TURBINE_NA )
    {};

    EEJumpRequestMessage(
        turbine_ix_t        acceptorIx,
        turbine_ix_t        donorIx,
        const ParticleEnergyEval&   energyEval,
        const EEJumpParams& params
    ) : acceptorIx( acceptorIx ), donorIx( donorIx )
      , energyEval( energyEval ), params( params )
    {}

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( acceptorIx );
        ar & BOOST_SERIALIZATION_NVP( donorIx );
        ar & BOOST_SERIALIZATION_NVP( energyEval );
        ar & BOOST_SERIALIZATION_NVP( params );
    }
};

BOOST_MPI_SERIALIZATION_TRAITS( EEJumpParams )
BOOST_MPI_SERIALIZATION_TRAITS_TEMPLATE( EEJumpRequestMessage, ParticleEnergyEval )
BOOST_IS_MPI_DATATYPE_TEMPLATE( EEJumpRequestMessage, ParticleEnergyEval )

/**
  *  Body of equi-energy jump response (MTAG_EEJumpResponse).
  *  Sent in response to MTAG_EEJumpRequest by hot turbine to cold one.
 **/
template<typename ParticleEnergyEval>
struct EEJumpResponseMessage: public EEJumpRequestMessage<ParticleEnergyEval> {
    typedef ParticleEnergyEval particle_energy_eval_type;
    typedef typename particle_energy_eval_type::particle_type particle_type;
    typedef CascadeParticle<particle_type> cascade_particle_type;
    typedef EEJumpRequestMessage<ParticleEnergyEval> request_message_type;

    boost::optional<cascade_particle_type> result; /** EE jump result, not defined,
                                                        if EE jump not available for given turbine */

    using request_message_type::donorIx;

    EEJumpResponseMessage() 
    : request_message_type()
    {}

    EEJumpResponseMessage( const request_message_type& request ) 
    : request_message_type( request )
    {}

    operator bool() const {
        return ( donorIx != TURBINE_NA );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::make_nvp( "request", boost::serialization::base_object<request_message_type>( *this ) );
        ar & BOOST_SERIALIZATION_NVP( result );
    }
};

//BOOST_MPI_SERIALIZATION_TRAITS_TEMPLATE( EEJumpResponseMessage, Particle )
// NO BOOST_IS_MPI_DATATYPE_TEMPLATE( EEJumpResponseMessage, Particle ) // NO MPI datatype, because of optional

/**
 *  Message with particles array.
 */
template<class Particle>
struct ParticleMessage {
    typedef Particle particle_type;

    turbine_ix_t    sourceIx;
    turbine_ix_t    destinationIx;
    particle_type   particle;

    ParticleMessage()
    : sourceIx( TURBINE_NA ), destinationIx( TURBINE_NA )
    {}

    ParticleMessage( turbine_ix_t sourceIx, turbine_ix_t destinationIx, const particle_type& particle )
    : sourceIx( sourceIx ), destinationIx( destinationIx ), particle( particle )
    {}

    operator bool() const {
        return ( sourceIx != TURBINE_NA );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( sourceIx );
        LOG_DEBUG1( "Serializing message source " << sourceIx );
        ar & BOOST_SERIALIZATION_NVP( destinationIx );
        ar & BOOST_SERIALIZATION_NVP( particle );
    }
};

//BOOST_MPI_SERIALIZATION_TRAITS_TEMPLATE( ParticleMessage, Particle )
//BOOST_IS_MPI_DATATYPE_TEMPLATE( ParticleMessage, Particle )

template<class Particle>
typename boost::enable_if< boost::mpi::is_mpi_datatype<ParticleMessage<Particle> >, std::vector< ParticleMessage<Particle> > >::type
receiveParticles(
    boost::mpi::communicator&   comm,
    const boost::mpi::status&   msgStatus
){
    typedef ParticleMessage<Particle> particle_message_type;
    std::vector<particle_message_type> particles( msgStatus.count<particle_message_type>().get_value_or(0) );
    LOG_DEBUG1( "#" << comm.rank() << ": got " << particles.size() << " particle(s)" );
    comm.recv( msgStatus.source(), msgStatus.tag(), particles.data(), particles.size() );
#if defined(_DEBUG)
    for ( size_t i = 0; i < particles.size(); i++ ) {
        const particle_message_type& particle = particles[ i ];
        BOOST_ASSERT( particle );
    }
#endif
    return ( particles );
}

template<class Particle>
typename boost::disable_if< boost::mpi::is_mpi_datatype<ParticleMessage<Particle> >, std::vector< ParticleMessage<Particle> > >::type
receiveParticles(
    boost::mpi::communicator&   comm,
    const boost::mpi::status&   msgStatus
){
    typedef ParticleMessage<Particle> particle_message_type;
    mpi::packed_iarchive iar( comm );
    std::vector<particle_message_type>  res;
    comm.recv( msgStatus.source(), msgStatus.tag(), iar );
    iar >> res;
    LOG_DEBUG1( "#" << comm.rank() << ": got " << res.size() << " particle(s)" );
    return ( res );
}

/**
    *  Updated EE energy ladder.
    *  Sent after any turbine has finished burning in to all units.
    */
template<class Particle>
struct EnergyLadderMessage {
    typedef Particle particle_type;

    turbine_ix_t            originIx;
    TurbineEnergyLadder     ladder;

    particle_type           iniParticle;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( originIx );
        ar & BOOST_SERIALIZATION_NVP( ladder );
        ar & BOOST_SERIALIZATION_NVP( iniParticle );
    }
};

BOOST_MPI_SERIALIZATION_TRAITS( EnergyTransform )
// BOOST_MPI_SERIALIZATION_TRAITS( TurbineEnergyLadder )  not MPI, since ladder is of variable size
// BOOST_MPI_SERIALIZATION_TRAITS_TEMPLATE( EnergyLadderMessage, Particle )  not MPI, since ladder is of variable size

template<class Particle>
class MPITurbineConnection;

/**
 *  Implementation of IUnitCommunicator interface over MPI.
 */
template<typename ParticleEnergyEval>
class MPIUnitCommunicator {
public:
    typedef ParticleEnergyEval energy_eval_type;
    typedef typename energy_eval_type::particle_type particle_type;
    typedef MPIUnitCommunicator<ParticleEnergyEval> communicator_type;
    typedef MPITurbineConnection<ParticleEnergyEval> turbine_connection_type;
    //typedef TurbineCascadeUnit<ParticleCollector, ParticleEnergyEval, DynamicParticleFactory, ParticleGenerator, communicator_type> cascade_unit_type;
    typedef boost::mpi::communicator mpi_communicator_type;
    typedef TurbineEnergyLadder::energy_landscape_type energy_landscape_type;
    typedef boost::mpi::timer timer_type;

private:
    /**
     *  Heap-allocated energy landscape request data.
     *  It's essentially non-copyable, since mpi requests need
     *  transiently access the memory of response.
     */
    class MPIEnergyLandscapeRequestStorage: public MPIIBroadcastStorage<turbine_ix_t, false> {
    public:
        typedef TurbineEnergyLadder::energy_landscape_type energy_landscape_type;

    private:
        typedef MPIIBroadcastStorage<turbine_ix_t, false> super;

        super::ranks_mask_type              _waitingUnits;
        energy_landscape_type               _cumulativeLandscape;
        std::vector<energy_landscape_type>  _landscapes;
        std::vector<mpi::request>           _landscapeResponseRequests;   /// MPI landscape response statuses
#if defined(MPI_TIMEOUT_CONTROL)
        double                              _lastReported;
#endif

    public:
        MPIEnergyLandscapeRequestStorage( communicator_type& comm, turbine_ix_t senderIx
#if defined(MPI_TIMEOUT_CONTROL)
            , double            timeout = 0
#endif
        );
        ~MPIEnergyLandscapeRequestStorage();
        bool complete();
        const energy_landscape_type& energyLandscape() const {
            return ( _cumulativeLandscape );
        }
    };

public:
    /**
     *  Stack-allocated wrapper for MPIEnergyLandscapeRequestStorage.
     */
    class MPIEnergyLandscapeRequestStatus {
    private:
        boost::shared_ptr<MPIEnergyLandscapeRequestStorage> _p;

    public:
        typedef TurbineEnergyLadder::energy_landscape_type energy_landscape_type;
        typedef typename MPIEnergyLandscapeRequestStorage::ranks_mask_type ranks_mask_type;

        MPIEnergyLandscapeRequestStatus( communicator_type& comm, turbine_ix_t senderIx )
        : _p( new MPIEnergyLandscapeRequestStorage( comm, senderIx
#if defined(MPI_TIMEOUT_CONTROL)
                    , MPI_LANDSCAPE_BROADCAST_TIMEOUT
#endif
        ) ) {};
        bool complete() { return ( _p->complete() ); }
        const energy_landscape_type& energyLandscape() const { return ( _p->energyLandscape() ); }
    };

    typedef MPIEnergyLandscapeRequestStatus landscape_request_status_type;

private:
    mpi_communicator_type&              _comm;
    int                                 _sampleCollectorRank;

    typedef typename boost::mpi::is_mpi_datatype<particle_type>::type is_particle_mpi_datatype;
    typedef boost::shared_ptr<turbine_connection_type> turbine_connection_pointer_type;
    typedef boost::unordered_map<turbine_ix_t, turbine_connection_pointer_type> turbine_connection_map_type;
    typedef boost::unordered_map<turbine_ix_t, TurbineStatus> turbine_status_map_type;
    typedef MPIIBroadcastStorage< std::vector<TurbineStatusMessage>, false > turbine_status_ibroadcast_type;
    typedef MPIISendStorage< std::vector<TurbineStatusMessage>, false > turbine_status_isend_type;
    typedef MPIIBroadcastStorage< EnergyLadderMessage<particle_type>, true > energy_ladder_status_ibroadcast_type;
    typedef MPIISendStorage< EEJumpResponseMessage<energy_eval_type>, true > jump_response_isend_request_type;
    typedef MPIISendStorage< std::vector< ParticleMessage<particle_type> >, !is_particle_mpi_datatype::value > particle_isend_request_type;
    typedef MPIISendStorage<energy_landscape_type, true> landscape_response_isend_request_type;
    typedef boost::ptr_vector<particle_isend_request_type> particle_isend_request_container_type;
    typedef boost::ptr_vector<jump_response_isend_request_type> jump_response_isend_request_container_type;
    typedef boost::ptr_vector<landscape_response_isend_request_type> landscape_response_isend_request_container_type;

    particle_isend_request_container_type           _particleISendRequests;
    jump_response_isend_request_container_type      _jumpResponseISendRequests;
    landscape_response_isend_request_container_type _landscapeResponseISendRequests;
    boost::scoped_ptr<turbine_status_ibroadcast_type>   _unitStatusesIBroadcast;
    boost::scoped_ptr<turbine_status_isend_type>        _unitStatusISend;
    boost::scoped_ptr<energy_ladder_status_ibroadcast_type> _energyLadderIBroadcastRequest;

    turbine_connection_map_type         _turbineConnections;
    turbine_status_map_type             _externalTurbineStatuses;

protected:
    template<class CascadeUnit>
    bool processIncomingMessage( CascadeUnit& unit, const mpi::status& msgStatus );
    template<class CascadeUnit>
    void processJumpRequest( const CascadeUnit& unit, const mpi::status& msgStatus );
    template<class CascadeUnit>
    void updateExternalStatuses( CascadeUnit& unit, const mpi::status& msgStatus );

    typedef ParticleMessage<particle_type> particle_message_type;

public:
    MPIUnitCommunicator( mpi_communicator_type& comm, int sampleCollectorRank )
    : _comm( comm ), _sampleCollectorRank( sampleCollectorRank )
    {
    }

    mpi_communicator_type& communicator() {
        return ( _comm );
    }
    int rank() const {
        return ( _comm.rank() );
    }
    boost::optional<TurbineStatus> lastTurbineStatus( turbine_ix_t turbineIx ) const;

    template<class CascadeUnit>
    void sendUnitStatus( const CascadeUnit& unit );
    template<class CascadeUnit>
    void broadcastUnitStatuses( const CascadeUnit& unit );
    void broadcastShutdown( turbine_ix_t originIx );
    void broadcastEnergyLadder( turbine_ix_t originIx, const TurbineEnergyLadder& ladder, const particle_type& iniParticle );
    template<class CascadeUnit>
    bool processExternalEvents( CascadeUnit& unit, bool wait );
    template<class CascadeUnit>
    landscape_request_status_type requestEnergyLandscape( const CascadeUnit& unit, turbine_ix_t senderIx ) {
            return ( landscape_request_status_type( *this, senderIx ) );
    }
    turbine_connection_type* createTurbineConnection( const TurbineNode& externalTurbine ) {
        return ( new turbine_connection_type( *this, externalTurbine ) );
    }
    void sendSample( turbine_ix_t originIx, const particle_type& sample );
    template<class ParticlesIterator>
    void sendParticles( int destinationRank, int tag, turbine_ix_t destinationTurbineIx, turbine_ix_t originIx, 
                        const ParticlesIterator& particlesBegin, const ParticlesIterator& particlesEnd );
};

/**
 *  Implementation of ITurbineConnection interface over MPI.
 */
template<typename ParticleEnergyEval>
class MPITurbineConnection {
public:
    typedef MPIUnitCommunicator<ParticleEnergyEval> communicator_type;
    typedef MPITurbineConnection<ParticleEnergyEval> connection_type;

    //typedef typename communicator_type::cascade_unit_type cascade_unit_type;
    typedef ParticleEnergyEval particle_energy_eval_type;
    typedef typename particle_energy_eval_type::particle_type particle_type;
    typedef typename communicator_type::mpi_communicator_type mpi_communicator_type;
    typedef CascadeParticle<particle_type> cascade_particle_type;
    typedef EEJumpResponseMessage<particle_energy_eval_type> jump_response_message_type;

    friend class MPIUnitCommunicator<ParticleEnergyEval>;

    class MPIJumpRequestStatus;
    friend class MPIJumpRequestStatus;

private:
    /**
     *  Equi-energy jump request status.
     */
    class MPIJumpRequestStorage: boost::noncopyable {
    private:
        typedef EEJumpRequestMessage<ParticleEnergyEval> request_message_type;
        request_message_type        _request;
        jump_response_message_type  _response;
        mpi::request                _jumpResponseRequest;   /// MPI jump response status

    public:
        MPIJumpRequestStorage( connection_type& conn, turbine_ix_t initiatorIx,
                               const ParticleEnergyEval& energyEval, const EEJumpParams& params )
        : _request( initiatorIx, conn.externalTurbineIx(), energyEval, params )
        , _jumpResponseRequest( mpi::request( conn._comm.communicator().irecv( conn.externalUnitRank(), MTAG_EEJumpResponse, _response ) ) )
        {
            LOG_DEBUG1( "MPIJumpRequestStorage::ctor() from #" << initiatorIx << " to " << conn.externalTurbineIx() );
            conn._comm.communicator().isend( conn.externalUnitRank(), MTAG_EEJumpRequest, _request );
        }
        bool complete() {
            LOG_DEBUG2( "EE jump response testing" );
            boost::optional<mpi::status> status = _jumpResponseRequest.test();
            LOG_DEBUG1_IF( status, "Got EE jump response" );
            BOOST_ASSERT( !status || _response ); // if we had answer, response is set
            return ( status );
        }
        const boost::optional<cascade_particle_type> particle() const {
            return ( _response.result );
        }
        ~MPIJumpRequestStorage() {
            LOG_DEBUG1( "MPIJumpRequestStorage::dtor()" );
            if ( !_response && !_jumpResponseRequest.test() ) {
                LOG_DEBUG1( "MPIJumpRequestStorage::dtor() cancelling request" );
                _jumpResponseRequest.cancel();
            }
        }
    };

public:
    class MPIJumpRequestStatus {
    private:
        boost::shared_ptr<MPIJumpRequestStorage> _p;
    public:
        MPIJumpRequestStatus( connection_type& conn, turbine_ix_t initiatorIx,
                              const particle_energy_eval_type& energyEval,
                              const EEJumpParams& params )
        : _p( new MPIJumpRequestStorage( conn, initiatorIx, energyEval, params ) )
        {}
        bool complete() { return ( _p->complete() ); }
        const boost::optional<cascade_particle_type> particle() const { return ( _p->particle() ); }
    };
    typedef MPIJumpRequestStatus jump_request_status_type;

private:
    communicator_type&          _comm;
    TurbineNode                 _externalTurbineInfo;

public:
    MPITurbineConnection(
        communicator_type&      comm,
        const TurbineNode&      externalTurbine
    ) : _comm( comm )
    , _externalTurbineInfo( externalTurbine )
    {
    }

    turbine_ix_t externalTurbineIx() const { return ( _externalTurbineInfo.turbineIx ); }
    int externalUnitRank() const { return ( _externalTurbineInfo.unitIx ); }
    size_t level() const { return ( _externalTurbineInfo.level ); }

    boost::optional<TurbineStatus> externalTurbineStatus() const {
        return ( _comm.lastTurbineStatus( externalTurbineIx() ) );
    }

    size_t prisonersCount() const {
        return ( externalTurbineStatus() ? externalTurbineStatus()->prisonersCount : 0 );
    }

    jump_request_status_type requestJump( turbine_ix_t initiatorIx,
                                          const particle_energy_eval_type& energyEval,
                                          const EEJumpParams& params
    ){
        return ( jump_request_status_type( *this, initiatorIx, energyEval, params ) );
    }

    template<class ParticleIterator>
    bool sendParticles( turbine_ix_t originIx, const ParticleIterator& particlesBegin, const ParticleIterator& particlesEnd )
    {
        _comm.sendParticles( externalUnitRank(), MTAG_GeneratedParticles, externalTurbineIx(), originIx,
                             particlesBegin, particlesEnd );
        return ( false );
    }

    template<class CascadeUnit>
    bool processExternalEvents( CascadeUnit& unit, bool wait ) {
        return ( _comm.processExternalEvents( unit, wait ) );
    }
};

/**
 *  Utility function to get rank of samples collector.
 */
template<class ParticleCollector>
int MPIGetCollectorRank(
    boost::mpi::communicator& comm,
    const TurbineCascadeParams& params,
    ParticleCollector* pCollector
){
    // find who is collecting and broadcast
    // this might be stupid, but otherwise checks that MPI is ok
    int collectorRank = -1;
    std::vector<int>   collectorMask( comm.size() );
    boost::mpi::all_gather( comm, (int)(pCollector != NULL), collectorMask );
    std::vector<int>::const_iterator collectorIt = std::find( collectorMask.begin(), collectorMask.end(), true );
    if ( collectorIt == collectorMask.end() ) {
        return ( -1 ); // not found
    }
    else {
        collectorRank = collectorIt - collectorMask.begin();
    }
    LOG_DEBUG1( "Unit #" << collectorRank << " is a collector" );

    return ( collectorRank );
}

template<class Container>
void erase_completed( Container& container )
{
   for ( size_t i = 0; i < container.size(); ) {
        if ( container[ i ].complete() ) {
            container.erase( container.begin() + i );
        }
        else {
            i++;
        }
    }
}

/******************************************************************************
 * MPIUnitCommunicator implementation
 ******************************************************************************/

template<class Particle>
template<class CascadeUnit>
bool MPIUnitCommunicator<Particle>::processExternalEvents( 
    CascadeUnit&    unit,
    bool            wait
){
    // housekeeping
    erase_completed( _particleISendRequests );
    erase_completed( _jumpResponseISendRequests );
    erase_completed( _landscapeResponseISendRequests );
    if ( _energyLadderIBroadcastRequest && _energyLadderIBroadcastRequest->complete() ) _energyLadderIBroadcastRequest.reset();
    if ( _unitStatusISend && _unitStatusISend->complete() ) _unitStatusISend.reset();
    if ( _unitStatusesIBroadcast && _unitStatusesIBroadcast->complete() ) _unitStatusesIBroadcast.reset();

    // probe for new events
    if ( wait ) {
        LOG_DEBUG1( "#" << rank() << ": probe()..." );
        mpi::status status = _comm.probe( mpi::any_source, mpi::any_tag );
        return ( processIncomingMessage( unit, status ) );
    }
    else {
        LOG_DEBUG2( "#" << rank() << ": iprobe()..." );
        boost::optional<mpi::status> oStatus = _comm.iprobe( mpi::any_source, mpi::any_tag );
        if ( oStatus ) {
            return ( processIncomingMessage( unit, *oStatus ) );
        }
    }
    return ( false );
}

template<class Particle>
template<class CascadeUnit>
bool MPIUnitCommunicator<Particle>::processIncomingMessage(
    CascadeUnit&                unit,
    const boost::mpi::status&   msgStatus
){
    typedef CascadeUnit cascade_unit_type;
    bool res = false;

    LOG_DEBUG2_IF( msgStatus.tag() != MTAG_UnitStatus, "#" << rank() << ": got tag " << msgStatus.tag() << " from #" << msgStatus.source() );
    switch ( msgStatus.tag() ) {
        case MTAG_Sample:
        {
            LOG_DEBUG2( "Getting sample" );
            std::vector< ParticleMessage<particle_type> > samples = receiveParticles<particle_type>( _comm, msgStatus );
            for ( size_t i = 0; i < samples.size(); i++ ) {
                bool res = unit.sampleParticle( samples[i].sourceIx, samples[i].particle );
                if ( res ) break;
            }
        }
        break;

        case MTAG_EEJumpRequest:
            LOG_DEBUG2( "#" << rank() << ": got EE jump request" );
            processJumpRequest( unit, msgStatus );
            break;

        case MTAG_EEJumpResponse:
        {
            LOG_DEBUG2( "EEJump response skipped: processed by MPIEEJumpRequestStatus" );
#if 0
            _jumpStatuses.find( response.id );
            if ( _turbines.find( response.acceptorIx ) == _turbines.end() ) {
                throw std::runtime_error( "Acceptor turbine " << jump.turbineIx << " not found in sampling unit " << _comm.rank() );
            }
                ConnectionTurbine& turbine = _turbines.find( jump.turbineIx );
                if ( turbine->hotConnection ) {
                    turbine->hotConnection->processJumpResponse( jump );
                }
                else {
                    LOG_DEBUG0( "Receive jump response message for turbine " << jump.turbineIx 
                                << ", which is not connected to any hot turbine" );
                }
            }
            else {
            }
#endif
        }
        break;

        case MTAG_GeneratedParticles:
        {
            std::vector<particle_message_type> particles = receiveParticles<particle_type>( _comm, msgStatus );
            LOG_DEBUG1( "#" << rank() << ": got " << particles.size() << " generated particle(s)" );
            for ( size_t i = 0; i < particles.size(); i++ ) {
                const ParticleMessage<particle_type>& particle = particles[ i ];
                typename cascade_unit_type::turbine_iterator tIt = unit.turbines().find( particle.destinationIx );
                if ( tIt == unit.turbines().end() ) {
                    THROW_RUNTIME_ERROR( "Cannot find turbine #" << particle.destinationIx << " for incoming generated particles in unit #" << unit.index() );
                }
                tIt->second->detainParticle( unit.elapsed(), particle.particle );
            }
        }
        break;

        case MTAG_TurbineStatuses:
            updateExternalStatuses( unit, msgStatus );
            break;

        case MTAG_EnergyLandscapeRequest:
        {
            LOG_DEBUG2( "#" << rank() << ": got energy landscape request from unit #" << msgStatus.source() );
            turbine_ix_t originIx;
            _comm.recv( msgStatus.source(), msgStatus.tag(), originIx );
            typename cascade_unit_type::energy_landscape_type localLandscape = unit.localEnergyLandscape();
            LOG_DEBUG2( "#" << rank() << ": sending " << localLandscape.size() << " particle energies to turbine #" << originIx );
            _landscapeResponseISendRequests.push_back( new landscape_response_isend_request_type( _comm, msgStatus.source(), 
                                                                                                  MTAG_EnergyLandscapeResponse, localLandscape
#if defined(MPI_TIMEOUT_CONTROL)
                                                                                                    , MPI_LANDSCAPE_RESPONSE_TIMEOUT
#endif
                                                                                                ) );
        }
        break;

        case MTAG_EnergyLandscapeResponse:
            LOG_DEBUG2( "Energy Landscape response from #" << msgStatus.source() << " skipped: processed by MPIEnergyLandscapeRequestStatus" );
            break;

        case MTAG_AdjustEnergyLadder:
        {
            LOG_DEBUG1( "#" << rank() << ": got energy ladder adjustment" );
            EnergyLadderMessage<particle_type> msg;
            _comm.recv( msgStatus.source(), msgStatus.tag(), msg );
            unit.adjustUnitEnergyLadder( msg.originIx, msg.ladder, msg.iniParticle );
        }
        break;

        case MTAG_Shutdown:
            LOG_DEBUG1( "#" << rank() << ": got shutdown signal" );
            res = true;
            break;

        default:
            THROW_RUNTIME_ERROR( "Unrecognized message tag " << msgStatus.tag() << " from rank #" << msgStatus.source() );
            break;
    }
    return ( res );
}

template<typename ParticleEnergyEval>
template<class CascadeUnit>
void MPIUnitCommunicator<ParticleEnergyEval>::processJumpRequest(
    const CascadeUnit&    unit,
    const mpi::status&    msgStatus
){
    typedef CascadeUnit cascade_unit_type;
    EEJumpRequestMessage<energy_eval_type>  request;
    _comm.recv( msgStatus.source(), MTAG_EEJumpRequest, request );
    EEJumpResponseMessage<energy_eval_type> response( request );
    typename cascade_unit_type::const_turbine_iterator donorIt = unit.turbines().find( request.donorIx );
    // process jump request in appropriate turbines of the unit
    if ( donorIt == unit.turbines().end() ) {
        LOG_DEBUG0( "Donor turbine " << request.donorIx << " not found in unit #" << unit.index() );
        //throw std::runtime_error( "Connected turbine is not created yet in given unit" );
    }
    else {
        response.result = donorIt->second->processJumpRequest( request.params, request.energyEval, false );
    }
    _jumpResponseISendRequests.push_back( new jump_response_isend_request_type( _comm, msgStatus.source(), MTAG_EEJumpResponse, response ) );
}

/**
 *  Updates local cache of external turbines statuses.
 */
template<typename ParticleEnergyEval>
template<class CascadeUnit>
void MPIUnitCommunicator<ParticleEnergyEval>::updateExternalStatuses(
    CascadeUnit&        unit,
    const mpi::status&  msgStatus
){
    typedef CascadeUnit cascade_unit_type;
    std::vector<TurbineStatusMessage> msgs( msgStatus.count<TurbineStatusMessage>().get_value_or(0) );
    LOG_DEBUG2( "#" << _comm.rank() << ": receiving unit #" << msgStatus.source() 
                    << " " << msgs.size() << " status(es)" );
    _comm.recv( msgStatus.source(), msgStatus.tag(), msgs.data(), msgs.size() );
    LOG_DEBUG2( "#" << _comm.rank() << ": received unit #" << msgStatus.source() <<
                " " << msgs.size() << " status(es), last iteration #" << msgs.back().status.iteration );
    // use only last status
    //const UnitStatusMessage& lastMsg = msgs.back();
    for ( size_t msgIx = 0; msgIx < msgs.size(); msgIx++ ){
        const TurbineStatusMessage& msg = msgs[ msgIx ];
        if ( unit.cascadeBurnIn() && !msg.status.cascadeBurnIn ) {
            unit.processBurnInOver();
        }
        turbine_status_map_type::iterator curIt = _externalTurbineStatuses.find( msg.turbineIx );
        if ( curIt == _externalTurbineStatuses.end() ) {
            unit.processTurbineCreated( msg.turbineIx );
        }
        _externalTurbineStatuses[ msg.turbineIx ] = msg.status;
    }
}

template<typename ParticleEnergyEval>
void MPIUnitCommunicator<ParticleEnergyEval>::broadcastEnergyLadder(
    turbine_ix_t                originIx,
    const TurbineEnergyLadder&  ladder,
    const particle_type&        iniParticle
){
    EnergyLadderMessage<particle_type> msg;
    msg.originIx = originIx;
    msg.ladder = ladder;
    msg.iniParticle = iniParticle;

    if ( _energyLadderIBroadcastRequest ) _energyLadderIBroadcastRequest.reset();
    _energyLadderIBroadcastRequest.reset( new energy_ladder_status_ibroadcast_type( _comm, MTAG_AdjustEnergyLadder, msg ) );
}

template<typename ParticleEnergyEval>
template<class CascadeUnit>
void MPIUnitCommunicator<ParticleEnergyEval>::sendUnitStatus(
    const CascadeUnit&  unit
){
    if ( _comm.rank() == _sampleCollectorRank ) {
        THROW_RUNTIME_ERROR( "Turbine Unit cannot send status to itself" );
    }
    typedef CascadeUnit cascade_unit_type;
    std::vector<TurbineStatusMessage> msgs( unit.turbines().size() );
    int i = 0;
    for ( typename cascade_unit_type::const_turbine_iterator tit = unit.turbines().begin(); tit != unit.turbines().end(); ++tit ) {
        msgs[ i++ ] = TurbineStatusMessage( tit->first, TurbineStatus( *tit->second, unit.cascadeBurnIn() ) );
    }
    LOG_DEBUG2( "#" << unit.index() << ": sending turbines statuses, it=" << unit.iteration() << ", sz=" << msgs.size() );

    if ( _unitStatusISend ) _unitStatusISend.reset();
    _unitStatusISend.reset( new turbine_status_isend_type( _comm, _sampleCollectorRank, MTAG_TurbineStatuses, msgs ) );
}

template<typename ParticleEnergyEval>
template<class CascadeUnit>
void MPIUnitCommunicator<ParticleEnergyEval>::broadcastUnitStatuses(
    const CascadeUnit&  unit
){
    typedef CascadeUnit cascade_unit_type;
    std::vector<TurbineStatusMessage> msgs( unit.turbines().size() + _externalTurbineStatuses.size() );
    int i = 0;
    // own turbines
    for ( typename cascade_unit_type::const_turbine_iterator tit = unit.turbines().begin(); tit != unit.turbines().end(); ++tit ) {
        msgs[ i++ ] = TurbineStatusMessage( tit->first, TurbineStatus( *tit->second, unit.cascadeBurnIn() ) );
    }
    // external turbines
    for ( turbine_status_map_type::const_iterator tsit = _externalTurbineStatuses.begin(); tsit != _externalTurbineStatuses.end(); ++tsit ) {
        msgs[ i++ ] = TurbineStatusMessage( tsit->first, tsit->second );
    }
    LOG_DEBUG2( "#" << unit.index() << ": sending turbines statuses, it=" << unit.iteration() << ", sz=" << msgs.size() );

    if ( _unitStatusesIBroadcast ) _unitStatusesIBroadcast.reset();
    _unitStatusesIBroadcast.reset( new turbine_status_ibroadcast_type( _comm, MTAG_TurbineStatuses, msgs ) );
}

template<typename ParticleEnergyEval>
void MPIUnitCommunicator<ParticleEnergyEval>::broadcastShutdown(
    turbine_ix_t originIx
){
    for ( int i = 0; i < _comm.size(); i++ ) {
        if ( i != _comm.rank() ) {
            /// @todo use of ibroadcast might be required
            _comm.isend( i, MTAG_Shutdown, originIx );
        }
    }
}

template<typename ParticleEnergyEval>
template<class ParticlesIterator>
void MPIUnitCommunicator<ParticleEnergyEval>::sendParticles(
    int                         destinationRank,
    int                         tag,
    turbine_ix_t                destinationTurbineIx,
    turbine_ix_t                originTurbineIx,
    const ParticlesIterator&    particlesBegin,
    const ParticlesIterator&    particlesEnd
){
    typedef ParticleMessage<particle_type> particle_message_type;
    std::vector<particle_message_type> msgs( particlesEnd - particlesBegin );
    LOG_DEBUG1( "#" << rank() << ": sending " << msgs.size() << " particles from turbine #" << originTurbineIx 
                << " to turbine #" << destinationTurbineIx << " in unit #" << destinationRank << ", tag " << tag );
    for ( ParticlesIterator it = particlesBegin; it != particlesEnd; ++it ) {
        msgs[ it - particlesBegin ] = particle_message_type( originTurbineIx, destinationTurbineIx, *it );
    }
    _particleISendRequests.push_back( new particle_isend_request_type( _comm, destinationRank, tag,
                                                                       msgs ) );
}

template<typename ParticleEnergyEval>
void MPIUnitCommunicator<ParticleEnergyEval>::sendSample(
    turbine_ix_t            originIx,
    const particle_type&    sample
){
    if ( _sampleCollectorRank != rank() ) {
        sendParticles( _sampleCollectorRank, MTAG_Sample, TURBINE_NA, originIx, &sample, &sample + 1 );
    }
}

template<typename ParticleEnergyEval>
boost::optional<TurbineStatus> MPIUnitCommunicator<ParticleEnergyEval>::lastTurbineStatus(
    turbine_ix_t    turbineIx
) const {
    turbine_status_map_type::const_iterator sit = _externalTurbineStatuses.find( turbineIx );
    boost::optional<TurbineStatus> res;
    if ( sit != _externalTurbineStatuses.end() ) {
        res = sit->second;
    }
    return ( res );
}

template<typename ParticleEnergyEval>
MPIUnitCommunicator<ParticleEnergyEval>::MPIEnergyLandscapeRequestStorage::MPIEnergyLandscapeRequestStorage(
    communicator_type&  comm,
    turbine_ix_t        senderIx
#if defined(MPI_TIMEOUT_CONTROL)
    , double            timeout
#endif
) : super( comm.communicator(), MTAG_EnergyLandscapeRequest, senderIx
#if defined(MPI_TIMEOUT_CONTROL)
    , timeout
#endif
) , _waitingUnits( super::waitingUnits() )
  , _landscapes( _comm.size() )
  , _landscapeResponseRequests( _comm.size() )
#if defined(MPI_TIMEOUT_CONTROL)
  , _lastReported( 0 )
#endif
{
    _waitingUnits.set( _comm.rank(), false );
    foreach_bit( size_t, i, _waitingUnits ) {
        _landscapeResponseRequests[ i ] = _comm.irecv( i, MTAG_EnergyLandscapeResponse, _landscapes[ i ] );
    }
}

template<typename ParticleEnergyEval>
MPIUnitCommunicator<ParticleEnergyEval>::MPIEnergyLandscapeRequestStorage::~MPIEnergyLandscapeRequestStorage()
{
    LOG_DEBUG2( "MPIEnergyLandscapeRequestStorage::dtor()" );
    foreach_bit( size_t, i, _waitingUnits ) {
        _landscapeResponseRequests[i].cancel();
    }
    LOG_DEBUG2( "All landscape requests are shut down" );
}

template<typename ParticleEnergyEval>
bool MPIUnitCommunicator<ParticleEnergyEval>::MPIEnergyLandscapeRequestStorage::complete()
{
    LOG_DEBUG2( "#" << _comm.rank() << ": testing landscape status" );
    if ( !super::complete() )   return ( false );
    if ( _waitingUnits.none() ) return ( true );

    foreach_bit( size_t, i, _waitingUnits ) {
        boost::optional<mpi::status> status = _landscapeResponseRequests[ i ].test();
        if ( status ) {
            LOG_DEBUG1( "#" << _comm.rank() << ": got response from " << status->source() << ", " << _landscapes[i].size() << " particles" );
            TurbineEnergyLadder::AppendLandscape( _cumulativeLandscape, _landscapes[i] );
            LOG_DEBUG2( "#" << _comm.rank() << ": erasing rank #" << status->source() << " from waiting" );
            _waitingUnits.set( i, false );
        }
    }
    if ( _waitingUnits.none() ) {
        LOG_DEBUG1( "#" << _comm.rank() << ": landscapes collected" );
        return ( true );
    } else {
#if defined(MPI_TIMEOUT_CONTROL)
        if ( _timeout > 0 && elapsed() > _lastReported + _timeout ) {
            std::ostringstream strstr;
            foreach_bit( size_t, i, _waitingUnits ) {
                strstr << ' ' << i;
            }
            LOG_WARN( "#" << _comm.rank() << ": MPIEnergyLandscapeRequestStorage timeout " 
                          << boost::format("%.0f") % elapsed() << "s still waiting for " << strstr.str() << "..." );
            _lastReported = elapsed();
        }
#endif
        return ( false );
    }
}

/******************************************************************************
 * MPITurbineConnection implementation
 ******************************************************************************/
