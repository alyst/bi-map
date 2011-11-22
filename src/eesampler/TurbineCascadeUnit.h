#pragma once

#include "../BasicTypedefs.h"

#include <boost/unordered_map.hpp>

#include <boost/format.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/timer.hpp>
#include <boost/ptr_container/ptr_unordered_map.hpp>
#include <boost/math/special_functions/fpclassify.hpp>

#include "Turbine.h"
#include "ITurbineConnection.h"

typedef size_t tcascadeunit_ix_t;

typedef struct {
    turbine_ix_t        turbineIx;
    tcascadeunit_ix_t   unitIx;
    size_t              level;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( turbineIx );
        ar & BOOST_SERIALIZATION_NVP( unitIx );
        ar & BOOST_SERIALIZATION_NVP( level );
    }
} TurbineNode;

typedef std::vector<TurbineNode> tcascade_structure_type;

tcascade_structure_type createTurbineCascadeStructure( const std::vector<size_t>& turbinesPerLevel, size_t units );
tcascade_structure_type createTurbineCascadeStructure( size_t levels, size_t turbines, size_t units );

#define UNIT_LOG_PREFIX( unit ) "#" << (unit).index() << "(" << ( boost::format("%.0f") % (unit).elapsed() ) << "s): "

/**
    Equi-energy sampler parameters.
 */
struct TurbineCascadeParams {
    size_t                  turbinesCount;          /** number of turbines in the cascade */
    std::vector<size_t>     turbinesPerLevel;       /** number of turbines per level */

    // energy ladder related
    size_t      levelsCount;            /** number of levels in the cascade */
    energy_type temperatureMultiplier;  /** ratio between temperatures of subsequent
                                            hotter and colder turbines */

    size_t      burnInIterations;       /** cascade burn-in iterations (iterations of the leading turbine) */
    size_t      ladderAdjustPeriod;     /** period (iterations of the leading turbine) the energy */
    size_t      broadcastStatusPeriod;  /** in how many iterations to broadcast unit's status */

    // single-turbine parameters
    TurbineParams   turbineParams;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( turbinesCount );
        ar & BOOST_SERIALIZATION_NVP( levelsCount );
        ar & BOOST_SERIALIZATION_NVP( turbinesPerLevel );
        ar & BOOST_SERIALIZATION_NVP( temperatureMultiplier );
        ar & BOOST_SERIALIZATION_NVP( burnInIterations );
        ar & BOOST_SERIALIZATION_NVP( ladderAdjustPeriod );
        ar & BOOST_SERIALIZATION_NVP( broadcastStatusPeriod );
        ar & BOOST_SERIALIZATION_NVP( turbineParams );
    }

    TurbineCascadeParams();
};

struct TurbineEnergyLadder {
    /// @todo move out, together with separate namespace creation
    typedef std::vector<energy_type> energy_vector_type;
    typedef std::map<size_t, energy_vector_type> energy_landscape_type;

    std::vector<EnergyTransform>    steps;          /// steps of energy ladder

    TurbineEnergyLadder() {};

    TurbineEnergyLadder( const TurbineCascadeParams& params,
                         const energy_landscape_type& energyLandscape );

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( steps );
    }

    template<typename PreEnergyType>
    static void AppendLandscape( std::map<size_t, std::vector<PreEnergyType> >& dest,
                                 const std::map<size_t, std::vector<PreEnergyType> >& src )
    {
        typedef std::map<size_t, std::vector<PreEnergyType> > pre_energy_landscape_type;
        // merge local and external landscapes
        for ( typename pre_energy_landscape_type::const_iterator sit = src.begin();
              sit != src.end(); ++sit
        ){
            typename pre_energy_landscape_type::iterator dit = dest.find( sit->first );
            if ( dit != dest.end() ) {
                dit->second.insert( dit->second.end(),
                                    sit->second.begin(), sit->second.end() );
            } else {
                dest.insert( dit, *sit );
            }
        }
    }
};

#if 0
class IParticleCollector {
public:
    virtual bool storeSample( double time, turbine_ix_t originIx, const IParticle& sample ) = 0;
};
#endif

/**
 *  Do-nothing implementation of UnitCommunicator concept.
 **/
template<class ParticleEnergyEval>
struct DummyUnitCommunicator {
    typedef ParticleEnergyEval energy_eval_type;
    typedef typename energy_eval_type::particle_type particle_type;

    typedef typename particle_type::pre_energy_type pre_energy_type;
    typedef std::vector<pre_energy_type> pre_energy_vector_type;
    typedef std::map<size_t, pre_energy_vector_type> pre_energy_landscape_type;

    typedef boost::timer timer_type;
    typedef TurbineEnergyLadder::energy_landscape_type energy_landscape_type;

    struct DummyLandscapeRequestStatus {
        bool complete() const { return ( true ); }
        pre_energy_landscape_type energyLandscape() const {
            return ( pre_energy_landscape_type() );
        }
    };

    typedef DummyLandscapeRequestStatus landscape_request_status_type;
    typedef DummyTurbineConnection<energy_eval_type> turbine_connection_type;

    template<class CascadeUnit>
    void sendUnitStatus( const CascadeUnit& unit ) {};
    template<class CascadeUnit>
    void broadcastUnitStatuses( const CascadeUnit& unit ) {};
    void broadcastShutdown( turbine_ix_t originIx ) {};
    void broadcastEnergyLadder( turbine_ix_t originIx, const TurbineEnergyLadder& ladder,
                                const energy_eval_type& energyEval, const particle_type& iniParticle ) {};
    template<class CascadeUnit>
    bool processExternalEvents( CascadeUnit& unit, bool wait ) { return ( false ); };
    template<class CascadeUnit>
    landscape_request_status_type requestEnergyLandscape( const CascadeUnit& unit, turbine_ix_t senderIx ) {
            return ( landscape_request_status_type() );
    }
    turbine_connection_type* createTurbineConnection( const TurbineNode& externalTurbine ) {
        return ( NULL );
    }
    void sendSample( turbine_ix_t originIx, const particle_type& sample ) {};
};

/**
    Execution unit of turbine cascade.
    Unit executes its turbines sequentially.
    Base class provides interface to implement
    messaging between parallel units (using e.g. MPI).
 */
template<class ParticleCollector, class DynamicParticleFactory, 
         class ParticleGenerator = SunkenParticleGenerator<typename DynamicParticleFactory::static_particle_energy_eval_type>, 
         class UnitCommunicator = DummyUnitCommunicator<typename DynamicParticleFactory::static_particle_energy_eval_type> >
class TurbineCascadeUnit {
public:
    typedef boost::dynamic_bitset<> units_mask_type;
    typedef DynamicParticleFactory factory_type;
    typedef typename DynamicParticleFactory::dynamic_particle_type dynamic_particle_type;
    typedef typename DynamicParticleFactory::static_particle_type particle_type;
    typedef boost::scoped_ptr<dynamic_particle_type> dynamic_particle_pointer_type;
    typedef UnitCommunicator communicator_type;
    typedef typename communicator_type::landscape_request_status_type landscape_request_status_type;
    typedef CascadeParticle<particle_type> cascade_particle_type;
    typedef typename communicator_type::timer_type timer_type;

    typedef ParticleGenerator generator_type;
    typedef ParticleCollector collector_type;

    class GenericTurbineConnector;
    friend class GenericTurbineConnector;

    typedef boost::dynamic_bitset<>  turbine_ix_set_type;
    typedef GenericTurbineConnector connector_type;
    typedef ParticleTurbine<factory_type, generator_type, connector_type> turbine_type;

    typedef typename connector_type::jump_request_status_type jump_request_status_type;

    typedef typename DynamicParticleFactory::static_particle_energy_eval_type particle_energy_eval_type;
    typedef typename particle_energy_eval_type::energy_type energy_type;
    typedef std::vector<energy_type> energy_vector_type;
    typedef std::map<size_t, energy_vector_type> energy_landscape_type;

    typedef typename particle_type::pre_energy_type pre_energy_type;
    typedef std::vector<pre_energy_type> pre_energy_vector_type;
    typedef std::map<size_t, pre_energy_vector_type> pre_energy_landscape_type;

    typedef TurbineCascadeUnit<collector_type, factory_type, generator_type, communicator_type> cascade_unit_type;

    typedef boost::ptr_unordered_map<turbine_ix_t, turbine_type> turbine_container_type;
    typedef typename turbine_container_type::iterator turbine_iterator;
    typedef typename turbine_container_type::const_iterator const_turbine_iterator;

private:
    const TurbineCascadeParams          _params;
    const tcascade_structure_type       _cascadeStructure;  /// structure of turbines cascade communication
    units_mask_type                     _coldestUnitsMask;   /// bitset with units containing 0-level turbines

    const tcascadeunit_ix_t             _index;
    turbine_ix_t                        _leadTurbineIx;
    const gsl_rng*                      _rng;

    const factory_type&                 _dynamicParticleFactory;
    const generator_type*               _pGenerator;
    collector_type*                     _pCollector; /** Instance of collector concept to collect the samples.
                                                         Only one unit can have collector interface */
    TurbineCascadeExecutionMonitor*     _pMonitor;

    communicator_type*                  _pUnitComm;

    timer_type                          _timer;
    size_t                              _iteration;

    typedef boost::unordered_map<size_t, std::vector<turbine_type*> > turbines_by_level_map_type;

    turbine_container_type              _turbines;          /// turbines of the unit
    turbines_by_level_map_type          _turbinesByLevel;   /// turbines grouped by level

    typedef typename communicator_type::turbine_connection_type connection_type;
    typedef boost::ptr_unordered_map<turbine_ix_t, connection_type> connection_container_type;
    typedef boost::unordered_map<size_t, std::vector<connection_type*> > connections_by_level_map_type;

    connection_container_type               _connections;           /// connections to external turbines
    connections_by_level_map_type           _connectionsByLevel;    /// connections grouped by level

    typedef std::vector<turbine_ix_set_type>   turbine_ixs_by_level_map_type;
    turbine_ixs_by_level_map_type           _jumpSourcesByLevel;    /// turbines that be used for EE jumps, by level
    turbine_ixs_by_level_map_type           _particleDestinationsByLevel;    /// turbines that be used for as destination for generated particles, by level

    bool                                    _cascadeBurnIn; /// cascade is burning in
    bool                                    _shuttingDown;  /// cascade is shutting down

    bool processExternalEvents( bool wait );
    connection_type* createTurbineConnection( const TurbineNode& externalTurbine );

    turbine_type* appendTurbine( const TurbineNode& node, const EnergyTransform& transform );

    TurbineCascadeUnit( const cascade_unit_type& that );

public:
    TurbineCascadeUnit( tcascadeunit_ix_t index, const TurbineCascadeParams& params,
                        const tcascade_structure_type& cascadeStructure, const gsl_rng* rng, 
                        const factory_type& dynParticleFactory, 
                        const generator_type* pGenerator, collector_type* pCollector, communicator_type* unitComm = NULL );
    void init( const particle_type& start );
    void run();
    void iterate();

    void processTurbineCreated( turbine_ix_t turbineIx );
    void processBurnInOver();
    void adjustCascadeEnergyLadder( turbine_ix_t turbineIx );
    void adjustUnitEnergyLadder( turbine_ix_t originIx, const TurbineEnergyLadder& ladder,
                                 const particle_energy_eval_type& energyEval,
                                 const particle_type& iniParticle );
    bool sampleParticle( turbine_ix_t turbineIx, const particle_type& particle );

    bool storesSamples() const {
        return ( _pCollector != NULL );
    }

    tcascadeunit_ix_t index() const {
        return ( _index );
    };
    const tcascade_structure_type& cascadeStructure() const {
        return ( _cascadeStructure );
    }
    const units_mask_type& coldestUnitsMask() const {
        return ( _coldestUnitsMask );
    }
    bool isLeading() const {
        return ( _cascadeStructure[ _leadTurbineIx ].unitIx == _index );
    }
    bool cascadeBurnIn() const {
        return ( _cascadeBurnIn );
    }

    const turbine_container_type& turbines() const {
        return ( _turbines );
    }
    turbine_container_type& turbines() {
        return ( _turbines );
    }
    size_t iteration() const {
        return ( _iteration );
    }
    double elapsed() const {
        return ( _timer.elapsed() );
    }
    pre_energy_landscape_type localEnergyLandscape() const;
};

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator>
TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, 
                   DummyUnitCommunicator<typename DynamicParticleFactory::static_particle_energy_eval_type> >*
CreateSingleUnitCascade( const TurbineCascadeParams& params, const gsl_rng* rng, 
                         const DynamicParticleFactory& dynParticleFactory, 
                         const ParticleGenerator* pGenerator, ParticleCollector& pCollector )
{
    return ( new TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, 
                                DummyUnitCommunicator<typename DynamicParticleFactory::static_particle_energy_eval_type> >( 0, params,
                                 createTurbineCascadeStructure( params.levelsCount, params.turbinesCount, 1 ),
                                 rng, dynParticleFactory, pGenerator, &pCollector ) );
}

/**
 *  Probability to use turbine of the same or lower levels for EE-jump
 *  during burn-in.
 */
#define CONNECTOR_LOWER_TURBINE_JUMP_PROB 0.1

/**
 *  Minimal number of jumpable nodes on level to all jumping between them.
 */
#define CONNECTOR_MIN_BROTHERS_LOOP 4

template<class ParticleCollector, class DynamicParticleFactory,
         class ParticleGenerator, class UnitCommunicator>
class TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::
GenericTurbineConnector
{
private:
    typedef TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator> cascade_unit_type;
    typedef ParticleTurbineBase<factory_type> base_turbine_type;
    typedef typename base_turbine_type::static_particle_energy_eval_type static_particle_energy_eval_type;
    typedef typename cascade_unit_type::turbine_ix_set_type turbine_ix_set_type;
    typedef typename communicator_type::turbine_connection_type connection_type;
    typedef typename connection_type::jump_request_status_type external_jump_request_status_type;

    const TurbineNode&          _nodeInfo; /// node that connects to (colder, one that requests EE jumps)
    cascade_unit_type&          _unit;     /// parent turbine cascade unit

public:
    struct GenericJumpRequestStatus {
        typedef typename cascade_unit_type::cascade_particle_type cascade_particle_type;
        typedef typename cascade_unit_type::connection_type connection_type;

        const EEJumpParams                                      _params;
        boost::optional<cascade_particle_type>                  _particle;
        boost::optional<external_jump_request_status_type>      _external;

        bool complete() {
            return ( _external ? _external->complete() : true );
        }
        const boost::optional<cascade_particle_type> particle() const {
            return ( _external ? _external->particle() : _particle );
        }

        GenericJumpRequestStatus( connector_type& conn, turbine_ix_t initiatorIx, turbine_ix_t processorIx,
                                  const static_particle_energy_eval_type& energyEval,
                                  const EEJumpParams& params )
        : _params( params )
        {
            typename turbine_container_type::const_iterator tit = conn._unit._turbines.find( processorIx );
            if ( tit != conn._unit._turbines.end() ) {
                // try local ee-jump
                _particle = tit->second->processJumpRequest( _params, energyEval, false );
            }
            else {
                typename connection_container_type::iterator cit = conn._unit._connections.find( processorIx );
                if ( cit != conn._unit._connections.end() ) {
                    // create external request
                    _external.reset( cit->second->requestJump( initiatorIx, energyEval, params ) );
                }
                else {
                    THROW_RUNTIME_ERROR( "Turbine #" << processorIx << " not found in unit or its connections" );
                }
            }
        }
    };
    typedef GenericJumpRequestStatus jump_request_status_type;

    turbine_ix_set_type eeJumpSources() const {
        size_t jumpLevel = _unit.cascadeBurnIn() /** during burn-in can jump to the turbines of the same or lower levels */
                         ? _nodeInfo.level + 1 - std::min( _nodeInfo.level + 1, 
                                                           (size_t)gsl_ran_geometric( _unit._rng,
                                                                                      1.0 - CONNECTOR_LOWER_TURBINE_JUMP_PROB ) - 1 )
                         : _nodeInfo.level + 1;
        turbine_ix_set_type res = _unit._jumpSourcesByLevel[ jumpLevel ];
        // don't let self-jumps
        if ( jumpLevel == _nodeInfo.level ) {
            res.set( _nodeInfo.turbineIx, false );
            // introduce ring topology of jumps between turbines of the same level
            // if there's enough turbines to avoid autocorrelation loop
            const turbine_ix_set_type& brotherNodes = _unit._jumpSourcesByLevel[ _nodeInfo.level ];
            if ( brotherNodes.count() >= CONNECTOR_MIN_BROTHERS_LOOP )
            {
                // find next jumpable turbine
                turbine_ix_t nextTurbineIx = brotherNodes.find_next( _nodeInfo.turbineIx );
                if ( nextTurbineIx == turbine_ix_set_type::npos ) nextTurbineIx = brotherNodes.find_first();
                res.set( nextTurbineIx, true );
            }
        }
        return ( res );
    }

    const turbine_ix_set_type& particleDestinations() const {
        return ( _unit._particleDestinationsByLevel[ _nodeInfo.level + 1 ] );
    }

    turbine_ix_t randomTurbineIx(
        const turbine_ix_set_type& turbineIxs
    ) const {
        if ( turbineIxs.none() ) return ( TURBINE_NA );
        turbine_ix_t res = turbineIxs.find_first();
        if ( turbineIxs.count() == 1 ) return ( res );
        size_t ix = gsl_rng_uniform_int( _unit._rng, turbineIxs.count() );
        while ( ix-- ) {
            res = turbineIxs.find_next( res );
            BOOST_ASSERT( res != turbine_ix_set_type::npos );
        }
        return ( res );
    }

    turbine_ix_t randomEeJumpSource() const {
        return ( randomTurbineIx( eeJumpSources() ) );
    }

    turbine_ix_t randomParticleDestination() const {
        return ( randomTurbineIx( particleDestinations() ) );
    }

    /**
     *  Counts of prisoners in the connected turbine (one that of higher temperature)
     */
    size_t prisonersCount( turbine_ix_t turbineIx ) const {
        // unit-local turbines
        typename turbine_container_type::const_iterator ltit = _unit._turbines.find( turbineIx );
        if ( ltit != _unit._turbines.end() ) {
            BOOST_ASSERT( ltit->second->level() == _nodeInfo.level + 1 );
            return ( ltit->second->prisonersCount() );
        }
        // external turbines
        typename connection_container_type::const_iterator ctit = _unit._connections.find( turbineIx );
        if ( ctit != _unit._connections.end() ) {
            BOOST_ASSERT( ctit->second->level() == _nodeInfo.level + 1 );
            return ( ctit->second->prisonersCount() );
        }
        else {
            THROW_RUNTIME_ERROR( "turbine #" << turbineIx << " not found" );
        }
    }

    jump_request_status_type requestJump( turbine_ix_t initiatorIx, turbine_ix_t processorIx,
                                          const static_particle_energy_eval_type& energyEval,
                                          const EEJumpParams& params ) {
        return ( jump_request_status_type( *this, initiatorIx, processorIx,
                                           energyEval, params ) );
    }

    template<class ParticleIterator>
    bool sendParticles( turbine_ix_t originIx, turbine_ix_t destinationIx,
                        const ParticleIterator& particlesBegin,
                        const ParticleIterator& particlesEnd )
    {
        typename turbine_container_type::iterator ltit = _unit._turbines.find( destinationIx );
        if ( ltit != _unit._turbines.end() ) {
            BOOST_ASSERT( ltit->second->level() >= _nodeInfo.level );
            return ( ltit->second->particleCache()->push_many( particlesBegin, particlesEnd,
                                                            _unit._timer.elapsed(), _unit._iteration ) );
        }
        typename connection_container_type::iterator ctit = _unit._connections.find( destinationIx );
        if ( ctit != _unit._connections.end() ) {
            BOOST_ASSERT( ctit->second->level() >= _nodeInfo.level );
            return ( ctit->second->sendParticles( originIx, particlesBegin, particlesEnd ) );
        }
        else {
            THROW_RUNTIME_ERROR( "turbine #" << destinationIx << " not found" );
        }
    }

    bool processExternalEvents( cascade_unit_type& unit, bool wait ) {
        BOOST_ASSERT( &unit == &_unit );
        return ( _unit.processExternalEvents( wait ) );
    }

    GenericTurbineConnector( cascade_unit_type& unit, const TurbineNode& nodeInfo )
    : _nodeInfo( nodeInfo ), _unit( unit )
    {
    }
};

/******************************************************************************
 * TurbineCascadeUnit Implementation
 ******************************************************************************/

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::TurbineCascadeUnit(
    tcascadeunit_ix_t               index,
    const TurbineCascadeParams&     params,
    const tcascade_structure_type&  cascadeStructure,
    const gsl_rng*                  rng,
    const factory_type&             factory,
    const generator_type*           pGenerator,
    collector_type*                 pCollector,
    communicator_type*              pUnitComm
) : _params( params )
  , _cascadeStructure( cascadeStructure )
  , _index( index )
  , _leadTurbineIx( TURBINE_NA )
  , _rng( rng )
  , _dynamicParticleFactory( factory )
  , _pGenerator( pGenerator )
  , _pCollector( pCollector )
  , _pMonitor( NULL )
  , _pUnitComm( pUnitComm )
  , _iteration( 0 )
  , _cascadeBurnIn( true )
  , _shuttingDown( false )
{
    // find how many levels in the cascade
    size_t maxLevel = 0;
    size_t maxUnit = _index;
    for ( turbine_ix_t turbineIx = 0; turbineIx < _cascadeStructure.size(); ++turbineIx ) {
        const TurbineNode& node = _cascadeStructure[ turbineIx ];
        if ( node.unitIx > maxUnit ) maxUnit = node.unitIx;
        if ( node.level > maxLevel ) maxLevel = node.level;
    }
    _coldestUnitsMask = units_mask_type( maxUnit + 1 );
    for ( turbine_ix_t turbineIx = 0; turbineIx < _cascadeStructure.size(); ++turbineIx ) {
        const TurbineNode& node = _cascadeStructure[ turbineIx ];
        if ( node.level == 0 ) _coldestUnitsMask.set( node.unitIx );
    }
    if ( _coldestUnitsMask.none() ) {
        THROW_EXCEPTION( std::invalid_argument, "Turbine cascade doesn't contain any T=1 turbines" );
    }
    _leadTurbineIx = _coldestUnitsMask.find_first(); // lead turbine is the first coldest turbine
    // initialize available jump sources cache (with one extra level)
    _jumpSourcesByLevel.resize( maxLevel + 2, turbine_ix_set_type( _cascadeStructure.size() ) );
    _particleDestinationsByLevel.resize( maxLevel + 2, turbine_ix_set_type( _cascadeStructure.size() ) );
    LOG_DEBUG1( UNIT_LOG_PREFIX( *this ) << "Turbine Cascade Unit #" << index << " created" );
}

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
void TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::init(
    const particle_type& start
){
    if ( !_turbines.empty() ) {
        THROW_RUNTIME_ERROR( UNIT_LOG_PREFIX( *this ) << "unit turbines already initialized" );
    }
    // create all turbines of this unit
    for ( turbine_ix_t turbineIx = 0; turbineIx < _cascadeStructure.size(); ++turbineIx ) {
        const TurbineNode& nodeInfo = _cascadeStructure[ turbineIx ];
        if ( nodeInfo.unitIx == _index ) {
            turbine_type* pTurbine = appendTurbine( nodeInfo, EnergyTransform() );
            pTurbine->setMovingParticle( start );
        }
    }
    if ( _turbines.empty() ) {
        LOG_WARN( UNIT_LOG_PREFIX( *this ) << "cascade structure had no turbines for this unit" );
    }
}

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
void TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::run()
{
    while ( !_shuttingDown ) {
        iterate();
        LOG_DEBUG1_IF( _shuttingDown, UNIT_LOG_PREFIX( *this ) << ": shutting down" );
    }
}

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
bool TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::processExternalEvents(
    bool wait
){
    bool res = _shuttingDown;
    if ( !_shuttingDown && _pUnitComm ) {
        if ( ( res = _pUnitComm->processExternalEvents( *this, wait ) ) ) {
            _shuttingDown = true;
        }
    }
    return ( res );
}

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
typename TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::connection_type*
TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::createTurbineConnection(
    const TurbineNode&  externalTurbine
){
    LOG_INFO( UNIT_LOG_PREFIX( *this ) << "creating connection to turbine #" << externalTurbine.turbineIx << ", level " << externalTurbine.level );
    BOOST_ASSERT( externalTurbine.unitIx != _index );
    if ( !_pUnitComm ) THROW_RUNTIME_ERROR( "Cannot connect to external turbine #" << externalTurbine.turbineIx 
                                             << ": no communicator" );
    std::auto_ptr<connection_type> connection( _pUnitComm->createTurbineConnection( externalTurbine ) );
    std::pair<typename connection_container_type::iterator, bool> connectionAddRes = _connections.insert( externalTurbine.turbineIx, connection );
    if ( !connectionAddRes.second ) THROW_RUNTIME_ERROR( "Connection to turbine #" << externalTurbine.turbineIx 
                                             << " not added to connections map of unit #" << index() );
    if ( _connectionsByLevel.find( externalTurbine.level ) == _connectionsByLevel.end() ) {
        _connectionsByLevel[ externalTurbine.level ] = std::vector<connection_type*>();
    }
    _connectionsByLevel[ externalTurbine.level ].push_back( connectionAddRes.first->second );
    return ( connectionAddRes.first->second );
}

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
typename TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::turbine_type*
TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::appendTurbine(
    const TurbineNode&      nodeInfo,
    const EnergyTransform&  energyTransform
){
    LOG_INFO( UNIT_LOG_PREFIX( *this ) << "creating turbine #" << nodeInfo.turbineIx << ", level " << nodeInfo.level );
    BOOST_ASSERT( nodeInfo.unitIx == _index );
    std::auto_ptr<turbine_type> turbine( new turbine_type( nodeInfo.level, nodeInfo.turbineIx, _params.turbineParams, energyTransform, 
            _rng, _dynamicParticleFactory, _pGenerator, GenericTurbineConnector( *this, nodeInfo ),
            gsl_rng_uniform_int( _rng, _params.turbineParams.particleSnapshotPeriod ),
            gsl_rng_uniform_int( _rng, _params.turbineParams.particleSamplingPeriod ),
            NULL, NULL ) );
    std::pair<turbine_iterator, bool> turbineAddRes = _turbines.insert( nodeInfo.turbineIx, turbine );
    if ( !turbineAddRes.second ) THROW_RUNTIME_ERROR( "Turbine #" << nodeInfo.turbineIx << " not added to turbines map of unit #" << index() );
    if ( _turbinesByLevel.find( nodeInfo.level ) == _turbinesByLevel.end() ) {
        _turbinesByLevel[ nodeInfo.level ] = std::vector<turbine_type*>();
    }
    _turbinesByLevel[ nodeInfo.level ].push_back( turbineAddRes.first->second );
    processTurbineCreated( nodeInfo.turbineIx );
    return ( turbineAddRes.first->second );
}

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
void TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::iterate()
{
    if ( _pMonitor ) _pMonitor->notifyIteration( _iteration );
    if ( _turbines.empty() ) {
        if ( !_pUnitComm ) {
            THROW_RUNTIME_ERROR( UNIT_LOG_PREFIX( *this ) << "no turbines! Have you forgotten to call init()?" );
        } else {
            LOG_WARN( UNIT_LOG_PREFIX( *this ) << "unit has no turbines! Have you forgotten to call init()?" );
        }
    }
    // iterate each turbine
    for ( turbine_iterator turbineIt = _turbines.begin(); turbineIt != _turbines.end(); ++turbineIt ) {
        // process external events
        processExternalEvents( false );
        if ( _shuttingDown ) break;

        turbine_ix_t turbineIx = turbineIt->first;
        turbine_type& turbine = *turbineIt->second;
        turbine.iterate( *this );

        if ( turbineIx == _leadTurbineIx && _cascadeBurnIn ) {
            if ( turbine.iteration() > 0
                 && ( turbine.iteration() % _params.ladderAdjustPeriod == 0 )
            ){
                adjustCascadeEnergyLadder( turbineIx );
            }
            if ( turbine.iteration() >= _params.burnInIterations ) {
                processBurnInOver();
            }
            else if ( turbine.iteration() % TURBINE_LOG_REPORT_PERIOD == 0 ) {
                LOG_INFO( UNIT_LOG_PREFIX( *this ) << "leading turbine #" << turbineIx << ": " 
                          << _params.burnInIterations - turbine.iteration() << " iterations to go..." );
            }
        }
        if ( _pUnitComm && ( _iteration % _params.broadcastStatusPeriod ) == 0 ) {
            if ( isLeading() ) {
                _pUnitComm->broadcastUnitStatuses( *this );
            } else {
                _pUnitComm->sendUnitStatus( *this );
            }
        }

        // send particles to the sample storage
        if ( ( turbine.level() == 0 ) && !cascadeBurnIn()
            && ( ( ( turbine.iteration() - 1 ) + turbine.particleSamplingOffset() ) % _params.turbineParams.particleSamplingPeriod == 0 )
        ){
            if ( sampleParticle( turbineIx, turbine.movingParticle() ) ) {
                break;
            }
        }
    }
    processExternalEvents( _turbines.empty() );
    _iteration++;
}

/**
 *  Initiates adjustments of the whole cascade energy ladder.
 *  Collects energies of the samples from all the turbines.
 */
template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
void TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::
adjustCascadeEnergyLadder(
    turbine_ix_t turbineIx  /** turbine initiating adjustment,
                                will collect info from other turbines */
){
    if ( _pUnitComm && !isLeading() ) _pUnitComm->sendUnitStatus( *this );
    const turbine_type& turbine = *_turbines.find( turbineIx )->second;
    LOG_INFO( UNIT_LOG_PREFIX( *this ) << "leading turbine #" << turbineIx << ": "
              << turbine.iteration() << " iterations done, starting ladder adjustment" );

    pre_energy_landscape_type preEnergyLandscape = localEnergyLandscape();
    if ( _pUnitComm ) {
        // send energy landscape request
        LOG_INFO( UNIT_LOG_PREFIX( *this ) << "turbine #" << turbineIx << " requesting energy landscape" );
        landscape_request_status_type status = _pUnitComm->requestEnergyLandscape( *this, turbineIx );
        double firstNotifyTime = elapsed();
        double lastNotifyTime = firstNotifyTime;
        // wait until landscape request complete
        while ( !status.complete() ) {
            processExternalEvents( false );
            // check for stuck landscape
            if ( elapsed() > lastNotifyTime + 60 ) {
                lastNotifyTime = elapsed();
                LOG_WARN( UNIT_LOG_PREFIX( *this ) << "waiting for energy landscape, " 
                                                   << boost::format("%.0f") % ( lastNotifyTime - firstNotifyTime ) <<  "s elapsed..." );
            }
            if ( _shuttingDown ) return;
        }
        // merge local and external landscapes
        TurbineEnergyLadder::AppendLandscape( preEnergyLandscape, status.energyLandscape() );
    }
    // convert to energies and sort
    energy_landscape_type energyLandscape;
    for ( typename pre_energy_landscape_type::iterator lit = preEnergyLandscape.begin();
          lit != preEnergyLandscape.end(); ++lit
    ){
        energy_vector_type levelLandscape;
        for ( size_t i = 0; i < lit->second.size(); i++ ) {
            levelLandscape.push_back( turbine.energyEval()( lit->second[i] ) );
        }
        std::sort( levelLandscape.begin(), levelLandscape.end() );
        energyLandscape[ lit->first ] = levelLandscape;
        LOG_INFO( UNIT_LOG_PREFIX( *this ) << " landscape level #" << lit->first << ": " << levelLandscape.size() 
                    << " energy values in range [" << ( levelLandscape.size() > 0 ? levelLandscape.front() : std::numeric_limits<energy_type>::quiet_NaN() )
                    << ", " << ( levelLandscape.size() > 0 ? levelLandscape.back() : std::numeric_limits<energy_type>::quiet_NaN() ) << "]" );
    }
    if ( energyLandscape.size() == 0 ) {
        LOG_WARN( UNIT_LOG_PREFIX( *this ) << "skipping energy ladder formation, landscape is empty" );
        return;
    }

    // build energy ladder steps
    TurbineEnergyLadder ladder( _params, energyLandscape );

    particle_energy_eval_type newEnergyEval = turbine.energyEval();
    if ( preEnergyLandscape.count( 0 ) ) {
        LOG_INFO( UNIT_LOG_PREFIX( *this ) << "Updating energy evaluator" );
        newEnergyEval = _dynamicParticleFactory.updateEnergyEval( turbine.energyEval(), preEnergyLandscape[0] );
    }

    const particle_type& iniParticle = turbine.movingParticle();
    adjustUnitEnergyLadder( turbineIx, ladder, newEnergyEval, iniParticle );
    if ( _pUnitComm ) {
        LOG_DEBUG1( UNIT_LOG_PREFIX( *this ) << "broadcasting adjusted ladder" );
        // broadcast adjusted turbine params to all units
        _pUnitComm->broadcastEnergyLadder( turbineIx, ladder, newEnergyEval, iniParticle );
    }
}

/**
 *  Manage local copy of cascade probe
 */
template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
void TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::processTurbineCreated(
    turbine_ix_t turbineIx
){
    LOG_DEBUG1( UNIT_LOG_PREFIX( *this ) << (_turbines.find( turbineIx ) != _turbines.end() ? "local" : "external") << " turbine #" << turbineIx << " created" );
    const TurbineNode& newTurbine = _cascadeStructure[ turbineIx ];
    if ( newTurbine.unitIx != index() ) {
        LOG_DEBUG2( UNIT_LOG_PREFIX( *this ) << "checking if new turbine connection should be created" );
        for ( turbine_ix_t turbineIx = 0; turbineIx < _cascadeStructure.size(); ++turbineIx ) {
            const TurbineNode& nodeProps = _cascadeStructure[ turbineIx ];
            if ( ( nodeProps.unitIx == index() )
                && ( newTurbine.level <= nodeProps.level + 1 )
            ){
                createTurbineConnection( newTurbine );
            };
        }
    }
    _jumpSourcesByLevel[ newTurbine.level ].set( turbineIx );
    _particleDestinationsByLevel[ newTurbine.level ].set( newTurbine.turbineIx );
}

/**
 *  Manage local copy of cascade probe
 */
template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
void TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::processBurnInOver()
{
    if ( _cascadeBurnIn ) {
        _cascadeBurnIn = false;
        if ( isLeading() && _pUnitComm ) _pUnitComm->broadcastUnitStatuses( *this );
        LOG_INFO( UNIT_LOG_PREFIX( *this ) << "all turbines ready, starting sampling" );
    }
}

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
typename TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::pre_energy_landscape_type
TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::localEnergyLandscape() const
{
    pre_energy_landscape_type energyLandscape;

    for ( typename turbine_container_type::const_iterator it = _turbines.begin(); it != _turbines.end(); ++it ) {
        pre_energy_vector_type turbineEnergies = it->second->particleCache()->preEnergies();
        if ( !turbineEnergies.empty() ) {
            typename pre_energy_landscape_type::iterator lit = energyLandscape.find( it->second->level() );
            if ( lit == energyLandscape.end() ) {
                energyLandscape[ it->second->level() ] = turbineEnergies;
            } else {
                lit->second.insert( lit->second.end(), turbineEnergies.begin(), turbineEnergies.end() );
            }
        }
    }
    return ( energyLandscape );
}

/**
 *  Adjusts sampling transforms of the unit's turbines.
 */
template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
void TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::adjustUnitEnergyLadder(
    turbine_ix_t                originIx,       /** turbine initiated adjustment */
    const TurbineEnergyLadder&  ladder,         /** adjusted energy transforms */
    const particle_energy_eval_type& energyEval,
    const particle_type&        iniParticle     /** particle to initialize new turbines */
){
    LOG_DEBUG1( UNIT_LOG_PREFIX( *this ) << "adjusting energy ladder" );
    // adjust parameters of existing turbines
    for ( turbine_iterator tit = _turbines.begin(); tit != _turbines.end(); ++tit ) {
        if ( tit->second->level() < ladder.steps.size() ) {
            turbine_type& adjTurbine = *tit->second;
            EnergyTransform newETransform= ladder.steps[ adjTurbine.level() ];
            if ( adjTurbine.level() > 0
              && boost::math::isfinite( adjTurbine.energyTransform().minEnergyThreshold )
              && newETransform.minEnergyThreshold > adjTurbine.energyTransform().minEnergyThreshold
            ){
                newETransform.minEnergyThreshold = adjTurbine.energyTransform().minEnergyThreshold;
            }
            adjTurbine.setEnergyEval( energyEval );
            adjTurbine.setEnergyTransform( newETransform );
            LOG_INFO( UNIT_LOG_PREFIX( *this ) << "adjusted turbine #" << adjTurbine.index()
                          << " params: E_min=" << newETransform.minEnergyThreshold
                          << " T=" << newETransform.temperature );
        }
    }
    // create hotter turbines and connect unit to origin turbine
    const TurbineNode& originNode = _cascadeStructure[ originIx ];
    LOG_DEBUG2( UNIT_LOG_PREFIX( *this ) << "checking if new turbine should be created" );
    for ( turbine_ix_t turbineIx = 0; turbineIx < _cascadeStructure.size(); ++turbineIx ) {
        const TurbineNode& nodeProps = _cascadeStructure[ turbineIx ];
        if ( nodeProps.unitIx != index() ) continue;
        // create turbine, if it's in this unit
        if ( ( nodeProps.level == originNode.level + 1 )
          && _turbines.find( nodeProps.turbineIx ) == _turbines.end()
        ){
            turbine_type* pNewTurbine = appendTurbine( nodeProps, ladder.steps[ nodeProps.level ] );
            BOOST_ASSERT( pNewTurbine->index() == nodeProps.turbineIx );
            BOOST_ASSERT( _turbines.find( nodeProps.turbineIx ) != _turbines.end() );
            pNewTurbine->setMovingParticle( iniParticle );
            pNewTurbine->setEnergyEval( energyEval );
        }
    }
}

template<class ParticleCollector, class DynamicParticleFactory, class ParticleGenerator, class UnitCommunicator>
bool TurbineCascadeUnit<ParticleCollector, DynamicParticleFactory, ParticleGenerator, UnitCommunicator>::sampleParticle(
    turbine_ix_t            originIx,
    const particle_type&    particle
){
    if ( _pCollector ) {
        bool collectorSaturated = _pCollector->storeSample( _timer.elapsed(), originIx, particle,
                                                            _turbines.begin()->second->movingParticle().energyEval() );
        if ( collectorSaturated && !_shuttingDown ) {
            LOG_DEBUG1( UNIT_LOG_PREFIX( *this ) << "sample collector is saturated, initiating shutdown" );
            _shuttingDown = true;
            if ( _pUnitComm ) _pUnitComm->broadcastShutdown( originIx );
        }
        return ( collectorSaturated );
    }
    else if ( _pUnitComm ) {
        _pUnitComm->sendSample( originIx, particle );
        return ( false );
    }
    else {
        THROW_RUNTIME_ERROR( UNIT_LOG_PREFIX( *this ) <<  "cannot store sample: disconnected unit without collector" );
    }
}
