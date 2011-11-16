#pragma once

#include "../BasicTypedefs.h"

#include <map>

#include <boost/format.hpp>
#include <boost/optional.hpp>
#include <boost/make_shared.hpp>
#include <boost/scoped_ptr.hpp>

#include "EnergyDisk.h"
#include "ITurbineConnection.h"

bool is_equienergy_jump_accepted( const gsl_rng* rng,
                                  energy_type curEnergy, const EnergyTransform& curTransform,
                                  energy_type jumpEnergy, const EnergyTransform& jumpTransform );

#if 0
class IParticleEnergyEval {
public:
    virtual energy_type operator()( const IParticle& particle ) const = 0;
};

class IDynamicParticle {
public:
    virtual ~IDynamicParticle() {};

    virtual energy_type energy() const = 0;
    virtual operator const IParticle&() const = 0;

    virtual void iterate() = 0;
    virtual void setParams( energy_type minEnergyThreshold, energy_type temperature ) = 0;
    virtual IDynamicParticle& operator=( const IParticle& particle ) = 0;
};

class IDynamicParticleFactory {
public:
    virtual IDynamicParticle* operator()( energy_type minEnergyThreshold, energy_type temperature ) const = 0;
};

class IParticleGenerator {
public:
    typedef ParticleCache particle_cache_type;
    typedef std::vector<const IParticle*> particle_container_type;

    virtual particle_container_type operator()( const gsl_rng* rng, const particle_cache_type& particles ) const = 0;
};
#endif

/**
    Equi-energy sampler parameters.
 */
struct TurbineParams {
    // ee-jump params
    prob_t      eeJumpRate;             /// rate at which equi-energy jumps are attempted
    energy_type eeJumpEnergyTolerance;  /** fraction of current particles cache energy range,
                                            to be used for selecting candidates of particles
                                            for jump to */
    prob_t      eeMaxEnergyQuantile;    /** quantile of particles' energies distribution 
                                            to use for calculation of energy range for
                                            turbines parameters adjusting and 
                                            equi-energy jumps */
    prob_t      eeLadderEnergyQuantile; /** quantile of particles' energies distribution 
                                            to use for calculation of energy ladder
                                            step threshold */

    size_t      maxParticles;           /** max number of particles in turbine's cache */
    size_t      particleSnapshotPeriod; /** period (# of iterations) at which particles are
                                            save to the cache of the turbine */
    size_t      particleSamplingPeriod; /** period (# of iterations) at which current particle 
                                            is collected for the final output */

    // generator-related
    prob_t      generateRate;           /** rate at which generation is initiated */
    prob_t      prisonerOnWalkSwitchRate;   /** rate at which active particle in the jail is being switch */
    size_t      detentionIterations;    /** number of cycles in prison a particle should make before it's released to
                                            the cache */

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( eeJumpRate );
        ar & BOOST_SERIALIZATION_NVP( eeJumpEnergyTolerance );
        ar & BOOST_SERIALIZATION_NVP( eeMaxEnergyQuantile );
        ar & BOOST_SERIALIZATION_NVP( maxParticles );
        ar & BOOST_SERIALIZATION_NVP( particleSnapshotPeriod );
        ar & BOOST_SERIALIZATION_NVP( particleSamplingPeriod );
        ar & BOOST_SERIALIZATION_NVP( generateRate );
        ar & BOOST_SERIALIZATION_NVP( prisonerOnWalkSwitchRate );
        ar & BOOST_SERIALIZATION_NVP( detentionIterations );
    }

    TurbineParams();
};


/**
    Implementation of ParticleGenerator concept that generates nothing.
 */
template<class ParticleEnergyEval>
class SunkenParticleGenerator{
public:
    typedef ParticleEnergyEval static_particle_energy_eval_type;
    typedef typename static_particle_energy_eval_type::particle_type particle_type;
    typedef std::vector<particle_type> particle_container_type;
    typedef ParticleCache<particle_type> particle_cache_type;

    virtual particle_container_type operator()( const gsl_rng* rng,
                                                const static_particle_energy_eval_type& energyEval,
                                                const particle_cache_type& particles )
    {
        return ( particle_container_type() );
    }

    virtual ~SunkenParticleGenerator()
    {}
};

/**
 *  Particle turbine missing iteration logic and connection to other turbines.
 */
template<class DynamicParticleFactory>
class ParticleTurbineBase {
public:
    typedef DynamicParticleFactory factory_type;
    typedef typename DynamicParticleFactory::dynamic_particle_type dynamic_particle_type;
    typedef typename DynamicParticleFactory::static_particle_type particle_type;
    typedef typename dynamic_particle_type::static_particle_energy_eval_type static_particle_energy_eval_type;
    typedef boost::scoped_ptr<dynamic_particle_type> dynamic_particle_pointer_type;
    typedef ParticleCache<particle_type> cache_type;
    typedef ParticleCacheEnergiesProxy<particle_type> energies_proxy_type;
    typedef typename energies_proxy_type::energy_vector_type energy_vector_type;
    typedef std::map<size_t, energy_vector_type> energy_landscape_type;
    typedef typename cache_type::cascade_particle_type cascade_particle_type;
    typedef boost::shared_ptr<cache_type>  cache_pointer_type;

private:
    /**
        Particle detained in a jail.
        Generated particles are put to the jail of higher-temperature turbine
        to be sampled there and become statistically independent from their parents,
        before being released to the cache of the turbine.
     */
    struct DetainedParticle {
        particle_type   particle;
        size_t          iterationsDone;     /// # of sampling iterations made in detention

        DetainedParticle( const particle_type& particle, size_t iterationsDone = 0 )
        : particle( particle ), iterationsDone( iterationsDone )
        {
        }
    };

    const size_t                _level;                     /** step level in temperature ladder */
    const turbine_ix_t          _index;                     /** unique index of turbine */
    const TurbineParams         _params;
    EnergyTransform             _energyTransform;

    dynamic_particle_pointer_type   _movingParticle;        /** sampling particle */
    cache_pointer_type              _pParticleCache;        /** cache where particles are saved */

    static const size_t PRISONER_NA = (size_t)(-1);
    dynamic_particle_pointer_type   _prisonerOnWalk;        /** detained particle that is sampling currently */
    size_t                          _prisonerOnWalkIx;      /** index of sampling detained particle in _jail array */
    size_t                          _prisonerIterations;    /** how many sampling steps prisoner on walk already made */
    std::vector<DetainedParticle>   _jail;                  /** array of detained particles */

protected:
    const gsl_rng*              _rng;
    size_t                      _iteration;                 /** # of iterations made */

protected:
    ParticleTurbineBase( size_t level, size_t index, 
                     const TurbineParams& params,
                     const EnergyTransform& energyTransform,
                     const gsl_rng* rng,
                     const factory_type& dynParticleFactory,
                     cache_type* particleCache = NULL,
                     TurbineCascadeExecutionMonitor* pMonitor = NULL
    ) : _level( level ), _index( index )
      , _params( params ), _energyTransform( energyTransform )
      , _movingParticle( dynParticleFactory( _energyTransform.minEnergyThreshold, _energyTransform.temperature ) )
      , _prisonerOnWalk( dynParticleFactory( _energyTransform.minEnergyThreshold, _energyTransform.temperature ) )
      , _prisonerOnWalkIx( PRISONER_NA )
      , _prisonerIterations( 0 )
      , _rng( rng )
      , _iteration( 0 )
    {
        if ( _params.maxParticles > 0 ) {
            _pParticleCache.reset( new cache_type( _params.maxParticles ) );
            if ( particleCache ) _pParticleCache->fill( *particleCache, 0, 0 );
        }
    }

    void moveParticle() {
        _movingParticle->iterate();
    }

    const cascade_particle_type* manageJail();

public:
    size_t level() const {
        return ( _level );
    }

    turbine_ix_t index() const {
        return ( _index );
    }

    const TurbineParams& params() const {
        return ( _params );
    }

    const EnergyTransform& energyTransform() const {
        return ( _energyTransform );
    }

    const dynamic_particle_type& movingParticle() const {
        return ( *_movingParticle );
    }

    const dynamic_particle_type& setMovingParticle( const particle_type& particle ) {
        *_movingParticle = particle;
        return ( *_movingParticle );
    }

    void setEnergyTransform( const EnergyTransform& energyTransform ) {
        _energyTransform = energyTransform;
        _movingParticle->setParams( _energyTransform.minEnergyThreshold, _energyTransform.temperature );
        _prisonerOnWalk->setParams( _energyTransform.minEnergyThreshold, _energyTransform.temperature );
    }

    size_t iteration() const {
        return ( _iteration );
    }

    boost::optional<const cascade_particle_type> processJumpRequest(
        const EEJumpParams&     params,
        const static_particle_energy_eval_type& energyEval,
        bool unboundMinEnergy = false
    ) const;

    const cache_pointer_type& particleCache() const {
        return ( _pParticleCache );
    }

    void removeCache()
    {
        _pParticleCache.reset();
    }

    size_t prisonersCount() const {
        return ( _jail.size() );
    }

    void detainParticle( double time, const particle_type& particle ) {
        if ( _params.detentionIterations == 0 ) {
            // no jail required, put directly to the cache
            if ( _pParticleCache ) {
                /// @todo serial
                _pParticleCache->push( particle, (particle_serial_t)0, time, _iteration );
            }
            else {
                throw std::runtime_error("No cache to add particles to");
            }
        }
        else {
            // detain particles to make them statistically independent from original particles
            _jail.push_back( DetainedParticle( particle ) );
        }
    }
    template<class Container>
    void detainParticles( double time, const Container& particles ) {
        if ( _params.detentionIterations == 0 ) {
            // no jail required, put directly to the cache
            if ( _pParticleCache ) {
                _pParticleCache->push_many( particles, time, _iteration );
            }
            else {
                throw std::runtime_error("No cache to add particles to");
            }
        }
        else {
            // detain particles to make them statistically independent from original particles
            for ( typename Container::const_iterator it = particles.begin(); it != particles.end(); it++ ) {
                _jail.push_back( particle_type( *it ) );
            }
        }
    }
};

/**
 *  Particle turbine.
 *  Turbine does dynamic particle iteration (i.e. MCMC sampling),
 *  EE jumps via connection to hotter turbine
 *  and generation of particles to be sent to hotter turbine.
 */
template<class DynamicParticleFactory, class ParticleGenerator, class TurbineConnector>
class ParticleTurbine: public ParticleTurbineBase<DynamicParticleFactory> {
public:
    typedef ParticleTurbineBase<DynamicParticleFactory> super;
    typedef TurbineConnector turbine_connector_type;
    typedef ParticleGenerator generator_type;
    typedef typename TurbineConnector::jump_request_status_type jump_request_status_type;
    typedef typename super::factory_type factory_type;
    typedef typename super::cache_type cache_type;
    typedef typename super::cascade_particle_type cascade_particle_type;

private:
    using super::_rng;
    using super::_iteration;

    const size_t            _particleSnapshotOffset;
    const size_t            _particleSamplingOffset;
    const generator_type*   _pGenerator;                /** (optional) generator of particle to be sent to hotter turbines */
    turbine_connector_type  _conn;                      /** connection to other turbines of cascade */

protected:
    template<class CascadeUnit>
    bool processExternalEvents( CascadeUnit& unit, bool wait )
    {
        return ( _conn.processExternalEvents( unit, wait ) );
    }

    using super::moveParticle;
    using super::manageJail;

public:
    ParticleTurbine( size_t level, size_t index, 
                     const TurbineParams& params,
                     const EnergyTransform& energyTransform,
                     const gsl_rng* rng,
                     const factory_type& dynParticleFactory,
                     const generator_type* pGenerator,
                     const turbine_connector_type& connector,
                     size_t particleSnapshotOffset,
                     size_t particleSamplingOffset,
                     cache_type* particleCache = NULL,
                     TurbineCascadeExecutionMonitor* pMonitor = NULL
    ) : super( level, index, params, energyTransform, rng, dynParticleFactory,
               particleCache, pMonitor )
      , _particleSnapshotOffset( particleSnapshotOffset )
      , _particleSamplingOffset( particleSamplingOffset )
      , _pGenerator( pGenerator )
      , _conn( connector )
    {
    }

    using super::index;
    using super::params;
    using super::energyTransform;
    using super::movingParticle;
    using super::setMovingParticle;
    using super::particleCache;
    using super::iteration;

    size_t particleSamplingOffset() const {
        return ( _particleSnapshotOffset ); 
    }

    prob_t dynamicGenerationRate( turbine_ix_t destinationIx ) const {
        return ( params().generateRate > 0
                 ? params().generateRate / ( 1.0 + _conn.prisonersCount( destinationIx ) )
                 : 0.0 );
    }

    template<class CascadeUnit>
    void iterate( CascadeUnit& unit );
};

/******************************************************************************
 * ParticleTurbineBase Implementation
 ******************************************************************************/

/**
 *  Processes equi-energy jump request.
 */
template<class DynamicParticleFactory>
boost::optional<const typename ParticleTurbineBase<DynamicParticleFactory>::cascade_particle_type>
ParticleTurbineBase<DynamicParticleFactory>::processJumpRequest(
    const EEJumpParams&     params,         /** i parameters of the jump */
    const static_particle_energy_eval_type& energyEval,     /** i energy evaluator used to calculate energies for particles in the cache */
    bool                unboundMinEnergy    /** i if bound the low energy of the jump alternatives,
                                                  when true, acts like energy optimization
                                              */
) const {
    typedef boost::optional<const cascade_particle_type> optional_particle;

    const cascade_particle_type* pJumpParticle = NULL;
    if ( particleCache() ) {
        typename cache_type::energies_proxy_type energies = particleCache()->energies( energyEval );
        energy_type energyDelta = _params.eeJumpEnergyTolerance * energies.energyRange( 0, _params.eeMaxEnergyQuantile );
        pJumpParticle = energies.pickRandomParticle( _rng, unboundMinEnergy
                                                  ? -std::numeric_limits<energy_type>::infinity()
                                                  : params.energy - energyDelta,
                                                  params.energy + energyDelta );
    }
    if ( pJumpParticle == NULL ) {
        // no suitable particles
        return ( optional_particle() );
    }
    else {
        const energy_type donEnergy = pJumpParticle->energy();

        bool acceptJump = is_equienergy_jump_accepted( _rng, params.energy, params.acceptorETransform,
                                                        donEnergy, _energyTransform );
#if 0
        if ( _pMonitor && turbineIx < _params.turbinesCount ) {
                _pMonitor->notifyEquiEnergyJump( turbineIx, pJumpParticle->released,
                                                    request.energy, donEnergy, acceptJump );
        }
#endif
        return ( acceptJump ? optional_particle( *pJumpParticle ) : optional_particle() );
    }
}

template<class DynamicParticleFactory>
const typename ParticleTurbineBase<DynamicParticleFactory>::cascade_particle_type*
ParticleTurbineBase<DynamicParticleFactory>::manageJail()
{
    bool switchWalker = ( _jail.size() > 0 && _prisonerOnWalkIx == PRISONER_NA )
            || ( _jail.size() > 1 && ( gsl_rng_uniform( _rng ) < _params.prisonerOnWalkSwitchRate ) );
    if ( switchWalker ) {
        if ( _prisonerOnWalkIx != PRISONER_NA ) {
            LOG_DEBUG2( "Prisoner #" << _prisonerOnWalkIx << " goes back to jail" );
            // put walker back to jail
            _jail[ _prisonerOnWalkIx ].particle = *_prisonerOnWalk; 
            _jail[ _prisonerOnWalkIx ].iterationsDone = _prisonerIterations; 
        }
        // get new prisoner to walk
        _prisonerOnWalkIx = gsl_rng_uniform_int( _rng, _jail.size() );
        LOG_DEBUG2( "Prisoner #" << _prisonerOnWalkIx << " starts walking" );
        *_prisonerOnWalk = _jail[ _prisonerOnWalkIx ].particle;
        _prisonerIterations = _jail[ _prisonerOnWalkIx ].iterationsDone;
    }
    const cascade_particle_type* pIt = NULL;
    if ( _prisonerOnWalkIx != PRISONER_NA ) {
        LOG_DEBUG2( "Prisoner #" << _prisonerOnWalkIx << " on walk..." );
        _prisonerOnWalk->iterate();
        _prisonerIterations++;
        if ( _prisonerIterations >= _params.detentionIterations ) {
            // you are free, go to the cache
            pIt = _pParticleCache->push( *_prisonerOnWalk, 0, 0, _iteration );
            _prisonerIterations = 0;
            _jail.erase( _jail.begin() + _prisonerOnWalkIx );
            _prisonerOnWalkIx = PRISONER_NA;
        }
    }
    return ( pIt );
}

/******************************************************************************
 * ParticleTurbine Implementation
 ******************************************************************************/

template<class DynamicParticleFactory, class ParticleGenerator, class TurbineConnection>
template<class CascadeUnit>
void ParticleTurbine<DynamicParticleFactory, ParticleGenerator, TurbineConnection>::iterate(
    CascadeUnit&    unit
){
    // decide if equi-energy jump is required
    if ( _iteration % TURBINE_LOG_REPORT_PERIOD == 0 ) {
        if ( particleCache() && !particleCache()->isEmpty() ) {
            const cascade_particle_type& cpt = (*particleCache())[ particleCache()->size() - 1 ];
            // output last stored particle information
            LOG_INFO( "Turbine #" << super::index() << ": " << _iteration << " iteration(s) made, last particle E=" 
                        << cpt.energy() << ":\n" << cpt.particle );
        }
        else {
            LOG_INFO( "Turbine #" << super::index() << ": " << _iteration << " iteration(s) made, no particles cache" );
        }
    }

    turbine_ix_t jumpTurbineIx = gsl_rng_uniform( _rng ) <= params().eeJumpRate
                               ? _conn.randomEeJumpSource() : TURBINE_NA;
    bool  eeJumpDone = false;
    if ( jumpTurbineIx != TURBINE_NA ) {
        LOG_DEBUG1( "Turbine #" << super::index() << ": requesting jump from "
                    "#" << jumpTurbineIx << ", "
                    "E=" << boost::format("%.2f") % movingParticle().energy() );
        // try to equi-energy jump
        jump_request_status_type status = _conn.requestJump( index(), jumpTurbineIx, 
                                                             EEJumpParams( energyTransform(), movingParticle().energy() ) );
        while ( !status.complete() ) {
            bool abortIteration = processExternalEvents( unit, true ); // process any incoming external events
            if ( abortIteration ) return;
        }
        if ( status.particle() ) {
            LOG_DEBUG1( "Turbine #" << super::index() << ": "
                        "EE jumping to particle "
                        "E=" << boost::format( "%.2f" ) % status.particle()->energy() );
            setMovingParticle( *status.particle() );
            eeJumpDone = true;
        } else {
            LOG_DEBUG1( "Turbine #" << super::index() << ": EE jump failed" );
        }
    }
    // if not jumping, iterate normally
    if ( !eeJumpDone ) {
        moveParticle();
    }
    processExternalEvents( unit, false ); // process any incoming external events
    // optionally generate particles for higher-energy turbines
    turbine_ix_t destinationIx =  particleCache()
                               ? _conn.randomParticleDestination()
                               : TURBINE_NA;
    if ( destinationIx != TURBINE_NA ) {
        prob_t generationRate = dynamicGenerationRate( destinationIx );
        if ( generationRate > 0 && gsl_rng_uniform( _rng ) <= generationRate ) {
            const typename generator_type::particle_container_type generated =
                (*_pGenerator)( _rng, particleCache()->energies( movingParticle().staticParticleEnergyEval() ) );
#if 0
            if ( _pMonitor && turbineIx < _params.turbinesCount ) {
                _pMonitor->notifyParticlesGenerated( index, generated.size() );
            }
#endif
            _conn.sendParticles( index(), destinationIx, generated.begin(), generated.end() );
        }
    }
    processExternalEvents( unit, false ); // process any incoming external events
    // process particles in the 'jail'
    //const cascade_particle_type* pIt = 
    manageJail();
#if 0
    if ( pIt != NULL && _pMonitor ) {
        _pMonitor->notifyParticleReleased( index, pIt->energy );
    }
#endif
    processExternalEvents( unit, false ); // process any incoming external events
    // put particle to the cache periodically
    if ( particleCache()
         && ( ( _iteration + _particleSnapshotOffset ) % params().particleSnapshotPeriod == 0)
    ){
        particleCache()->push( movingParticle(), 0, 0, _iteration );
    }
    _iteration++;
}
