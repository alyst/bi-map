#pragma once

#include "../BasicTypedefs.h"

#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/indexed_by.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>
#include <boost/multi_index/mem_fun.hpp>

#include "ExecutionMonitor.h"
#include "Particle.h"

template<class Particle>
class EnergyDisk {
public:
    typedef Particle particle_type;
    typedef CascadeParticle<particle_type> cascade_particle_type;
    typedef typename cascade_particle_type::energy_type energy_type;
    typedef std::vector<energy_type> energy_vector_type;
    typedef EnergyDisk<particle_type> disk_type;

private:
    typedef boost::multi_index::multi_index_container<
    cascade_particle_type,
    boost::multi_index::indexed_by<
        // by assigned serial number
        boost::multi_index::ordered_non_unique<boost::multi_index::const_mem_fun<
                cascade_particle_type, const energy_type, &cascade_particle_type::energy> >,
        // by entity contents
        boost::multi_index::sequenced<>
    > > particles_container_type;

    typedef typename particles_container_type::template nth_index<0>::type particles_energy_index;
    typedef typename particles_container_type::template nth_index<1>::type particles_sequential_index;
    typedef typename particles_energy_index::const_iterator const_particle_iterator;
    typedef typename particles_energy_index::const_reverse_iterator const_reverse_particle_iterator;
    typedef typename particles_energy_index::iterator particle_iterator;
    typedef typename particles_sequential_index::const_iterator const_particle_sequential_iterator;

    std::size_t                         _maxParticles;
    particles_container_type            _particles;
    TurbineCascadeExecutionMonitor*     _pMonitor;

    particles_energy_index& particlesEnergyIndex() {
        return ( _particles.template get<0>() );
    }

    particles_sequential_index& particlesSequentialIndex() {
        return ( _particles.template get<1>() );
    }

public:
    EnergyDisk(
        std::size_t maxParticles
    ) : _maxParticles( maxParticles ), _pMonitor( NULL )
    {
    }

    ~EnergyDisk()
    {
        for ( const_particle_iterator pIt = particlesEnergyIndex().begin();
              pIt != particlesEnergyIndex().end(); ++pIt ) {
            if ( _pMonitor ) {
                /// @todo: submit proper values, change signature
                _pMonitor->notifyParticlePopped( 0, pIt->energy(), 0, false );
            }
        }
    }

    const particles_energy_index& particlesEnergyIndex() const {
        return ( _particles.template get<0>() );
    }

    const particles_sequential_index& particlesSequentialIndex() const {
        return ( _particles.template get<1>() );
    }

    void fill( const disk_type& energyDisk, double time, size_t iteration )
    {
        // paste particles to new disk
        for ( const_particle_iterator pIt = energyDisk.particlesEnergyIndex().begin();
              pIt != energyDisk.particlesEnergyIndex().end(); ++pIt )
        {
            push( *pIt, pIt->serial, 0, 0 );
        }
    }

    std::size_t capacity() const {
        return ( _maxParticles );
    }

    std::size_t size() const {
        return ( _particles.size() );
    }

    bool isEmpty() const {
        return ( _particles.empty() );
    }

    bool isFull() const {
        return ( particlesEnergyIndex().size() >= capacity() );
    }

    std::pair<const_particle_iterator, const_particle_iterator> findParticlesRange( energy_type minEnergy, energy_type maxEnergy ) const;

    const cascade_particle_type* pickRandomParticle( const gsl_rng* rng, energy_type minEnergy, energy_type maxEnergy ) const;

    const cascade_particle_type* energyQuantileParticle( prob_t quantile ) const;

    energy_type energyRange( prob_t quantileMin = 0, prob_t quantileMax = 1 ) const;

    const cascade_particle_type* nthParticle( size_t rank ) const {
        if ( particlesEnergyIndex().empty() ) return ( NULL );
        const_particle_iterator pIt = particlesEnergyIndex().begin();
        std::advance( pIt, std::min( rank, size()-1 ) );
        return ( &*pIt );
    }

    const cascade_particle_type* particleNotFound() const {
        return ( NULL );
    }

    const cascade_particle_type* push( const cascade_particle_type& particle ) ;

    const cascade_particle_type* push( const particle_type& particle, particle_serial_t serial,
                                       double time, size_t iteration
    ){
        return ( push( cascade_particle_type( serial, time, iteration, particle ) ) );
    }

    template<class Iterator>
    bool push_many( const Iterator&  particlesBegin,
                    const Iterator&  particlesEnd,
                    double time,
                    size_t iteration
    ){
        for ( Iterator pit = particlesBegin; pit != particlesEnd; ++pit ) {
            /// @todo serial
            cascade_particle_type mp( 0, time, iteration, *pit );
            if ( push( mp ) == particleNotFound() ) {
                //return ( false );
            }
        }
        return ( true );
    }

    size_t particlesCount() const {
        return ( particlesEnergyIndex().size() );
    }

    energy_vector_type particleEnergies(
        energy_type minEnergy = -std::numeric_limits<energy_type>::infinity(),
        energy_type maxEnergy = std::numeric_limits<energy_type>::infinity() ) const ;

    template<class OutputStream>
    void print( OutputStream& out, size_t iteration, size_t diskId ) const
    {
        energy_vector_type  energies;
        energies.reserve( particlesCount() );

        for ( const_particle_iterator pit = particlesEnergyIndex().begin(); pit != particlesEnergyIndex().end(); ++pit ) {
            energies.push_back( pit->energy() );
        }
        for ( int i = 0; i < energies.size(); i++ ) {
            out << iteration << '\t' << diskId << '\t' << energies[ i ] << '\n';
        }
    }
};

/******************************************************************************
 * Implementation
 ******************************************************************************/

template<class Particle>
typename EnergyDisk<Particle>::energy_type EnergyDisk<Particle>::energyRange(
    prob_t  quantileMin,
    prob_t  quantileMax
) const {
    if ( quantileMin > quantileMax ) return ( std::numeric_limits<double>::signaling_NaN() );
    if ( particlesEnergyIndex().size() <= 1 ) return ( 0 );
    const cascade_particle_type* eMin = energyQuantileParticle( quantileMin );
    if ( eMin == NULL ) return ( std::numeric_limits<double>::signaling_NaN() );
    const cascade_particle_type* eMax = energyQuantileParticle( quantileMax );
    if ( eMax == NULL ) return ( std::numeric_limits<double>::signaling_NaN() );
    return ( eMax->energy() - eMin->energy() );
}

template<class Particle>
const typename EnergyDisk<Particle>::cascade_particle_type*
EnergyDisk<Particle>::energyQuantileParticle(
    prob_t  quantile
) const {
    if ( particlesEnergyIndex().empty() || quantile < 0 || quantile > 1 ) {
        return ( NULL );
    }
    size_t ix = (int)( std::floor( quantile * particlesEnergyIndex().size() ) );
    if ( ix < particlesEnergyIndex().size() / 2 ) {
        const_particle_iterator pIt = particlesEnergyIndex().begin();
        if ( ix > 0 ) std::advance( pIt, ix );
        return ( &*pIt );
    }
    else {
        const_reverse_particle_iterator pRit = particlesEnergyIndex().rbegin();
        if ( ix < particlesEnergyIndex().size() - 1 ) std::advance( pRit, particlesEnergyIndex().size() - 1 - ix );
        return ( &*pRit );
    }
}

template<class Particle>
const typename EnergyDisk<Particle>::cascade_particle_type*
EnergyDisk<Particle>::pickRandomParticle(
    const gsl_rng*  rng,
    energy_type     minEnergy,
    energy_type     maxEnergy
) const {
#if 0
    double rndEnergy = meanEnergy + gsl_ran_gaussian( rng, sigmaEnergy );
    const_particle_iterator pit = particlesEnergyIndex().upper_bound( rndEnergy );
    if ( pit == particlesEnergyIndex().begin() )        return ( NULL );
    pit--;
    return ( pit->energy <= rndEnergy ? &(*pit) : NULL );
#else
    std::pair<const_particle_iterator, const_particle_iterator> range = findParticlesRange( minEnergy, maxEnergy );
    if ( range.first == particlesEnergyIndex().end() ) return ( NULL );
    size_t rangeSize = 0;
    for ( const_particle_iterator pit = range.first; pit != range.second; ++pit ) {
        rangeSize++;
    }
    const_particle_iterator pit = range.first;
    size_t offset = rangeSize > 0 ? gsl_rng_uniform_int( rng, rangeSize ) : 0;
    if ( offset > 0 ) {
        std::advance( pit, offset );
    }
    return ( &*pit );
#endif
}

template<class Particle>
std::pair<typename EnergyDisk<Particle>::const_particle_iterator,
          typename EnergyDisk<Particle>::const_particle_iterator> 
EnergyDisk<Particle>::findParticlesRange(
    energy_type minEnergy,
    energy_type maxEnergy
) const {
    BOOST_ASSERT( minEnergy <= maxEnergy );
    return ( std::make_pair( particlesEnergyIndex().lower_bound( minEnergy ),
                             particlesEnergyIndex().upper_bound( maxEnergy ) ) );
}

template<class Particle>
const typename EnergyDisk<Particle>::cascade_particle_type*
EnergyDisk<Particle>::push(
    const cascade_particle_type& particle
){
    if ( !std::isfinite( particle.energy() ) ) {
        LOG_DEBUG2( "Particle energy E=" << particle.energy() << ": not finite, skipped" );
        return ( particleNotFound() );
    }
    if ( particle.energy() < 0 ) {
        LOG_DEBUG2( "Particle E=" << particle.energy() );
    }
    const cascade_particle_type* pushed = NULL;
    if ( isFull() ) {
#if 0
        if ( _pMonitor ) {
            const_particle_sequential_iterator oldPit = particlesSequentialIndex().begin();
            _pMonitor->notifyParticlePopped( _index, oldPit->energy, oldPit->outerIteration, oldPit->released );
        }
#endif
        // remove the earliest particle stored
        particlesSequentialIndex().pop_front();
    }
    std::pair<const_particle_iterator, bool> insRes = particlesEnergyIndex().insert( particle );
    if ( insRes.second ) {
        pushed = &*insRes.first;
#if 0
        if ( _pMonitor ) _pMonitor->notifyParticlePushed( _index, pushed->energy );
#endif
    } else {
        LOG_WARN( "Particle E=" << particle.energy() << " not pushed" );
    }
    return ( pushed );
}

/**
 *  Gets a vector of energies of all particles in the disk
 *  within the specified energy range.
 */
template<class Particle>
typename EnergyDisk<Particle>::energy_vector_type
EnergyDisk<Particle>::particleEnergies(
    energy_type minEnergy,   /** minimal energy */
    energy_type maxEnergy    /** maximal energy */
) const {
    energy_vector_type  res;
    res.reserve( particlesCount() );

    for ( const_particle_iterator pit = particlesEnergyIndex().begin();
          pit != particlesEnergyIndex().end(); ++pit ) {
        if ( pit->energy() > minEnergy && pit->energy() <= maxEnergy ) {
            res.push_back( pit->energy() );
        }
    }
    LOG_DEBUG1( res.size() << " particles E=(" << minEnergy
                << " (" << ( res.size() > 0 ? res[0] : std::numeric_limits<double>::quiet_NaN() )
                << "), " << maxEnergy
                << "])" );
    return ( res );
}

