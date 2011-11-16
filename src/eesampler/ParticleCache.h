#pragma once

#include "../BasicTypedefs.h"

#include <map>
#include <vector>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "ExecutionMonitor.h"
#include "Particle.h"

template<class Particle>
class ParticleCacheEnergiesProxy;

template<class Particle>
class ParticleCache {
public:
    typedef Particle particle_type;
    typedef CascadeParticle<particle_type> cascade_particle_type;
    typedef ParticleCache<particle_type> cache_type;
    typedef ParticleCacheEnergiesProxy<particle_type> energies_proxy_type;

private:
    typedef std::vector<cascade_particle_type> particles_container_type;

    typedef typename particles_container_type::const_iterator const_particle_iterator;
    typedef typename particles_container_type::iterator particle_iterator;

    std::size_t                         _maxParticles;
    particles_container_type            _particles;
    TurbineCascadeExecutionMonitor*     _pMonitor;

public:
    ParticleCache(
        std::size_t maxParticles
    ) : _maxParticles( maxParticles ), _pMonitor( NULL )
    {
    }

    ~ParticleCache()
    {
        for ( const_particle_iterator pIt = _particles.begin();
              pIt != _particles.end(); ++pIt ) {
            if ( _pMonitor ) {
                /// @todo: submit proper values, change signature
                _pMonitor->notifyParticlePopped( 0, 0, 0, false );
            }
        }
    }

    void fill( const cache_type& cache, double time, size_t iteration )
    {
        // paste particles to the new cache
        for ( size_t ix = 0; ix < cache.size(); ++ix )
        {
            push( cache[ ix ], cache[ ix ].serial, 0, 0 );
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
        return ( size() >= capacity() );
    }

    const cascade_particle_type& operator[]( size_t index ) const {
        return ( _particles[ index ] );
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
            if ( push( mp ) == NULL ) {
                //return ( false );
            }
        }
        return ( true );
    }

    size_t particlesCount() const {
        return ( size() );
    }

    const cascade_particle_type* particleNotFound() const {
        return ( NULL );
    }

    template<typename EnergyEval>
    energies_proxy_type energies( const EnergyEval& eval ) const
    {
        return ( ParticleCacheEnergiesProxy<Particle>( *this, eval ) );
    }
};

template<class Particle>
class ParticleCacheEnergiesProxy {
public:
    typedef Particle particle_type;
    typedef float energy_type;
    typedef CascadeParticle<particle_type> cascade_particle_type;
    typedef ParticleCache<particle_type> cache_type;
    typedef ParticleCacheEnergiesProxy<particle_type> energies_proxy_type;
    typedef std::vector<energy_type> energy_vector_type;

private:
    typedef size_t sequence_index;
    typedef std::multimap<energy_type, sequence_index> energy_to_sequence_map_type;

    friend class ParticleCache<Particle>;

    class ParticleEnergyIterator {
        typedef typename energy_to_sequence_map_type::const_iterator internal_iterator;

        const energies_proxy_type&   _energies;
        internal_iterator           _it;

    protected:
        ParticleEnergyIterator( const energies_proxy_type& energies )
        : _energies( energies ), _it( energies._energyMap.end() )
        {
        }

        ParticleEnergyIterator( const energies_proxy_type& energies, const internal_iterator& it )
        : _energies( energies ), _it( it )
        {
        }

    public:
        bool operator==( const ParticleEnergyIterator& that ) const
        {
            if ( &_energies != &that._energies ) throw std::runtime_error( "Iterator belong to different containers" );
            return ( _it == that._it );
        }
        bool operator!=( const ParticleEnergyIterator& that ) const
        {
            return ( !operator==( that ) );
        }
        const cascade_particle_type& operator*() const {
            return ( _energies._cache[ _it->second ] );
        }
        const cascade_particle_type& operator->() const {
            return ( _energies._cache[ _it->second ] );
        }
        ParticleEnergyIterator operator++() {
            ParticleEnergyIterator res( *this );
            _it++;
            return ( res );
        }
        ParticleEnergyIterator& operator++(int) {
            _it++;
            return ( *this );
        }
    };

    typedef ParticleEnergyIterator const_particle_iterator;
    typedef energy_to_sequence_map_type::const_iterator const_energy_iterator;
    typedef energy_to_sequence_map_type::const_reverse_iterator const_reverese_energy_iterator;

protected:
    const cache_type& _cache;
    energy_to_sequence_map_type _energyMap;

    template<typename EnergyEval>
    ParticleCacheEnergiesProxy( const ParticleCache<Particle>& cache, const EnergyEval& eval )
    : _cache( cache )
    {
        for ( size_t ix = 0; ix < cache.size(); ix++ ) {
            _energyMap.insert( std::make_pair( eval( cache[ ix ] ), ix ) );
        }
    };

    const_energy_iterator energyQuantile( prob_t quantile ) const;

public:
    static const sequence_index ParticleNA = (sequence_index)( -1 );

    energy_vector_type particleEnergies(
        energy_type minEnergy = -std::numeric_limits<energy_type>::infinity(),
        energy_type maxEnergy = std::numeric_limits<energy_type>::infinity() ) const;

    const cascade_particle_type* pickRandomParticle( const gsl_rng* rng, energy_type minEnergy, energy_type maxEnergy ) const;

    energy_type energyRange( prob_t quantileMin = 0, prob_t quantileMax = 1 ) const;

    const_particle_iterator begin() const {
        return ( const_particle_iterator( *this, _energyMap.begin() ) );
    }
    const_particle_iterator end() const {
        return ( const_particle_iterator( *this ) );
    }

    size_t size() const {
        return ( _energyMap.size() );
    }

    const cascade_particle_type* nthParticle( size_t rank ) const {
        if ( _energyMap.empty() ) return ( NULL );
        const_energy_iterator pIt = _energyMap.begin();
        std::advance( pIt, std::min( rank, _energyMap.size()-1 ) );
        return ( &_cache[ pIt->second ] );
    }
    const cascade_particle_type* particleNotFound() const {
        return ( NULL );
    }

    template<class OutputStream>
    void print( OutputStream& out, size_t iteration, size_t cacheId ) const
    {
        energy_vector_type  energies;
        energies.reserve( _energyMap.size() );

        for ( energy_to_sequence_map_type::const_iterator pit = _energyMap.begin(); pit != _energyMap.end(); ++pit ) {
            out << iteration << '\t' << cacheId << '\t' << pit->first << '\n';
        }
    }

};

/******************************************************************************
 * Implementation
 ******************************************************************************/

template<class Particle>
typename ParticleCacheEnergiesProxy<Particle>::energy_type ParticleCacheEnergiesProxy<Particle>::energyRange(
    prob_t  quantileMin,
    prob_t  quantileMax
) const {
    if ( quantileMin > quantileMax ) return ( std::numeric_limits<double>::signaling_NaN() );
    if ( _energyMap.size() <= 1 ) return ( 0 );
    const const_energy_iterator eMinIt = energyQuantile( quantileMin );
    if ( eMinIt == _energyMap.end() ) return ( std::numeric_limits<double>::signaling_NaN() );
    const const_energy_iterator eMaxIt = energyQuantile( quantileMax );
    if ( eMaxIt == _energyMap.end() ) return ( std::numeric_limits<double>::signaling_NaN() );
    return ( eMaxIt->first - eMinIt->first );
}

template<class Particle>
typename ParticleCacheEnergiesProxy<Particle>::const_energy_iterator
ParticleCacheEnergiesProxy<Particle>::energyQuantile(
    prob_t  quantile
) const {
    if ( _energyMap.empty() || quantile < 0 || quantile > 1 ) {
        return ( _energyMap.end() );
    }
    size_t ix = (int)( std::floor( quantile * _energyMap.size() ) );
    if ( ix < _energyMap.size() / 2 ) {
        const_energy_iterator eIt = _energyMap.begin();
        if ( ix > 0 ) std::advance( eIt, ix );
        return ( eIt );
    }
    else {
        energy_to_sequence_map_type::const_iterator eIt = _energyMap.end();
        for ( ; ix < _energyMap.size(); ix++ ) {
            eIt--;
        }
        return ( eIt );
    }
}

template<class Particle>
const typename ParticleCacheEnergiesProxy<Particle>::cascade_particle_type*
ParticleCacheEnergiesProxy<Particle>::pickRandomParticle(
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
    const_energy_iterator lowerBound = _energyMap.lower_bound( minEnergy );
    const_energy_iterator upperBound = _energyMap.upper_bound( maxEnergy );
    if ( lowerBound == _energyMap.end() ) return ( NULL );
    size_t rangeSize = 0;
    for ( const_energy_iterator pit = lowerBound; pit != upperBound; ++pit ) {
        rangeSize++;
    }
    const_energy_iterator pit = lowerBound;
    size_t offset = rangeSize > 0 ? gsl_rng_uniform_int( rng, rangeSize ) : 0;
    if ( offset > 0 ) {
        std::advance( pit, offset );
    }
    return ( &_cache[ pit->second ] );
#endif
}

template<class Particle>
const typename ParticleCache<Particle>::cascade_particle_type*
ParticleCache<Particle>::push(
    const cascade_particle_type& particle
){
    if ( !std::isfinite( particle.energy() ) ) {
        LOG_DEBUG2( "Particle energy E=" << particle.energy() << ": not finite, skipped" );
        return ( particleNotFound() );
    }
    if ( particle.energy() < 0 ) {
        LOG_DEBUG2( "Particle E=" << particle.energy() );
    }
    if ( isFull() ) {
#if 0
        if ( _pMonitor ) {
            const_particle_sequential_iterator oldPit = particlesSequentialIndex().begin();
            _pMonitor->notifyParticlePopped( _index, oldPit->energy, oldPit->outerIteration, oldPit->released );
        }
#endif
        // remove the earliest particle stored, all ParticleCacheEnergiesProxy<> should be discarded as indices changed
        _particles.erase( _particles.begin() );
    }
    _particles.push_back( particle );
    return ( &_particles.back() );
}

/**
 *  Gets a vector of energies of all particles in the cache
 *  within the specified energy range.
 */
template<class Particle>
typename ParticleCacheEnergiesProxy<Particle>::energy_vector_type
ParticleCacheEnergiesProxy<Particle>::particleEnergies(
    energy_type minEnergy,   /** minimal energy */
    energy_type maxEnergy    /** maximal energy */
) const {
    energy_vector_type  res;
    res.reserve( _energyMap.size() );

    const const_energy_iterator ubIt = _energyMap.upper_bound( maxEnergy );
    for ( const_energy_iterator eit = _energyMap.lower_bound( minEnergy );
          eit != ubIt; ++eit ) {
        res.push_back( eit->first );
    }
    LOG_DEBUG1( res.size() << " particles E=(" << minEnergy
                << " (" << ( res.size() > 0 ? res[0] : std::numeric_limits<double>::quiet_NaN() )
                << "), " << maxEnergy
                << "])" );
    return ( res );
}
