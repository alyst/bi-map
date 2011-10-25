#pragma once

#include "../BasicTypedefs.h"

#include <vector>
#include <boost/utility/enable_if.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/export.hpp>

typedef size_t turbine_ix_t;

#define TURBINE_NA  ((turbine_ix_t)(-1))

typedef  size_t particle_serial_t;

typedef float energy_type;

enum ParticleEventType {
    Sampled,
    Copied,
    Generated,
    JumpedTo
};

struct ParticleEvent {
    ParticleEventType   type;
    double              time;
    energy_type         energy;
    turbine_ix_t        tubineIndex;
    size_t              iteration;
    particle_serial_t   precursor1_serial;
    particle_serial_t   precursor2_serial;
};

typedef std::vector<ParticleEvent> ParticleHistory;

#if 0
/**
 *  Particle interface.
 *  Object, that you want to sample, must implement IParticle interface.
 */
class IParticle {
public:
    /** Clone object, i.e. return deep copy */
    virtual IParticle* clone() const = 0;
    virtual ~IParticle() {};

    virtual void serialize(boost::mpi::packed_iarchive& ar, const unsigned int version) = 0;
    virtual void serialize(boost::mpi::packed_oarchive& ar, const unsigned int version) = 0;
};
#endif

/**
 *  Particle + EE sampler internal information.
 */
template<class Particle>
struct CascadeParticle {
    typedef Particle particle_type;
    typedef ::energy_type energy_type;
    typedef CascadeParticle<particle_type> cascade_particle_type;

    /// sampler-unique serial # of particle
    particle_serial_t                       serial;
    /// particle's energy
    particle_type                           particle;

    CascadeParticle()
    : serial( 0 )
    {}

    CascadeParticle( particle_serial_t serial,
                     double time, size_t iteration, const particle_type& particle )
    : serial( serial ), particle( particle )
    {
    }

    friend bool operator<( const cascade_particle_type& a, const cascade_particle_type& b )
    {
        return ( a.particle.energy() < b.particle.energy() );
    }

    operator const particle_type&() const {
        return ( particle );
    }

    const energy_type energy() const {
        return ( particle.energy() );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( serial );
        ar & BOOST_SERIALIZATION_NVP( particle );
    }
};

template<class Particle>
typename boost::enable_if<boost::is_same<Particle, boost::shared_ptr< typename Particle::value_type > >, std::ostream>::type&
operator<<( std::ostream& out, const CascadeParticle<Particle>& p )
{
    return ( out << ( *(const Particle&)p ) );
}

template<class Particle>
typename boost::disable_if<boost::is_same<Particle, boost::shared_ptr< typename Particle::value_type > >, std::ostream>::type&
operator<<( std::ostream& out, const CascadeParticle<Particle>& p )
{
    return ( out << ( (const Particle&)p ) );
}
