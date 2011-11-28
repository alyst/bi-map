#pragma once

#include "../BasicTypedefs.h"

#include <vector>
#include <gsl/gsl_randist.h>

#include <boost/optional.hpp>

#include "Particle.h"

/**
    Parameters of turbine's energy evaluation.
 */
struct EnergyTransform {
    energy_type         minEnergyThreshold;     /** the minimal energy the particle could have,
                                                    everything below is clamped, when doing sampling steps */
    energy_type         temperature;            /** temperature of the particles (in Bolzmann's sense) */

    EnergyTransform( energy_type minEnergyThreshold = -std::numeric_limits<energy_type>::infinity(), energy_type temperature = 1.0 )
    : minEnergyThreshold( minEnergyThreshold )
    , temperature( temperature )
    {}

    bool isNonTransforming() const
    {
        return ( ( temperature == 1.0 )
                 && ( minEnergyThreshold < 0 && std::isinf( minEnergyThreshold ) ) );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( minEnergyThreshold );
        ar & BOOST_SERIALIZATION_NVP( temperature );
    }
};

/**
 *  Parameters of equi-energy jump.
 */
struct EEJumpParams {
    EnergyTransform acceptorETransform;     /** energy transform of the turbine requesting EE jump */
    energy_type     energy;                 /** current particle energy of acceptor */

    EEJumpParams( EnergyTransform acceptorETransform, energy_type energy )
    : acceptorETransform( acceptorETransform ), energy( energy )
    {}

    EEJumpParams()
    : energy( std::numeric_limits<energy_type>::quiet_NaN() )
    {}

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( acceptorETransform );
        ar & BOOST_SERIALIZATION_NVP( energy );
    }
};

/**
 *  Dummy implementation of TurbineConnection concept.
 *  TurbineConnection is an adapter for communication between colder and hotter turbines.
 */
template<class ParticleEnergyEval>
struct DummyTurbineConnection {
    typedef ParticleEnergyEval energy_eval_type;
    typedef typename energy_eval_type::particle_type particle_type;

    /**
     *  Equi-energy jump request status.
     */
    struct DummyJumpRequestStatus {
        bool complete() const { return ( true ); }
        boost::optional<CascadeParticle<particle_type> > particle() const {
            return ( boost::optional<CascadeParticle<particle_type> >() );
        }
    };
    typedef DummyJumpRequestStatus jump_request_status_type;

    size_t level() const { return ( 0 ); }
    size_t prisonersCount() const { return ( 0 ); }
    jump_request_status_type requestJump( turbine_ix_t initiatorIx,
                                          const energy_eval_type& eval,
                                          const EEJumpParams& params ) {
        return ( jump_request_status_type() );
    }
    template<class ParticleIterator>
    bool sendParticles( turbine_ix_t originIx,
                        const ParticleIterator& particlesBegin,
                        const ParticleIterator& particlesEnd ) { return ( true ); }
    template<class CascadeUnit>
    bool processExternalEvents( CascadeUnit& unit, bool wait ) { return ( true ); }
};
