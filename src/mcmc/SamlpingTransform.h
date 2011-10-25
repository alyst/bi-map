#pragma once

#include "../BasicTypedefs.h"

/// floating-point type for temperature
typedef log_prob_t temperature_t;

struct SamplingTransform {
    log_prob_t      maxLnP;     /** ln of maximal probability */
    temperature_t   temperature;/** temperature (of Boltzmann distribution */

    SamplingTransform( log_prob_t maxLnP = std::numeric_limits<log_prob_t>::infinity(),
                       temperature_t temperature = 1.0 )
    : maxLnP( maxLnP ), temperature( temperature )
    {}

    /**
        Model Log ratio for MH step in case maximal probability is clamped.
    */
    static log_prob_t ClampedDeltaLnP( 
        log_prob_t deltaLnP,    /** i  difference between current and new unclamped probabilities */
        log_prob_t currentLnP,  /** i  logarithm of current unclamped probability */
        log_prob_t maxLnP       /** i  clamping threshold */
    ){
        if ( deltaLnP >= 0 ) {
            // both probabilities are clamped to maxLnP
            if ( currentLnP >= maxLnP )                 return ( 0 );
            // the new one would be clamped
            if ( currentLnP + deltaLnP >= maxLnP )      return ( maxLnP - currentLnP );
            // none clamped
                                                        return ( deltaLnP );
        }
        else {
            // both clamped
            if ( currentLnP + deltaLnP >= maxLnP )      return ( 0 );
            // the old one is clamped, correct the ratio
            if ( currentLnP >= maxLnP )                 return ( deltaLnP + ( currentLnP - maxLnP ) );
            // none clamped
                                                        return ( deltaLnP );
        }
    }

    log_prob_t operator()( log_prob_t lnP ) const {
        return ( std::min( lnP, maxLnP ) / temperature );
    }

    log_prob_t clampedDeltaLnP(
        log_prob_t deltaLnP,    /** i  difference between current and new unclamped probabilities */
        log_prob_t currentLnP   /** i  logarithm of current unclamped probability */
    ) const {
        return ( ClampedDeltaLnP( deltaLnP, currentLnP, maxLnP ) );
    }
};
