#pragma once

#include "../BasicTypedefs.h"

#include <math.h>
#include <gsl/gsl_rng.h>

#include "../math/logmath.h"

#include "GibbsSample.h"
#include "SamlpingTransform.h"

/**
    Use Metropolis-Hastings rule to accept new sample
    of the modelled distribution.
 */
inline bool is_accept_step_by_mh(
    const gsl_rng*  r,          /**< @param[in] random number generator */
    prob_t          new_factor, /**< @param[in] probability x tran.prob(old->new) of the new sample */
    prob_t          cur_factor  /**< @param[in] probability x tran.prob(new->old) of the old sample */
){
    return ( ( new_factor >= cur_factor ) 
             || ( gsl_rng_uniform( r ) <= ( new_factor / cur_factor ) ) );
}

/**
 *  Calculate odds ratio for MH step.
 */
template<class ValueType, class DataLogLikelihood, class PriorLogProbability>
OddsRatio MetropolisHastringsPosteriorOddsRatio(
    const DataLogLikelihood&        llh,    /** likelihood calculator */
    const PriorLogProbability&      lpp,    /** prior probability calculator */
    const ValueType&                value,  /** current value */
    const ValueType&                newValue    /** output: new value */
){
    return ( OddsRatio( llh( newValue ), llh( value ), 
                        lpp( newValue ), lpp( value ) ) );
}

template<class ValueType, class TransitionDistribution, 
         class DataLogLikelihood, class PriorLogProbability>
ValueType MetropolisHastringsPosteriorSample(
    const gsl_rng*                  rng,
    const TransitionDistribution&   transDistr,
    const DataLogLikelihood&        dataLLH,
    const PriorLogProbability&      priorLP,
    const ValueType&                value
){
    ValueType   newValue = transDistr.generate( rng, value );
    if ( newValue == value )    return ( value );

    log_prob_t cur2newTransLP = transDistr.lnPdf( value, newValue );
    log_prob_t new2curTransLP = transDistr.lnPdf( newValue, value );
    OddsRatio logRatios = MetropolisHastringsPosteriorOddsRatio( dataLLH, priorLP, value, newValue );

    return ( is_accept_step_by_mh( rng, exp( logRatios.llhRatio + logRatios.lppRatio
                                             + new2curTransLP - cur2newTransLP ), 1 )
             ? newValue : value );
}

/// temperature-aware version
template<typename ValueType, class TransitionDistribution, 
         class DataLogLikelihood, class PriorLogProbability>
GibbsSample<ValueType> MetropolisHastringsPosteriorSample(
    const gsl_rng*                  rng,
    const TransitionDistribution&   transDistr,
    const DataLogLikelihood&        dataLLH,
    const PriorLogProbability&      priorLP,
    const ValueType&                value,
    const log_prob_t                totalLnP,
    const SamplingTransform&        transform
){
    ValueType   newValue = transDistr.generate( rng, value );
    if ( newValue == value )    return ( GibbsSample<ValueType>( value, 0, 0 ) );

    log_prob_t cur2newTransLP = transDistr.lnPdf( value, newValue );
    log_prob_t new2curTransLP = transDistr.lnPdf( newValue, value );

    const OddsRatio oddsRatio = MetropolisHastringsPosteriorOddsRatio( dataLLH, priorLP, value, newValue );

    // actual log prob ratio to use for sampling
    log_prob_t clampedLogRatio = transform.clampedDeltaLnP( oddsRatio.llhRatio + oddsRatio.lppRatio, totalLnP );
    BOOST_ASSERT( !std::isnan( clampedLogRatio ) );

    log_prob_t fullLogRatio = ( clampedLogRatio + new2curTransLP - cur2newTransLP ) / transform.temperature;
    if ( ( fullLogRatio > 0 ) || ( is_finite( fullLogRatio )
           && is_accept_step_by_mh( rng, __builtin_exp( fullLogRatio ), 1.0 ) )
    ){
        LOG_DEBUG3( "Accepted MH sample: logRatio=" << fullLogRatio
                    << " (llh=" << oddsRatio.llhRatio
                    << " ,lpp=" << oddsRatio.lppRatio << ")"
                    << " T=" << transform.temperature
                    << " totalLnP=" << totalLnP
                    << " maxLnP=" << transform.maxLnP );
        BOOST_ASSERT( !std::isinf( clampedLogRatio ) || clampedLogRatio > 0 );
        //BOOST_ASSERT( gsl_finite( logRatio.llh ) );
        //BOOST_ASSERT( gsl_finite( logRatio.priorLP ) );
        // modify the total probability
        return ( GibbsSample<ValueType>( newValue, oddsRatio ) );
    }
    else {
        return ( GibbsSample<ValueType>( value, 0, 0 ) );
    }
}
