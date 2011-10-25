#pragma once

#include "../BasicTypedefs.h"

/**
 *  Odds ratio of two alternatives.
 */
struct OddsRatio {
    log_prob_t  llhRatio;   /** logarithm of likelihood odds ratio */
    log_prob_t  lppRatio;   /** logarithm of prior probability odds ratio */ 

    OddsRatio()
    : llhRatio( unset() ), lppRatio( unset() )
    {
    }

    OddsRatio( log_prob_t llhRatio, log_prob_t lppRatio )
    : llhRatio( llhRatio ), lppRatio( lppRatio )
    {
        BOOST_ASSERT( !is_unset( llhRatio ) );
        BOOST_ASSERT( !is_unset( lppRatio ) );
    }

    OddsRatio( log_prob_t  llh1, log_prob_t llh2, 
               log_prob_t priorLP1, log_prob_t priorLP2 )
    : llhRatio( llh1 - llh2 ), lppRatio( priorLP1 - priorLP2 )
    {
        if ( std::isinf( llh1 ) !=0 && llh1 == llh2 ) {
            llhRatio = 0;
        }
        BOOST_ASSERT( !std::isnan( llhRatio ) );
        if ( std::isinf( priorLP1 ) && priorLP1 == priorLP2 ) {
            lppRatio = 0;
        }
        BOOST_ASSERT( !std::isnan( lppRatio ) );
    }

    /**
     *  Probability to select the first alternative.
     *   a/(a+b)  = 1/(1+b/a)
     */
    prob_t prob() const {
        return ( 1.0/ ( 1 + __builtin_exp( -llhRatio - lppRatio ) ) );
    }

    /**
     *  Log of probability to select the first alternative.
     *  ln a/(a+b) = -ln ( 1 + b/a )
     */
    log_prob_t lnProb() const {
        return ( ln_odds_to_prob( llhRatio + lppRatio ) );
    }
};

/**
 *  Result of Gibbs sampling step.
 */
template<typename ValueType>
struct GibbsSample {
    ValueType   value;      /** value of variable after Gibbs sampling */
    OddsRatio   oddsRatio;  /** odds ratio of newValue vs oldValue */

    GibbsSample( const ValueType& value, log_prob_t llhRatio, log_prob_t lppRatio )
    : value( value )
    , oddsRatio( llhRatio, lppRatio )
    {}

    GibbsSample( const ValueType& value, const OddsRatio& oddsRatio )
    : value( value ), oddsRatio( oddsRatio )
    {}

    operator const ValueType&() const
    {
        return ( value );
    }
};
