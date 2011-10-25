/**
    Various transition distributions for implementing MCMC samplers.
 */

#pragma once

#include "../BasicTypedefs.h"
#include "../math/Distributions.h"

#include <utility>

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

struct BinaryTransition
{
    bool generate( const gsl_rng* rng, bool value ) const
    {
        return ( !value );
    }

    double lnPdf( bool oldValue, bool newValue ) const
    {
        return ( oldValue != newValue ? 0.0 : -std::numeric_limits<double>::infinity() );
    }
};

class GaussianTransitionDistribution
{
    signal_t  sigma;

public:
    GaussianTransitionDistribution( signal_t sigma )
    : sigma( sigma )
    {
    }

    signal_t generate( const gsl_rng* rng, signal_t value ) const
    {
        return ( value + gsl_ran_gaussian( rng, sigma ) );
    }

    log_prob_t lnPdf( signal_t oldValue, signal_t newValue ) const
    {
        return ( 0.0 );     // since transition is symmetric, no pdf is required - nominator and denominator terms are cancelled out
        //return ( gsl_ran_gaussian( newValue - oldValue, Sigma ) );
    }
};

class BivariateGaussianTransitionDistribution
{
public:
    typedef std::pair<signal_t, signal_t> value_type;

private:
    value_type sigma;
    signal_t ratio;

public:
    BivariateGaussianTransitionDistribution( signal_t sigma1, signal_t sigma2, signal_t ratio )
    : sigma( value_type( sigma1, sigma2 ) )
    , ratio( ratio )
    {
    }

    value_type generate( const gsl_rng* rng, const value_type& value ) const
    {
        double offs1, offs2;
        gsl_ran_bivariate_gaussian( rng, sigma.first, sigma.second, ratio, &offs1, &offs2 );
        return ( value_type( value.first + offs1, value.second + offs2 ) );
    }

    log_prob_t lnPdf( const value_type& oldValue, const value_type& newValue ) const
    {
        return ( 0.0 ); // symmetric
    }
};

class DiscreteGaussianTransitionDistribution
{
    signal_t  sigma;
    signal_t  step;

public:
    DiscreteGaussianTransitionDistribution( signal_t sigma, signal_t step = 0.01 )
    : sigma( sigma ), step( step )
    {
    }

    signal_t generate( const gsl_rng* rng, signal_t value ) const
    {
        return ( floor( ( value + gsl_ran_gaussian( rng, sigma ) ) / step ) * step );
    }

    log_prob_t lnPdf( signal_t oldValue, signal_t newValue ) const
    {
        return ( 0.0 );     // since transition is symmetric, no pdf is required - nominator and denominator terms are cancelled out
        //return ( gsl_ran_gaussian( newValue - oldValue, Sigma ) );
    }
};

class PositiveGaussianTransitionDistribution
{
    signal_t  sigma;

public:
    PositiveGaussianTransitionDistribution( signal_t sigma )
    : sigma( sigma )
    {
    }

    signal_t generate( const gsl_rng* rng, signal_t value ) const
    {
        signal_t newVal;
        do {
            newVal = value + gsl_ran_gaussian( rng, sigma );
        } while ( newVal < 0 );
        return ( newVal );
    }

    log_prob_t lnPdf( signal_t oldValue, signal_t newValue ) const
    {
        return ( gaussian_pdf_ln( newValue - oldValue, sigma )
                 - log( gsl_cdf_gaussian_Q( -oldValue, sigma ) ) );
    }
};

struct ExponentialTransitionDistribution
{
    signal_t generate( const gsl_rng* rng, signal_t value ) const
    {
        return ( gsl_ran_exponential( rng, M_LN2 / value ) );
    }

    log_prob_t lnPdf( signal_t oldValue, signal_t newValue ) const
    {
        return ( log( gsl_ran_exponential_pdf( newValue, M_LN2 / oldValue ) ) );
    }
};

class GammaTransitionDistribution
{
private:
    signal_t scale;

public:
    GammaTransitionDistribution( signal_t scale )
    : scale( scale )
    {}

    signal_t generate( const gsl_rng* rng, signal_t value ) const
    {
        double logVal = log( value );
        return ( exp( logVal + gsl_ran_gaussian( rng, scale ) ) );
    }

    log_prob_t lnPdf( signal_t oldValue, signal_t newValue ) const
    {
        return ( gaussian_pdf_ln( log( newValue ) - log( oldValue ), scale ) );
    }
};

class UniformTransitionDistribution
{
    signal_t  a;
    signal_t  b;

public:
     UniformTransitionDistribution( signal_t a, signal_t b )
     : a( a ), b( b )
     {}

    signal_t generate( const gsl_rng* rng, signal_t value ) const
    {
        return ( gsl_ran_flat( rng, a, b ) );
    }

    log_prob_t lnPdf( signal_t oldValue, signal_t newValue ) const
    {
        return ( log( gsl_ran_flat_pdf( newValue, a, b ) ) );
    }
};


#if 0
class ExponentialTransitionDistribution
{
    double rate;
public:
    ExponentialTransitionDistribution( double rate )
    : rate( rate )
    {}

    double generate( const gsl_rng* rng, double value ) const
    {
        return ( gsl_ran_exponential( rng, rate ) );
    }

    double pdf( double oldValue, double newValue ) const
    {
        return ( gsl_ran_exponential_pdf( newValue, rate ) );
    }
};
#endif
