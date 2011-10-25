#include <vector>

#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>

#include "logmath.h"

ln_cache_t __ln_factorial_cache;

ln_cache_t __log_int_cache;

#define SERIES_EPS 1E-5

/**
    Expansion at z -> 0.
    @see http://functions.wolfram.com/GammaBetaErf/GammaRegularized/06/01/05/01/01/
 */
log_prob_t lngammareg_inc_z0_series( log_prob_t a, log_prob_t z )
{
    log_prob_t res = 0.0;
    log_prob_t term = 1.0;
    size_t k = 0;
    while ( fabs( term ) > SERIES_EPS ) {
        k++;
        term *= ( -z ) / k;
        res += a / ( a + k ) * term;
    }
    return ( res );
}

/**
    Expansion at z -> infty.
    @see http://functions.wolfram.com/GammaBetaErf/GammaRegularized/06/02/02/
 */
log_prob_t lngammareg_inc_zinf_series( log_prob_t a, log_prob_t z )
{
    log_prob_t res = 0.0;
    log_prob_t term = 1.0;
    size_t k = 0;
    while ( fabs( term ) > SERIES_EPS ) {
        k++;
        term *= -( k - a ) / z;
        res += term;
    }
    return ( res );
}

/**
    Expansion at a -> infty.
    @see http://functions.wolfram.com/GammaBetaErf/GammaRegularized/06/02/01/
 */
inline log_prob_t lngammareg_inc_ainf_expansion( log_prob_t a, log_prob_t z )
{
    return ( ( 12 * z - 1 ) / ( 12 * a ) 
             + ( 288 * z * z - 312 * z + 1 ) / ( 288 * a * a ) );
}

log_prob_t ln_gammareg_inc_P( int a, log_prob_t z, log_prob_t lnZ )
{
    if ( lnZ < 0 ) {
        /* expand at z-> 0 */
        return ( a * lnZ - ln_factorial( a ) + fast_log1p( lngammareg_inc_z0_series( a, z ) ) );
    }
    else if ( a > 100 && a > 3 * z ) {
       ainf:
        /* expand at a-> inf */
        return ( ( -a - 0.5 ) * log( a ) + ( a - z ) + a * lnZ - 0.5 * ( M_LN2 + M_LNPI ) 
                 + gsl_log1p( lngammareg_inc_ainf_expansion( a, z ) ) );
    }
    else {
        log_prob_t res = gsl_sf_gamma_inc_P( a, z );
        if ( res == 0 ) goto ainf;
        return ( log( res ) );
    }
}

log_prob_t ln_gammareg_inc_Q( int a, log_prob_t z, log_prob_t lnZ )
{
    if ( a == 0 ) {
        return ( z - lnZ );
    }
    else if ( z > 10 * a ) {
        /* expand at z-> infty */
        return ( -z + ( a - 1 ) * lnZ - ln_factorial( a - 1 ) 
                 + fast_log1p( lngammareg_inc_zinf_series( a, z ) ) );
    }
    else {
        return ( log( gsl_sf_gamma_inc_Q( a, z ) ) );
    }
}
