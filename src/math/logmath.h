#pragma once

#include "../BasicTypedefs.h"

#include <cmath>

#include <vector>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>

typedef std::vector<log_prob_t> ln_cache_t;

#define MAX_CACHED_FACTORIAL 100000

extern ln_cache_t __ln_factorial_cache;
extern ln_cache_t __log_int_cache;

/**
 *  Memoize ln(n!) for n upto maxn.
 */
inline void ln_factorial_cache_fill(unsigned int maxn)
{
    size_t oldSize = __ln_factorial_cache.size();
    if ( oldSize <= maxn ) {
        __ln_factorial_cache.resize( maxn + 1, std::numeric_limits<log_prob_t>::quiet_NaN() );
    }
    for ( unsigned int i = oldSize; i <= maxn; i++ ) {
        log_prob_t& val = __ln_factorial_cache[ i ];
        val = gsl_sf_lnfact( i );
    }
}

/**
 *  Lookup into ln(n!) memoize table without any checking.
 */
inline log_prob_t fast_ln_factorial( unsigned int n )
{
    BOOST_ASSERT( __ln_factorial_cache.size() > n 
                  && !std::isnan(__ln_factorial_cache[ n ] ) );
    return ( __ln_factorial_cache[ n ] );
}

inline log_prob_t log_int( unsigned int n )
{
    if ( n >= __log_int_cache.size() ) {
        unsigned int i = __log_int_cache.size();
        __log_int_cache.resize( n + 1, std::numeric_limits<log_prob_t>::quiet_NaN() );
        for ( ; i < __log_int_cache.size(); i++ ) {
            __log_int_cache[ i ] = log( i );
        }
    }
    return ( __log_int_cache[ n ] );
}

inline log_prob_t ln_factorial(unsigned int n)
{
    if ( n >= MAX_CACHED_FACTORIAL ) {
        return ( gsl_sf_lnfact( n ) );
    }
    else if ( __ln_factorial_cache.size() <= n ) {
        ln_factorial_cache_fill( n );
    }
    return ( __ln_factorial_cache[ n ] );
}

inline log_prob_t fast_log1p( log_prob_t x )
{
    volatile log_prob_t y, z;
    y = 1 + x;
    z = y - 1;
    return ( __builtin_log( y ) - (z-x) / y );
}

/**
    log(1+exp(x))
 */
inline log_prob_t ln_1pexp( log_prob_t x )
{
    return ( x > -40 ? fast_log1p( __builtin_exp( x ) ) : 0 );
}

/**
    Converts log of odds ratio into probability:
    log(a/(a+b)) = log(1/(1+b/a)) = -log(1+exp(-ln a/b)) = log(a/b / (a/b + 1)) = log(a/b) - log(1+a/b)
 */
inline log_prob_t ln_odds_to_prob(
    log_prob_t odds_ratio   /** i ln a/b */
){
    return ( odds_ratio < 0 
             ? ( odds_ratio < -40 ? odds_ratio : -ln_1pexp( -odds_ratio ) )
             : ( odds_ratio > 40 ? 0 : odds_ratio - ln_1pexp( odds_ratio ) ) );
}

/**
    log(1-exp(x))
 */
inline log_prob_t ln_1mexp( log_prob_t x )
{
    if ( x > 0 ) {
        return ( std::numeric_limits<log_prob_t>::signaling_NaN() );
    }
    else if ( x == 0 ) {
        return ( -std::numeric_limits<log_prob_t>::infinity() );
    }
    else if ( x > -1E-2 ) {
        log_prob_t xx = x * x;
        return ( __builtin_log( -x ) + 0.5 * x + xx / 24 - ( xx * xx ) / 2880 );
    }
    else {
        return ( fast_log1p( -__builtin_exp( x ) ) );
    }
}

/**
 *  Log(1-x).
 */
inline log_prob_t ln_1m( prob_t x )
{
    return ( x > 0.5 ? __builtin_log( 1.0 - x ) : fast_log1p( -x ) );
}

/**
    Logarithm of the sum of a and b.
 */
inline log_prob_t logsumexp(
    log_prob_t  loga,   /**< @param[in] log(a) */
    log_prob_t  logb    /**< @param[in] log(b) */
){
    if ( -std::numeric_limits<log_prob_t>::infinity() == loga ) {
        return ( logb );
    } else if ( -std::numeric_limits<log_prob_t>::infinity() == logb ) {
        return ( loga );
    }
    else {
        return ( loga > logb
             ? loga + ln_1pexp( logb - loga )
             : logb + ln_1pexp( loga - logb ) );
    }
}

/**
    Logarithm of the difference of a and b.
 */
inline log_prob_t ln_minusexp(
    log_prob_t  loga,   /**< @param[in] log(a) */
    log_prob_t  logb    /**< @param[in] log(b) */
){
    if ( loga < logb ) throw std::invalid_argument( "loga < logb" );
    return ( loga == logb
             ? -std::numeric_limits<log_prob_t>::infinity()
             : loga + ln_1pexp( logb - loga ) );
}

log_prob_t ln_gammareg_inc_P( int a, log_prob_t z, log_prob_t lnZ );
log_prob_t ln_gammareg_inc_Q( int a, log_prob_t z, log_prob_t lnZ );

/**
 *  Logarithm of binomial coefficient.
 *  \binom( n, k ) = n!/k!/(n-k)!
 */
inline log_prob_t ln_binomial_coef( size_t n, size_t k )
{
    BOOST_ASSERT( k <= n );
    ln_factorial_cache_fill( n );
    return ( fast_ln_factorial( n ) - fast_ln_factorial( k ) - fast_ln_factorial( n - k ) );
}

#if 0
/**
 *  Fast exp() evaluation, based on
 *  "A Fast, Compact Approximation of the Exponential Function",
 *  Nicol N. Schraudolph
 *  @see http://cnl.salk.edu/~schraudo/pubs/Schraudolph99.pdf
 */
inline log_prob_t fast_exp( log_prob_t x )
{
    if ( x < -50 || x > 50 ) return ( __builtin_exp( x ) );
    double d;
    *((int*)(&d) + 0) = 0;
    *((int*)(&d) + 1) = (int)(1512775 * x + 1072632447);
    return ( d );
}
#else
inline log_prob_t fast_exp( log_prob_t x )
{
    return ( __builtin_exp( x ) );
}
#endif
