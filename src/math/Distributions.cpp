#include <numeric>

#include "Distributions.h"

#include <math.h>

#include "logmath.h"

#if 0
/**
    Generate discrete random value according to the vector of its probabilities.
 */
int gsl_ran_probability_vector(
    const gsl_rng*              r,          /**< @param[in] random number generator */
    const probability_vector_t& prob_vec    /**< @param[in] vector of discrete probabilities (may be unnormalized) */
){
    if ( prob_vec.size() == 0 )     return ( -1 );

    double  sum = 0.0;
    for ( int i = 0; i < prob_vec.size(); i++ ) {
        sum += prob_vec[ i ];
    }

    double  rnd_val = gsl_ran_flat( r, 0.0, sum );
    for ( int i = 0; i < prob_vec.size(); i++ ) {
        const double pdf = prob_vec[ i ];
        if ( rnd_val <= pdf )   return ( i );
        rnd_val -= pdf;
    }
    
    /// @warning should never reach here
    return ( prob_vec.size() - 1 );
}
#endif

log_prob_t BinomialDistribution::operator()( size_t successes ) const
{
    return ( ln_binomial_coef( total, successes )
             + logSuccessRate * successes
             + logFailureRate * ( total - successes ) );
}

/**
 *  Posterior distribution for lambda parameter of classic Poisson.
 *  shape' = shape + \sum ob
 */
GammaDistribution GammaDistribution::poissonPosterior(
    size_t totalCounts,     /** sum of counts in all observations */
    size_t nObservations    /** number of observations */
) const {
    return ( GammaDistribution( shape + totalCounts,
                                scale / ( 1.0 + scale * nObservations ) ) );
}

log_prob_t HypergeometricDistribution::lnPdf( size_t group1_draw_size ) const
{
    return ( ln_binomial_coef( group1_size, group1_draw_size )
            + ln_binomial_coef( group2_size, draw_size - group1_draw_size )
            - _draw_combinations );
}

log_prob_t HypergeometricDistribution::lnCdf_P(
    size_t group1_draw_size,
    bool include_bound
) const {
    const size_t min_succ_draw = draw_size < group1_size ? 0 : draw_size - group1_size;
    const size_t max_succ_draw = group1_draw_size + ( include_bound ? 1 : 0 );
    if ( min_succ_draw == max_succ_draw ) return ( -std::numeric_limits<log_prob_t>::infinity() );
    log_prob_t maxPdf = lnPdf( std::max( mode(), group1_draw_size ) );
    double res = 0;
    for ( size_t succ_group1_draw_size = min_succ_draw;
          succ_group1_draw_size < max_succ_draw;
          succ_group1_draw_size++
    ){
        res += exp( lnPdf( succ_group1_draw_size ) - maxPdf );
    }
    return ( log( res ) + maxPdf );
}

log_prob_t HypergeometricDistribution::lnCdf_Q(
    size_t group1_draw_size,
    bool include_bound
) const {
    const size_t min_succ_draw = group1_draw_size + ( include_bound ? 0 : 1 );
    const size_t max_succ_draw = std::min( draw_size, group1_size ) + 1;
    if ( min_succ_draw == max_succ_draw ) return ( -std::numeric_limits<log_prob_t>::infinity() );
    log_prob_t maxPdf = lnPdf( std::max( mode(), group1_draw_size ) );
    double res = 0;
    for ( size_t succ_group1_draw_size = min_succ_draw;
          succ_group1_draw_size < max_succ_draw;
          succ_group1_draw_size++
    ){
        res += exp( lnPdf( succ_group1_draw_size ) - maxPdf );
    }
    return ( log( res ) + maxPdf );
}

void HypergeometricDistribution::precalc()
{
    _draw_combinations = ln_binomial_coef( group1_size + group2_size, draw_size );
}

void FisherHypergeometricDistribution::check_group1_draw_size(
    size_t group1_draw_size
) const {
    if ( group1_draw_size > draw_size ) {
        THROW_EXCEPTION( std::invalid_argument, "Drawing more (" << group1_draw_size
            << ") from group1 than drawing in total (" << draw_size << ")" );
    }
    if ( group1_draw_size > group1_size ) {
        THROW_EXCEPTION( std::invalid_argument, "Drawing more (" << group1_draw_size 
            << ") from group1 than available (" << group1_size << ")" );
    }
    if ( draw_size - group1_draw_size > group2_size ) {
        THROW_EXCEPTION( std::invalid_argument, "Drawing more (" << draw_size - group1_draw_size
            << ") from group2 than available (" << group2_size << ")" );
    }
}

void FisherHypergeometricDistribution::precalc()
{
    size_t minSucc = draw_size > group2_size ? draw_size - group2_size : 0;
    size_t maxSucc = std::min( draw_size, group1_size );
    double maxLnTerm = lnTerm( std::min( std::max( minSucc, mode() ), maxSucc ) );
    double res = 0.0;
    for ( size_t successes = minSucc; successes <= maxSucc; successes++ ){
        res += exp( lnTerm( successes ) - maxLnTerm );
    }
    _normalizer = maxLnTerm + log( res );
}

/**
 *  Mode of Fisher Noncentral Hypergeometric Distribution.
 *  See Wikipedia.
 */
size_t FisherHypergeometricDistribution::mode() const
{
    double ratio = exp( ln_group1_odds );
    double A = gsl_expm1( ln_group1_odds );
    double B = draw_size - group2_size - ( draw_size + group1_size + 2 ) * ratio;
    double C = ( draw_size + 1 ) * ( group1_size + 1 ) * ratio;
    size_t res = (size_t)floor( 2 * C / ( sqrt( B * B - 4 * A * C ) - B ) );
    return ( std::min( std::max( draw_size > group2_size ? draw_size - group2_size : (size_t)0, res ),
                       std::min( draw_size, group1_size ) ) );
}

log_prob_t FisherHypergeometricDistribution::lnTerm(
    size_t successes
) const {
    return (   ln_binomial_coef( group1_size, successes )
             + ln_binomial_coef( group2_size, draw_size - successes )
             + ln_group1_odds * successes );
}

log_prob_t FisherHypergeometricDistribution::lnCdf_P(
    size_t group1_draw_size,
    bool include_bound
) const {
    check_group1_draw_size( group1_draw_size );
    const size_t min_succ_draw = draw_size < group1_size ? 0 : draw_size - group1_size;
    const size_t max_succ_draw = group1_draw_size + ( include_bound ? 1 : 0 );
    if ( min_succ_draw == max_succ_draw ) return ( -std::numeric_limits<log_prob_t>::infinity() );
    log_prob_t maxLnTerm = lnTerm( std::max( mode(), group1_draw_size ) );
    double res = 0;
    for ( size_t succ_group1_draw_size = min_succ_draw;
          succ_group1_draw_size < max_succ_draw;
          succ_group1_draw_size++
    ){
        res += exp( lnTerm( succ_group1_draw_size ) - maxLnTerm );
    }
    return ( log( res ) + maxLnTerm - _normalizer );
}

log_prob_t FisherHypergeometricDistribution::lnCdf_Q(
    size_t  group1_draw_size,
    bool    include_bound
) const {
    check_group1_draw_size( group1_draw_size );
    const size_t min_succ_draw = group1_draw_size + ( include_bound ? 0 : 1 );
    const size_t max_succ_draw = std::min( draw_size, group1_size ) + 1;
    if ( min_succ_draw == max_succ_draw ) return ( -std::numeric_limits<log_prob_t>::infinity() );
    log_prob_t maxLnTerm = lnTerm( std::max( mode(), group1_draw_size ) );
    double res = 0;
    for ( size_t succ_group1_draw_size = min_succ_draw;
          succ_group1_draw_size < max_succ_draw;
          succ_group1_draw_size++
    ){
        res += exp( lnTerm( succ_group1_draw_size ) - maxLnTerm );
    }
    return ( log( res ) + maxLnTerm - _normalizer );
}

log_prob_t MultinomialHypergeometricDistribution::lnPdf(
    const std::vector<size_t>&  group_draw_size
) const {
    log_prob_t  res = 0;
    for ( size_t i = 0; i < group_draw_size.size(); i++ ) {
        res += ln_binomial_coef( group_sizes[ i ], group_draw_size[ i ] );
    }
    return ( res - _draw_combinations );
}

void MultinomialHypergeometricDistribution::precalc()
{
    size_t  total = 0;
    for ( size_t i = 0; i < group_sizes.size(); i++ ) {
        total += group_sizes[ i ];
    }
    _draw_combinations = ln_binomial_coef( total, draw_size );
}

/**
 *  Build new Bernoulli trial, where all n of given
 *  Bernoulli trials either succeed or fail.
 */
BernoulliDistribution BernoulliDistribution::all_or_nothing(
    size_t n
) const {
    if ( lnSuccessRate >= lnFailureRate ) {
        log_prob_t lnSucc = -ln_1pexp( n * ( lnFailureRate - lnSuccessRate ) );
        return ( BernoulliDistribution( exp( lnSucc ), lnSucc, ln_1mexp( lnSucc ) ) );
    } else {
        log_prob_t lnFail = -ln_1pexp( n * ( lnSuccessRate - lnFailureRate ) );
        log_prob_t lnSucc = ln_1mexp( lnFail );
        return ( BernoulliDistribution( exp( lnSucc ), lnSucc, lnFail ) );
    }
}

/**
 *  Posterior Beta distribution for geometrically-distributed variable.
 */
BetaDistribution BetaDistribution::posterior(
    const data_vector_t& data,
    signal_t base
) const {
    if ( data.size() == 0 ) return ( *this );

    double dataSum = 0.0;
    for ( std::size_t i = 0; i < data.size(); i++ ) {
        dataSum += data[i] - base;
    }
    return ( BetaDistribution( successes + data.size(), failures + dataSum ) );
}

/**
    Gets posterior parameters for the inverse-gamma distributed 
    variance of Normally-distributed value.
 */
InverseGammaDistribution InverseGammaDistribution::posterior(
    signal_t                priorDataMean,  /**< @param[in] prior mean value of data */
    const data_vector_t&    data            /**< @param[in] observed normally-distributed data vector with known mean */
) const {
    double var = 0.0;
    for ( std::size_t i = 0; i < data.size(); i++ ) {
        const double delta = data[ i ] - priorDataMean;
        var += delta * delta;
    }
    InverseGammaDistribution res = InverseGammaDistribution( shape + 0.5 * data.size(), scale + 0.5 * var );
    LOG_DEBUG2( "IGamma(" << shape << "," << scale << ") posterior for " << data.size() << "," << var 
                << " is IGamma(" << res.shape << "," << res.scale << ")" );
    return ( res );
}

NormalScaledInverseGammaDistribution NormalScaledInverseGammaDistribution::posterior(
    const data_vector_t& data
) const {
    if ( data.size() == 0 ) return ( *this );

    double dataSqrSum = 0.0;
    double dataSum = 0.0;
    for ( std::size_t i = 0; i < data.size(); i++ ) {
        dataSqrSum += gsl_pow_2( data[i] );
        dataSum += data[i];
    }
    return ( NormalScaledInverseGammaDistribution( 
        ( meanMean * meanVarScale + dataSum ) / ( meanVarScale + data.size() ),
        meanVarScale + data.size(),
        varDistrib.shape + 0.5 * data.size(),
        varDistrib.scale + 0.5 * ( dataSqrSum - gsl_pow_2( dataSum ) / data.size() ) 
                         + 0.5 * data.size() * meanVarScale / ( data.size() + meanVarScale ) * gsl_pow_2( dataSum / data.size() - meanMean )
    ) );
}
