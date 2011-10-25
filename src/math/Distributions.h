/**
    Various distributions, suitable for using in MH sampler.
 */

#pragma once

#include "../BasicTypedefs.h"

#include <vector>

#include <math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>

#include <boost/serialization/version.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/is_bitwise_serializable.hpp>

#include "logmath.h"

inline log_prob_t gaussian_pdf_ln( signal_t x, signal_t sigma )
{
    return ( -0.5 * ( gsl_pow_2( x / sigma ) + M_LNPI + M_LN2 ) - log( sigma ) );
}

inline log_prob_t gamma_pdf_ln( signal_t x, signal_t shape, signal_t scale )
{
    return ( x > 0 ? ( shape - 1 ) * log( x ) - x / scale - shape * log( scale ) - gsl_sf_lngamma( shape ) : gsl_neginf() );
}

/**
    Generates random inverse-Gamma distributed value.
 */
inline signal_t ran_inverse_gamma( const gsl_rng* r, signal_t shape, signal_t scale )
{
    return ( 1.0 / gsl_ran_gamma( r, shape, 1.0 / scale ) );
}

/**
 *  Logarithm of PDF and CDF.
 */
struct LnCdfPdf {
    log_prob_t  lnPdf;      // X == n
    log_prob_t  lnCdfP;     // X <= n
    log_prob_t  lnCdfQ;     // X >= n

    LnCdfPdf()
    : lnPdf( 0 ), lnCdfP( 0 ), lnCdfQ( 0 )
    {}

    LnCdfPdf( log_prob_t lnPdf, log_prob_t lnCdfP, log_prob_t lnCdfQ )
    : lnPdf( 0 ), lnCdfP( 0 ), lnCdfQ( 0 )
    {}

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( lnPdf );
        ar & BOOST_SERIALIZATION_NVP( lnCdfP );
        ar & BOOST_SERIALIZATION_NVP( lnCdfQ );
    }
};

// declare that fast memory-dump type serialization is possible for LnCdfPdf
BOOST_IS_BITWISE_SERIALIZABLE( LnCdfPdf )

template<class Value>
struct DegeneratedDistribution
{
    double operator()( const Value& value ) const
    {
        return ( 0 );
    }
};

struct GaussianGammaJointDistribution
{
    bool        useGaussian;
    signal_t    gaussianMean;
    signal_t    gaussianSigma;
    signal_t    gammaShape;
    signal_t    gammaScale;
    signal_t    curOther;

    GaussianGammaJointDistribution(
        bool    useGaussian,
        signal_t  gaussianMean,
        signal_t  gaussianSigma,
        signal_t  gammaShape,
        signal_t  gammaScale,
        signal_t  curOther
    ) : useGaussian( useGaussian )
      , gaussianMean( gaussianMean ), gaussianSigma( gaussianSigma )
      , gammaShape( gammaShape )
      , gammaScale( gammaScale )
      , curOther( curOther )
    {}

    log_prob_t operator()( signal_t x ) const
    {
        signal_t a = ( useGaussian ? x : curOther ) - gaussianMean;
        signal_t b = ( useGaussian ? x - curOther : curOther - x );
        return ( gaussian_pdf_ln( a, gaussianSigma )
                 + gamma_pdf_ln( b, gammaShape, gammaScale ) );
    }
};

struct GaussianDistribution
{
    signal_t  mean;
    signal_t  sigma;

    GaussianDistribution( signal_t mean, signal_t sigma )
    : mean( mean ), sigma( sigma )
    {}

    log_prob_t operator()( signal_t value ) const
    {
        return ( gaussian_pdf_ln( value - mean, sigma ) );
    }

    log_prob_t generate( const gsl_rng* rng ) const
    {
        return ( mean + gsl_ran_gaussian( rng, sigma ) );
    }

    bool operator==( const GaussianDistribution& that ) const {
        return ( mean == that.mean && sigma == that.sigma );
    }
    bool operator!=( const GaussianDistribution& that ) const {
        return ( !operator==( that ) );
    }

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( mean );
        ar & BOOST_SERIALIZATION_NVP( sigma );
    }

    GaussianDistribution posterior( signal_t samplesMean, size_t samplesCount ) const {
        if ( samplesCount == 0 ) return ( *this );
        return ( GaussianDistribution( 
                ( samplesCount * samplesMean + mean ) / ( 1 + samplesCount ),
                sigma / sqrt( 1 + samplesCount ) ) );
    }
};


class LogNormalDistribution
{
    signal_t  zeta;
    signal_t  sigma;

public:
    LogNormalDistribution( signal_t zeta, signal_t sigma )
    : zeta( zeta ), sigma( sigma )
    {
    }

    log_prob_t operator()( signal_t value ) const
    {
        return ( log( gsl_ran_lognormal_pdf( value, zeta, sigma ) ) );
    }
};

class WeibullDistribution
{
    signal_t  a;
    signal_t  b;

public:
    WeibullDistribution( signal_t a, signal_t b )
    : a( a ), b( b )
    {
    }

    log_prob_t operator()( signal_t value ) const
    {
        return ( log( gsl_ran_weibull_pdf( value, a, b ) ) );
    }
};

struct BernoulliDistribution
{
    prob_t      successRate;
    log_prob_t  lnSuccessRate;
    log_prob_t  lnFailureRate;

    BernoulliDistribution( prob_t successRate,
                           log_prob_t lnSuccessRate,
                           log_prob_t lnFailureRate )
    : successRate( successRate )
    , lnSuccessRate( lnSuccessRate ), lnFailureRate( lnFailureRate )
    {}

    BernoulliDistribution( prob_t successRate = 0.5 )
    : successRate( successRate )
    , lnSuccessRate( ln_1m( 1.0 - successRate ) )
    , lnFailureRate( ln_1m( successRate ) )
    {}

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP( successRate );
        lnSuccessRate = log( successRate );
        lnFailureRate = ln_1m( successRate );
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP( successRate );
    }

    log_prob_t operator()( bool value ) const
    {
        return ( value ? lnSuccessRate : lnFailureRate );
    }

    bool generate( const gsl_rng* rng ) const
    {
        return ( gsl_ran_bernoulli( rng, successRate ) );
    }

    bool operator==( const BernoulliDistribution& that ) const {
        return ( successRate == that.successRate );
    }

    bool operator!=( const BernoulliDistribution& that ) const {
        return ( successRate != that.successRate );
    }

    BernoulliDistribution all_or_nothing( size_t n ) const;
};


class BinomialDistribution
{
private:
    prob_t      successRate;
    log_prob_t  logSuccessRate;
    log_prob_t  logFailureRate;
    size_t      total;

public:
    BinomialDistribution( prob_t successRate, size_t total )
    : successRate( successRate )
    , logSuccessRate( __builtin_log( successRate ) )
    , logFailureRate( fast_log1p( -successRate ) )
    , total( total )
    {}

    log_prob_t operator()( size_t successes ) const;

    size_t generate( const gsl_rng* rng ) const
    {
        return ( gsl_ran_binomial( rng, successRate, total ) );
    }
};

/**
 *  Geometric distribution,
 *  p(k) = p (1-p)^k, k = base,base+1,base+2,...
 */
class GeometricDistribution
{
public:
    size_t base;
    prob_t successRate;

private:
    log_prob_t lnSuccessRate;
    log_prob_t lnFailureRate;

public:
    GeometricDistribution( prob_t successRate = 0.5, prob_t failureRate = 0.5, size_t base = 0 )
    : base( base ), successRate( successRate )
    , lnSuccessRate( log( successRate ) )
    , lnFailureRate( log( failureRate ) )
    {
        BOOST_ASSERT( fabs( failureRate + successRate - 1.0 ) <= std::numeric_limits<prob_t>::epsilon() );
    }

    log_prob_t lnPdf( size_t value ) const
    {
        return ( lnSuccessRate + lnFailureRate * ( value - base ) );
    }
    log_prob_t operator()( size_t value ) const
    {
        return ( lnPdf( value ) );
    }
    log_prob_t lnCdf_P( size_t value ) const
    {
        return ( ln_1mexp( lnCdf_Q( value ) ) );
    }
    log_prob_t lnCdf_Q( size_t value ) const
    {
        return ( lnFailureRate * ( value - base + 1 ) );
    }

    bool generate( const gsl_rng* rng ) const
    {
        return ( gsl_ran_geometric( rng, successRate ) + base - 1 );
    }

    bool operator==( const GeometricDistribution& that ) const {
        return ( successRate == that.successRate && base == that.base );
    }

    bool operator!=( const GeometricDistribution& that ) const {
        return ( successRate != that.successRate || base != that.base );
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP( base );
        ar >> BOOST_SERIALIZATION_NVP( successRate );
        lnSuccessRate = log( successRate );
        lnFailureRate = ln_1m( successRate );
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP( base );
        ar << BOOST_SERIALIZATION_NVP( successRate );
    }

    static GeometricDistribution BySuccessRate( prob_t successRate, size_t base = 0 )
    {
        return ( GeometricDistribution( successRate, 1.0 - successRate, base ) );
    }

    static GeometricDistribution ByFailureRate( prob_t failureRate, size_t base = 0 )
    {
        return ( GeometricDistribution( 1.0 - failureRate, failureRate, base ) );
    }
};

/**
 *  Hypergeometric distribution,
 *  Probability that in a draw of draw_size from group1 + group2 elements
 *  without replacement, there would be group1_draw_size elements from group1.
 * 
 *  @see FischerHypergeometricDistribution
 */
class HypergeometricDistribution
{
public:
    size_t group1_size;
    size_t group2_size;
    size_t draw_size;

public:
    HypergeometricDistribution( size_t group1_size = 0, size_t group2_size = 0, size_t draw_size = 0 )
    : group1_size( group1_size ), group2_size( group2_size ), draw_size( draw_size )
    {
        precalc();
    }

    log_prob_t lnPdf( size_t group1_draw_size ) const;

    log_prob_t operator()( size_t value ) const
    {
        return ( lnPdf( value ) );
    }
    log_prob_t lnCdf_P( size_t group1_draw_size, bool include_bound = false ) const;
    log_prob_t lnCdf_Q( size_t group1_draw_size, bool include_bound = false ) const;

    size_t mode() const {
        return ( ( draw_size + 1 ) * ( group1_size + 1 ) / ( group1_size + group2_size + 2 ) );
    }

    bool generate( const gsl_rng* rng ) const
    {
        return ( gsl_ran_hypergeometric( rng, group1_size, group2_size, draw_size ) );
    }

    bool operator==( const HypergeometricDistribution& that ) const {
        return ( group1_size == that.group1_size 
                 && group2_size == that.group2_size
                 && draw_size == that.draw_size
               );
    }

    bool operator!=( const HypergeometricDistribution& that ) const {
        return ( !operator==( that ) );
    }

     BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP( draw_size );
        ar >> BOOST_SERIALIZATION_NVP( group1_size );
        ar >> BOOST_SERIALIZATION_NVP( group2_size );
        precalc();
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP( draw_size );
        ar << BOOST_SERIALIZATION_NVP( group1_size );
        ar << BOOST_SERIALIZATION_NVP( group2_size );
    }

private:
    log_prob_t  _draw_combinations;

    void precalc();
};

/**
 *  Fischer Noncentral Hypergeometric distribution,
 *  Probability that in a draw of draw_size from group1 + group2 elements
 *  without replacement, there would be group1_draw_size elements from group1,
 *  when there's a bias in taking a ball of specific color.
 *
 *  @see HypergeometricDistribution
 */
struct FisherHypergeometricDistribution {
public:
    size_t      group1_size;
    size_t      group2_size;
    size_t      draw_size;
    log_prob_t  ln_group1_odds;  /** log of odds ratio of taking group1 element */

public:
    FisherHypergeometricDistribution( size_t group1_size = 0, size_t group2_size = 0, size_t draw_size = 0, 
                                       log_prob_t ln_group1_odds = 0.0 )
    : group1_size( group1_size ), group2_size( group2_size ), draw_size( draw_size )
    , ln_group1_odds( ln_group1_odds )
    {
        if ( draw_size > group1_size + group2_size ) {
            THROW_EXCEPTION( std::invalid_argument, "Drawing more (" << draw_size 
                << ") than available (" << group1_size << " and " << group2_size << ")" );
        }
        precalc();
    }

    log_prob_t lnPdf( size_t group1_draw_size ) const
    {
        check_group1_draw_size( group1_draw_size );
        return ( lnTerm( group1_draw_size ) - _normalizer );
    }

    log_prob_t operator()( size_t value ) const
    {
        return ( lnPdf( value ) );
    }
    log_prob_t lnCdf_P( size_t group1_draw_size, bool include_bound = false ) const;
    log_prob_t lnCdf_Q( size_t group1_draw_size, bool include_bound = false ) const;

    size_t mode() const;

    bool operator==( const FisherHypergeometricDistribution& that ) const {
        return ( group1_size == that.group1_size 
                 && group2_size == that.group2_size
                 && draw_size == that.draw_size
                 && ln_group1_odds == that.ln_group1_odds
               );
    }

    bool operator!=( const FisherHypergeometricDistribution& that ) const {
        return ( !operator==( that ) );
    }

     BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP( draw_size );
        ar >> BOOST_SERIALIZATION_NVP( group1_size );
        ar >> BOOST_SERIALIZATION_NVP( group2_size );
        ar >> BOOST_SERIALIZATION_NVP( ln_group1_odds );
        precalc();
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP( draw_size );
        ar << BOOST_SERIALIZATION_NVP( group1_size );
        ar << BOOST_SERIALIZATION_NVP( group2_size );
        ar << BOOST_SERIALIZATION_NVP( ln_group1_odds );
    }

private:
    log_prob_t  _normalizer;

    void check_group1_draw_size( size_t group1_draw_size ) const;

    log_prob_t lnTerm( size_t successes ) const;
    void precalc();
};

/**
 *  Multinomial Hypergeometric distribution,
 * 
 *  @see HypergeometricDistribution
 */
class MultinomialHypergeometricDistribution
{
public:
    std::vector<size_t> group_sizes;
    size_t draw_size;

public:
    MultinomialHypergeometricDistribution( const std::vector<size_t>& group_sizes, size_t draw_size = 0 )
    : group_sizes( group_sizes ), draw_size( draw_size )
    {
        precalc();
    }

    log_prob_t lnPdf( const std::vector<size_t>& group_draw_size ) const;

    bool operator==( const MultinomialHypergeometricDistribution& that ) const {
        return ( group_sizes == that.group_sizes
                 && draw_size == that.draw_size );
    }

    bool operator!=( const MultinomialHypergeometricDistribution& that ) const {
        return ( !operator==( that ) );
    }

     BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP( draw_size );
        ar >> BOOST_SERIALIZATION_NVP( group_sizes );
        precalc();
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP( draw_size );
        ar << BOOST_SERIALIZATION_NVP( group_sizes );
    }

private:
    log_prob_t  _draw_combinations;

    void precalc();
};

struct GammaDistribution {
    signal_t    shape;
    signal_t    scale;
    log_prob_t  lnScale;
    log_prob_t  lnGammaShape;

    GammaDistribution( signal_t shape, signal_t scale )
    : shape( shape ), scale( scale )
    {
        BOOST_ASSERT( shape > 0 );
        BOOST_ASSERT( scale > 0 );
        lnScale = log( scale );
        lnGammaShape = gsl_sf_lngamma( shape );
    }

    signal_t generate( const gsl_rng* rng ) const
    {
        return ( gsl_ran_gamma( rng, shape, scale ) );
    }

    signal_t mean() const {
        return ( shape * scale );
    }

    /**
     *  log( x^{shape-1}/exp(x/shape)/scale^shape/Gamma(shape) )
     */
    log_prob_t operator()( signal_t x ) const {
        return ( x > 0 ? ( shape - 1 ) * log( x ) - x / scale - shape * lnScale - lnGammaShape 
                       : gsl_neginf() );
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP( shape );
        ar >> BOOST_SERIALIZATION_NVP( scale );
        lnScale = log( scale );
        lnGammaShape = gsl_sf_lngamma( shape );
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP( shape );
        ar << BOOST_SERIALIZATION_NVP( scale );
    }

    GammaDistribution poissonPosterior( size_t totalCounts, size_t nObservations ) const;
};

/**
 *  Beta Distribution.
 *  It's specified in terms of observed successes and failures.
 *  The relation with classic beta definition:
 *     Beta(\alpha,\beta) = Beta(successes+1, failures+1)
 */
struct BetaDistribution {
    prob_t      successes;
    prob_t      failures;
    log_prob_t  lnBeta;

    /**
     *  Logarithm of Beta(\alpha,\beta) function.
     */
    static double LnBeta( prob_t alpha, prob_t beta )
    {
        return ( gsl_sf_lngamma( beta ) - gsl_sf_lnpoch( alpha, beta ) );
    }

    BetaDistribution( prob_t successes, prob_t failures )
    : successes( successes ), failures( failures )
    {
        BOOST_ASSERT( successes > 0 );
        BOOST_ASSERT( failures > 0 );
        lnBeta = LnBeta( successes + 1, failures + 1 );
    }

    prob_t generate( const gsl_rng* rng ) const
    {
        return ( gsl_ran_beta( rng, successes + 1, failures + 1 ) );
    }

    log_prob_t operator()( prob_t x ) const {
        return ( x > 0 && x < 1 ? successes * log( x ) + failures * gsl_log1p( -x ) - lnBeta
                                : gsl_neginf() );
    }

    typedef std::vector<signal_t> data_vector_t;

    BetaDistribution posterior( const data_vector_t& data, signal_t base = 0.0 ) const;

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar >> BOOST_SERIALIZATION_NVP( successes );
        ar >> BOOST_SERIALIZATION_NVP( failures );
        lnBeta = LnBeta( successes + 1, failures + 1 );
    }

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << BOOST_SERIALIZATION_NVP( successes );
        ar << BOOST_SERIALIZATION_NVP( failures );
    }
};

struct TransformedGammaDistribution: GammaDistribution {
    signal_t  mult;
    signal_t  add;

    TransformedGammaDistribution( const GammaDistribution& gamma, signal_t mult = 1.0, signal_t add = 0.0 )
    : GammaDistribution( gamma ), mult( mult ), add( add )
    {
        BOOST_ASSERT( mult != 0 );
    }

    signal_t generate( const gsl_rng* rng ) const
    {
        return ( mult * GammaDistribution::generate( rng ) + add );
    }

    log_prob_t operator()( signal_t value ) const {
        return ( GammaDistribution::operator()( ( value - add ) / mult ) );
    }
};

/**
 *  Inverse Gamma distribution.
 *  If X is Gamma-distributed, then 1/X is inverse-Gamma distributed.
 * 
 *  @see http://en.wikipedia.org/wiki/Inverse-gamma_distribution
 */
struct InverseGammaDistribution {
    signal_t    shape;
    signal_t    scale;

    InverseGammaDistribution( signal_t shape, signal_t scale )
    : shape( shape ), scale( scale )
    {
        BOOST_ASSERT( shape > 0 );
        BOOST_ASSERT( scale > 0 );
    }

    typedef std::vector<signal_t> data_vector_t;

    InverseGammaDistribution posterior( signal_t priorDataMean, const data_vector_t& data ) const;

    signal_t generate( const gsl_rng* rng ) const
    {
        signal_t res = ran_inverse_gamma( rng, shape, scale );
        LOG_DEBUG2( "IGamma(" << shape << "," << scale << ")~" << res );
        return ( res );
    }


    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( shape );
        ar & BOOST_SERIALIZATION_NVP( scale );
    }
};

/**
 *  Normal-scaled inverse Gamma distribution.
 *  Distribution over normal distributions.
 * 
 *  @see http://en.wikipedia.org/wiki/Normal-scaled_inverse_gamma_distribution
 */
struct NormalScaledInverseGammaDistribution {
    signal_t    meanMean;           /** lambda */
    signal_t    meanVarScale;       /** scale (denominator) of mean variance */
    InverseGammaDistribution  varDistrib; /** distribution of unscaled variance */

    NormalScaledInverseGammaDistribution( 
        signal_t meanMean, signal_t meanVarScale,
        signal_t varShape, signal_t varScale
    ) : meanMean( meanMean ), meanVarScale( meanVarScale )
    , varDistrib( varShape, varScale )
    {
        BOOST_ASSERT( meanVarScale > 0 );
    }

    typedef std::vector<signal_t> data_vector_t;

    NormalScaledInverseGammaDistribution posterior( const data_vector_t& data ) const;

    GaussianDistribution generate( const gsl_rng* rng ) const
    {
        signal_t var = varDistrib.generate( rng );
        signal_t mean = GaussianDistribution( meanMean, sqrt( var / meanVarScale ) ).generate( rng );
        return ( GaussianDistribution( mean, sqrt( var ) ) );
    }


    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( meanMean );
        ar & BOOST_SERIALIZATION_NVP( meanVarScale );
        ar & BOOST_SERIALIZATION_NVP( varDistrib );
    }
};
