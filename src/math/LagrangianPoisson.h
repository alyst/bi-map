#pragma once

#include "../BasicTypedefs.h"

#include <gsl/gsl_poly.h>

#include <boost/serialization/version.hpp>
#include <boost/serialization/nvp.hpp>

#include "logmath.h"

/**
    Lagrangian Poisson Distribution aka Poisson-Consul Distribution.

    pdf = mu e^(-x lambda-mu) (x lambda+mu)^(-1+x)/(x!),
    ln pdf = -x lambda -mu + (-1+x)log( x lambda + mu ) - ln_factorial(x)
    where mu = rate, lambda = shape.

    When shape = 0 becomes standard Poisson distribution.
    When shape < 0 is underdispersed (not implemented).
    When shape > 0 is overdispersed.

    See [1] http://reference.wolfram.com/mathematica/ref/PoissonConsulDistribution.html
    See [2] "Univariate discrete distributions"
        By Norman Lloyd Johnson, Adrienne W. Kemp, Samuel Kotz,
        3rd  Ed., Wiley, p.336
    See [3] "Lagrangian probability distributions"
        By P.C.Consul, Felix Famoye, Section 9
*/
struct LagrangianPoissonDistribution {
    log_prob_t  lnRate;    	/** log(_rate) */
    log_prob_t  shape;      /** pdf's shape, the more - the steepier */
    mutable log_prob_t _rate; /** exp(lnRate) */
    mutable int _mode;

    LagrangianPoissonDistribution( log_prob_t lnRate = 0, log_prob_t shape = 0 )
    : lnRate( lnRate ), shape( shape )
    , _rate( __builtin_exp( lnRate ) ), _mode( -1 )
    {
        if ( shape >= 1.0 ) throw std::invalid_argument( "Shape >= 1" );
        else if ( shape < 0 ) throw std::invalid_argument( "Shape < 0 not supported" );
        // shape <= 0 (i.e. truncated LPD) exists, but not supported here
    }

    void swap( LagrangianPoissonDistribution& that )
    {
        std::swap( lnRate, that.lnRate );
        std::swap( shape, that.shape );
        std::swap( _rate, that._rate );
    }

    bool operator==( const LagrangianPoissonDistribution& that ) const {
        return ( lnRate == that.lnRate && shape == that.shape );
    }

    inline double lnPdf( signal_t x ) const;
    double pdf( signal_t x ) const {
        return ( exp( lnPdf( x ) ) );
    }

    log_prob_t lnCdf_P( signal_t x ) const {
        BOOST_ASSERT( x >= -1 );
        if ( x == -1 ) {
            return ( gsl_neginf() );
        }
        else if ( x == 0 ) {
            return ( lnPdf( x ) );
        }
        else if ( shape == 0 ) {
            // if shape > 0 this is incorrect
            return ( ln_gammareg_inc_Q( x+1, _rate, lnRate ) );
        }
        else {
            throw std::runtime_error( "CDF of LPD with shape!=0 not implemented" );
        }
    }
    log_prob_t lnCdf_Q( signal_t x ) const {
        BOOST_ASSERT( x >= -1 );
        if ( x == -1 ) {
            return ( 0 );
        }
        else if ( shape == 0 ) {
            // if shape > 0 this is incorrect
            return ( ln_gammareg_inc_P( (int)(x+1), _rate, lnRate ) );
        }
        else {
            throw std::runtime_error( "CDF of LPD with shape!=0 not implemented" );
        }
    }

    signal_t lnMean() const {
        return ( lnRate - log( 1-shape ) );
    }
    signal_t mean() const {
        return ( _rate / ( 1.0-shape ) );
    }

    int mode() const {
        if ( shape == 0 ) {
            return ( _rate <= std::numeric_limits<int>::max()
                     ? _rate : std::numeric_limits<int>::max() );
        } else if ( _mode >= 0 ) {
            return ( _mode );
        } else {
            // see [3], p.170
            double expShape = exp( shape );
            double mode1 = ( _rate - 1.0 / expShape ) / ( expShape - 2 * shape );
            double mode2 = mode1;
            double tmp;
            gsl_poly_solve_quadratic( shape*shape,
                                      2 * shape * _rate - ( _rate + 2 * shape ) * expShape,
                                      _rate*_rate, &mode2, &tmp );
            log_prob_t maxLnPdf = -std::numeric_limits<log_prob_t>::infinity();
            for ( int i = (int)std::min(mode1, mode2); i <= (int)std::max(mode1, mode2); i++ ) {
                double lnpdf = lnPdf( i );
                if ( lnpdf > maxLnPdf ) {
                    _mode = i;
                    maxLnPdf = lnpdf;
                }
            }
            return ( _mode );
        }
    };

    double operator()( signal_t x ) const {
        return ( lnPdf( x ) );
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load( Archive & ar, const unsigned int version )
    {
        ar >> BOOST_SERIALIZATION_NVP( lnRate );
        ar >> BOOST_SERIALIZATION_NVP( shape );
        _rate = __builtin_exp( lnRate );
        _mode = -1;
    }

    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
    {
        ar << BOOST_SERIALIZATION_NVP( lnRate );
        ar << BOOST_SERIALIZATION_NVP( shape );
    }

    template<typename _CharT, typename _Traits>
    friend std::basic_ostream<_CharT,_Traits>& operator<<(
        std::basic_ostream<_CharT,_Traits>&     out,
        const LagrangianPoissonDistribution&    a
    ){
        return ( out << '(' << a.lnRate << ", " << a.shape << ')' );
    }

    int supportMin() const {
        return ( 0 );
    }
    int supportMax() const {
        return ( std::numeric_limits<int>::max() );
    }
};

inline log_prob_t LagrangianPoissonDistribution::lnPdf(
    signal_t x
) const {
    BOOST_ASSERT( x >= 0 );
    if ( x <= 1 ) {
        // ln[pdf](0) = -rate
        // ln[pdf](1) = ln(rate) - rate - shape
        return ( ( lnRate - shape ) * x - _rate );
    }
    else if ( shape == 0 ) {
        // Poisson distribution
        return ( lnRate <= 600
                    ? lnRate * x - _rate - ln_factorial( x )
                    // extremly large rate
                    : ( lnRate - __builtin_exp( lnRate - __builtin_log( x ) ) ) * x - ln_factorial( x ) );
    }
    else {
        double f = shape * x + _rate;
        return ( lnRate - f + ( x - 1 ) * __builtin_log( f )  - ln_factorial( x ) );
    }
}

struct ZeroInflatedLagrangianPoissonDistribution: public LagrangianPoissonDistribution {
    typedef LagrangianPoissonDistribution base_type;

    log_prob_t  _inflationFactor;
    log_prob_t  _lnPdfZero;
    log_prob_t  _lnPdfNorm;

    ZeroInflatedLagrangianPoissonDistribution( log_prob_t lnRate = 0, log_prob_t shape = 0,
                                               log_prob_t inflationFactor = 0 )
    : base_type( lnRate, shape )
    , _inflationFactor( inflationFactor )
    {
        updateConstants();
    }

    void updateConstants()
    {
        if ( _inflationFactor > 0 ) {
            double lnPdfZeroAdd = base_type::lnPdf( _inflationFactor * base_type::mode() );
            _lnPdfNorm = logsumexp( lnPdfZeroAdd, 0.0 );
            _lnPdfZero = logsumexp( lnPdfZeroAdd, base_type::lnPdf( 0 ) );
        } else {
            _lnPdfZero = base_type::lnPdf( 0 );
            _lnPdfNorm = 0;
        }
    }

    bool operator==( const ZeroInflatedLagrangianPoissonDistribution& that ) const {
        return ( base_type::operator==( that ) && _inflationFactor == that._inflationFactor );
    }

    double lnPdf( signal_t x ) const {
        return ( ( x > 0 ? base_type::lnPdf( x ) : _lnPdfZero ) - _lnPdfNorm );
    }
    double pdf( signal_t x ) const {
        return ( exp( lnPdf( x ) ) );
    }

    double operator()( signal_t x ) const {
        return ( lnPdf( x ) );
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER();

    template<class Archive>
    void load( Archive & ar, const unsigned int version )
    {
        ar >> boost::serialization::make_nvp( "lp", boost::serialization::base_object<base_type>( *this ) );
        ar >> boost::serialization::make_nvp( "inflationFactor", _inflationFactor );
        updateConstants();
    }

    template<class Archive>
    void save( Archive & ar, const unsigned int version ) const
    {
        ar << boost::serialization::make_nvp( "lp", boost::serialization::base_object<base_type>( *this ) );
        ar << boost::serialization::make_nvp( "inflationFactor", _inflationFactor );
    }

    template<typename _CharT, typename _Traits>
    friend std::basic_ostream<_CharT,_Traits>& operator<<(
        std::basic_ostream<_CharT,_Traits>&     out,
        const ZeroInflatedLagrangianPoissonDistribution&    a
    ){
        return ( out << '(' << a.lnRate << ", " << a.shape
                     << ", " << a._inflationFactor << ')' );
    }

    int supportMin() const {
        return ( 0 );
    }
    int supportMax() const {
        return ( std::numeric_limits<int>::max() );
    }
};
