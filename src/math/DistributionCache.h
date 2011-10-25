#pragma once

#include <vector>
#include <stdexcept>

#include "../BasicTypedefs.h"
#include "logmath.h"

/**
 *  Cache of one-parameter family of distributions.
 *  Supports storing ln(PDF) and ln(Q-value) of distribution.
 *
 *  @remark It's assumed (by qvalue calculation) that Distribution is unimodal.
 */
template<typename DistributionFactory, typename Param=double, typename Value=int>
class DistributionCache {
public:
    typedef DistributionFactory distrib_factory_type;
    typedef typename DistributionFactory::distribution_type distrib_type;
    typedef Param param_type;
    typedef Value value_type;
    typedef log_prob_t log_prob_type;
    typedef DistributionCache<DistributionFactory, Param, Value> distrib_cache_type;

    DistributionCache( const distrib_factory_type& distrFactory,
                       const value_type& minValue, const value_type& maxValue,
                       const param_type& paramStep )
    : _distrFactory( distrFactory )
    , _minValue( minValue ), _maxValue( maxValue )
    , _paramStep( paramStep )
    {}

    ~DistributionCache()
    {
        // free cached tables
        for ( size_t i = 0; i < _data.size(); i++ ) {
            if ( _data[i] ) delete [] _data[i];
        }
    }

    log_prob_type lnPdf( const param_type& param, const value_type& val );
    log_prob_type lnNormPdf( const param_type& param, const value_type& val );

private:
    typedef int param_ix_type;
    struct CacheItem {
        log_prob_type   lnPdf;
        log_prob_type   lnNormPdf;

        CacheItem( log_prob_type lnPdf = unset(), log_prob_type lnMaxPdf = unset() )
        : lnPdf( lnPdf ), lnNormPdf( is_finite( lnPdf ) && is_finite( lnMaxPdf ) ? lnPdf - lnMaxPdf : unset() )
        {};
    };
    typedef CacheItem cache_item_type;
    typedef std::vector<cache_item_type*> distrib_container;

    void evalTable( param_ix_type parmIx );
    void checkValue( const value_type& value ) const;
    const cache_item_type* ensureTable( const param_type& param );

    distrib_factory_type   _distrFactory;
    value_type      _minValue;
    value_type      _maxValue;
    param_type      _paramStep;
    param_type      _minParam;

    distrib_container    _data;

public:
    class CachedDistributionRef {
    private:
        friend class DistributionCache<DistributionFactory,Param,Value>;

        const param_type    _param;
        const distrib_cache_type&   _cache;
        const cache_item_type*      _distribTable;

    public:
        log_prob_t lnPdf( const value_type& value ) const {
            _cache.checkValue( value );
            return ( _distribTable[ value - _cache._minValue ].lnPdf );
        }
        log_prob_t lnNormPdf( const value_type& value ) const {
            _cache.checkValue( value );
            return ( _distribTable[ value - _cache._minValue ].lnNormPdf );
        }
    private:
        CachedDistributionRef( const distrib_cache_type& cache,
                               const param_type& param,
                               const cache_item_type* distribTable )
        : _param( param ), _cache( cache ), _distribTable( distribTable )
        {}
    };

    size_t valuesCount() const  {
        return ( _maxValue - _minValue + 1 );
    }
    typedef CachedDistributionRef const_distrib_ref_type;

    const_distrib_ref_type distrib( const param_type& param )
    {
        return ( const_distrib_ref_type( *this, param, ensureTable( param ) ) );
    }
};

template<typename DistributionFactory, typename Param, typename Value>
void DistributionCache<DistributionFactory,Param,Value>::evalTable(
    param_ix_type parmIx
){
    param_type param = _minParam + parmIx * _paramStep;
    distrib_type distr = _distrFactory( param );
    // fill lnPdf and inverse lnPdf multimap
    cache_item_type* distribTable = _data[ parmIx ];
    double mode = distr.mode();
    double maxLnPdf = distr.lnPdf( mode );
    BOOST_ASSERT( distribTable != NULL );
    for ( value_type val = _minValue; val <= _maxValue; ++val ) {
        log_prob_type lnpdf = distr.lnPdf( val );
        distribTable[ val - _minValue ].lnPdf = lnpdf;
        if ( lnpdf > maxLnPdf ) maxLnPdf = lnpdf;
    }
    for ( int i = 0; i <= _maxValue - _minValue; i++ ) {
        distribTable[i].lnNormPdf = distribTable[i].lnPdf - maxLnPdf;
    }
#if DEBUG_LEVEL >= 3
    LOG_DEBUG3( "Table for param p=" << param );
    for ( size_t i = 0; i <= _maxValue - _minValue; i++ ) {
        LOG_DEBUG3( "lnPdf[" << _minValue + i << "]=" << distribTable[i].lnPdf );
        LOG_DEBUG3( "lnNormPdf[" << _minValue + i << "]=" << distribTable[i].lnNormPdf );
    }
#endif
}

template<typename DistributionFactory, typename Param, typename Value>
const typename DistributionCache<DistributionFactory,Param,Value>::cache_item_type*
DistributionCache<DistributionFactory,Param,Value>::ensureTable(
    const param_type& param
){
    if ( is_unset( param ) ) {
        THROW_EXCEPTION( std::invalid_argument, "Parameter is NaN" );
    }
    // no data yet, initialize cache
    if ( _data.empty() ) {
        // initialize
        param_ix_type paramIx = ( param / _paramStep );
        _minParam = paramIx * _paramStep;
        _data.resize( 1, NULL );
    }
    param_ix_type paramIx = ( param - _minParam ) / _paramStep;
    // extend the array
    if ( paramIx < 0 ) {
        _data.insert( _data.begin(), -paramIx, NULL );
        _minParam += paramIx * _paramStep;
        paramIx = 0;
    }
    else if ( paramIx >= (int)_data.size() ) {
        _data.insert( _data.end(), paramIx + 1 - _data.size(), NULL );
    }
    if ( _data[ paramIx ] == NULL ) {
        // calculate if required
        LOG_DEBUG3( "Allocating memory for distribution"
                    << " param=" << param
                    << " minVal=" << _minValue
                    << " maxVal=" << _maxValue
                    << " datasize=" << _data.size()
                    << " allocsize=" << ( _maxValue - _minValue + 1 ) );
#if DEBUG_LEVEL >= 3
        for ( param_ix_type i = 0; i < (int)_data.size(); i++ ) {
            if ( _data[i] ) {
                LOG_DEBUG3( "data[" << _minParam + i * _paramStep << "]=" << _data[i] );
            }
        }
#endif
        _data[ paramIx ] = new cache_item_type[ _maxValue - _minValue + 1 ];
        evalTable( paramIx );
    }
    return ( _data[ paramIx ] );
}

template<typename DistributionFactory, typename Param, typename Value>
void DistributionCache<DistributionFactory,Param,Value>::checkValue(
    const value_type& value
) const {
    if ( value < _minValue || value > _maxValue ) {
        THROW_EXCEPTION( std::invalid_argument, value << " is outside of supported range" );
    }
}

template<typename DistributionFactory, typename Param, typename Value>
typename DistributionCache<DistributionFactory,Param,Value>::log_prob_type
DistributionCache<DistributionFactory,Param,Value>::lnPdf(
    const param_type& param,
    const value_type& value
){
    checkValue( value );
    return ( _data( ensureTable( param ), value - _minValue ).lnPdf );
}

template<typename DistributionFactory, typename Param, typename Value>
typename DistributionCache<DistributionFactory,Param,Value>::log_prob_type
DistributionCache<DistributionFactory,Param,Value>::lnNormPdf(
    const param_type& param,
    const value_type& value
){
    checkValue( value );
    return ( _data( ensureTable( param ), value - _minValue ).lnNormPdf );
}
