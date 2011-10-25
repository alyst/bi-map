#include "PrecomputedData.h"

#include <ostream>
#include <boost/format.hpp>
#include <gsl/gsl_min.h>

#include "dynamic_bitset_utils.h"

ENABLE_STATIC_TRACKING( PrecomputedDataParams )
ENABLE_STATIC_TRACKING( PrecomputedData )

PrecomputedDataParams::PrecomputedDataParams()
: objectFreqThreshold( 0.65 )
, probeFreqThreshold( 0.65 )
, objectNeighbourThreshold( 0.05 )
, probeNeighbourThreshold( 0.05 )
{
}

PrecomputedData::PrecomputedData(
    const OPAData&                  data,
    const PrecomputedDataParams&    params,
    const CellSignalParams&         signalParams
) : _data( data )
  , _scDistribCache( signalParams, 0, _data.maxMeasurement().sc, signalParams.scRateStep )
  , _llhMaximizer( data, signalParams, _scDistribCache, true )
  , _objectCoHitsDistances( ObjectsObservations( data ), params.probeFreqThreshold, true )
  , _probeCoHitsDistances( ProbesObservations( data ), params.objectFreqThreshold, true )
  , _objectCoSignalDistances( _llhMaximizer.evalObjectCoSignalDistances() )
  , _probeCoSignalDistances( _llhMaximizer.evalProbeCoSignalDistances() )
  , _objectsNeighbourhood( CreateCoOccurrenceGraph( _objectCoHitsDistances.allObservations(),
                                                    _objectCoHitsDistances.specificObservations(),
                                                     params.objectNeighbourThreshold ) )
  , _probesNeighbourhood( CreateCoOccurrenceGraph( _probeCoHitsDistances.allObservations(),
                                                   _probeCoHitsDistances.specificObservations(),
                                                   params.probeNeighbourThreshold ) )
{
    LOG_INFO( "Objects sampling bias ln_odds_ratio=" << boost::format( "%.3f" ) % _probeCoHitsDistances.observationsBias() );
    LOG_INFO( "Probes sampling bias ln_odds_ratio=" << boost::format( "%.3f" ) % _objectCoHitsDistances.observationsBias() );
}

/**
 *  Serialization-only version of ctor.
 */
PrecomputedData::PrecomputedData(
    const OPAData&                  data,
    const CellSignalParams&         signalParams
) : _data( data )
  , _scDistribCache( signalParams, 0, _data.maxMeasurement().sc, signalParams.scRateStep )
  , _llhMaximizer( data, signalParams, _scDistribCache, false )
{
}

/**
 *  Bitmask objects x probes, containing 1 if object
 *  is detected in any assay of given probe.
 */
std::vector<observations_mask_t> PrecomputedData::ObjectsObservations( const OPAData& data )
{
    std::vector<observations_mask_t> res;
    res.resize( data.objectsCount(), observations_mask_t( data.probesCount() ) );
    for ( object_index_t objIx = 0; objIx < data.objectsCount(); objIx++ ) {
        for ( probe_index_t probeIx = 0; probeIx < data.probesCount(); probeIx++ ) {
            const assay_container_t& assayIxs = data.probe( probeIx ).assayIndexes();
            for ( size_t aix = 0; aix < assayIxs.size(); aix++ ) {
                if ( data.measurement( objIx, assayIxs[ aix ] ).sc > 0 ) {
                    res[ objIx ].set( probeIx );
                    break;
                }
            }
        }
    }
    return ( res );
}

/**
 *  Bitmask objects x probes, containing 1 if object
 *  is detected in any assay of given probe.
 */
std::vector<observations_mask_t> PrecomputedData::ProbesObservations( const OPAData& data )
{
    std::vector<observations_mask_t> res;
    res.resize( data.probesCount(), observations_mask_t( data.objectsCount() ) );
    for ( object_index_t objIx = 0; objIx < data.objectsCount(); objIx++ ) {
        for ( probe_index_t probeIx = 0; probeIx < data.probesCount(); probeIx++ ) {
            const assay_container_t& assayIxs = data.probe( probeIx ).assayIndexes();
            for ( size_t aix = 0; aix < assayIxs.size(); aix++ ) {
                if ( data.measurement( objIx, assayIxs[ aix ] ).sc > 0 ) {
                    res[ probeIx ].set( objIx );
                    break;
                }
            }
        }
    }
    return ( res );
}

PrecomputedData::observations_mask_t PrecomputedData::objectsUnion(
    const probe_bitset_t& probes
) const {
    observations_mask_t res( _data.objectsCount() );
    foreach_bit( probe_index_t, probeIx, probes ) {
        res |= _probeCoHitsDistances.observations( probeIx );
    }
    return ( res );
}

PrecomputedData::observations_mask_t PrecomputedData::probesUnion(
    const object_set_t& objects
) const {
    observations_mask_t res( _data.probesCount() );
    for ( object_set_t::const_iterator oit = objects.begin(); oit != objects.end(); ++oit ) {
        res |= _objectCoHitsDistances.observations( *oit );
    }
    return ( res );
}

PrecomputedData::dist_to_obj_t PrecomputedData::nearestExternalObjectHitsMatch(
    const object_set_t& objects
) const {
    dist_to_obj_t res( OBJECT_NA, 0 );
    for ( object_set_t::const_iterator oit = objects.begin(); oit != objects.end(); ++oit ) {
        object_index_t objIx = *oit;
        // find minimal distance for given object
        for ( size_t i = 0; i < _objectCoHitsDistances.size(); i++ ) {
            const dist_to_obj_t& d2e = _objectCoHitsDistances( objIx, i );
            if ( d2e.second < res.second && objects.find( d2e.first ) == objects.end() ) {
                res = d2e;
            }
            if ( d2e.second >= res.second ) {
                // there would be no nearer object for objIx
                break;
            }
        }
    }
    BOOST_ASSERT( !is_unset( res.second ) );
    return ( res );
}

PrecomputedData::dist_to_probe_t PrecomputedData::nearestExternalProbeHitsMatch(
    const probe_bitset_t& probes
) const {
    dist_to_probe_t res( PROBE_NA, 0 );
    foreach_bit( probe_index_t, probeIx, probes ) {
        // find minimal distance for given object
        for ( size_t i = 0; i < _probeCoHitsDistances.size(); i++ ) {
            const dist_to_probe_t& d2e = _probeCoHitsDistances( probeIx, i );
            if ( d2e.second < res.second && !probes.test( d2e.first ) ) {
                res = d2e;
            }
            if ( d2e.second >= res.second ) {
                // there would be no nearer obj for objIx
                break;
            }
        }
    }
    BOOST_ASSERT( !is_unset( res.second ) );
    return ( res );
}
