#include "ChessboardBiclusteringStructureLLHEval.h"

ChessboardBiclusteringStructureLLHEval::ChessboardBiclusteringStructureLLHEval(
    const PrecomputedData& precomputed,
    prob_t priorCCEnabledProb,
    prob_t probeCluOffObjCluRate,
    prob_t objCluOffProbeCluRate,
    const ChessboardBiclustering& clustering
) : _lnPriorCCProbeRatio( log( ( 1.0 - priorCCEnabledProb ) / priorCCEnabledProb ) )
  , _probeCluPerObjClu( 1.0 - probeCluOffObjCluRate, probeCluOffObjCluRate, 1 )
  , _objCluPerProbeClu( 1.0 - objCluOffProbeCluRate, objCluOffProbeCluRate, 1 )
  , _precomputed( precomputed ), _clustering( clustering )
{
}

log_prob_t ChessboardBiclusteringStructureLLHEval::ClustersPatternMismatchLLH(
    const observations_mask_type&  mask1,
    const observations_mask_type&  mask2
){
    BOOST_ASSERT( mask1.size() == mask2.size() );

    size_t nCoOccurs = ( mask1 & mask2 ).count();
    size_t n1Occurs = mask1.count();
    size_t n2Occurs = mask2.count();

    return ( HypergeometricDistribution( n1Occurs, mask1.size() - n1Occurs, n2Occurs )
                .lnPdf( nCoOccurs ) );
}

log_prob_t ChessboardBiclusteringStructureLLHEval::objectsClusterMismatchLLH(
    const object_set_t& objects
) const {
    log_prob_t res = 0;
    for ( object_set_t::const_iterator oit = objects.begin(); oit != objects.end(); ++oit ) {
        PrecomputedData::dist_to_obj_t dist = _precomputed.rankedObject( *oit, objects, false, 0, true, false );
        if ( !is_unset( dist.second ) ) res += dist.second;
    }
    return ( res );
}

ChessboardBiclusteringStructureLLHEval::observations_mask_type ChessboardBiclusteringStructureLLHEval::compactObjectsMask(
    const ChessboardBiclusteringStructureLLHEval::observations_mask_type& objectsMask
) const {
    const observations_mask_type& spMask = _precomputed.specificObjectsMask();
    observations_mask_type res( _clustering.objectsClusters().size() );
    for ( object_clundex_t cluIx = 0; cluIx < res.size(); cluIx++ ) {
        const object_set_t& objs = _clustering.objectsCluster( cluIx ).items();
        for ( object_set_t::const_iterator oit = objs.begin(); oit != objs.end(); ++oit ) {
            if ( spMask.test( *oit ) && objectsMask.test( *oit ) ) {
                res.set( cluIx );
                break;
            }
        }
    }
    return ( res );
}

ChessboardBiclusteringStructureLLHEval::observations_mask_type ChessboardBiclusteringStructureLLHEval::compactProbesMask(
    const ChessboardBiclusteringStructureLLHEval::observations_mask_type& probesMask
) const {
    const observations_mask_type& spMask = _precomputed.specificProbesMask();
    observations_mask_type res( _clustering.probesClusters().size() );
    for ( probe_clundex_t cluIx = 0; cluIx < res.size(); cluIx++ ) {
        const probe_bitset_t& probes = _clustering.probesCluster( cluIx ).items();
        if ( ( probes & spMask & probesMask ).any() ) {
            res.set( cluIx );
        }
    }
    return ( res );
}

log_prob_t ChessboardBiclusteringStructureLLHEval::probesClusterMismatchLLH(
    const probe_bitset_t& probes
) const {
    log_prob_t res = 0;
    foreach_bit( probe_index_t, probeIx, probes ) {
        PrecomputedData::dist_to_probe_t dist = _precomputed.rankedProbe( probeIx, probes, false, 0, true, false );
        if ( !is_unset( dist.second ) ) res += dist.second;
    }
    return ( res );
}
