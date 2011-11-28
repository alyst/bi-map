#include <limits>

#include "dynamic_bitset_utils.h"

#include "ChessboardBiclusteringFitInternal.h"

#include "ProbesPartition.h"

std::string ProbesPartition::ProbesClusterProxy::label() const
{
    std::ostringstream out;
    out << "S";
    foreach_bit( probe_index_t, probeIx, cluster().items() ) {
        out << (int)probeIx;
    }
    return ( out.str() );
}

void ProbesPartitionEx::putToCluster(
    element_index_type  elmIx,
    cluster_index_type  cluIx
){
    unshare();
    wrapped().setProbeCluster( elmIx, cluIx );
}

/**
 *  Exchange given set of elements between 2 clusters.
 *  If second cluster index is PROBE_NA, new cluster is created.
 */
ProbesPartitionEx::cluster_index_type ProbesPartitionEx::exchangeElements(
    cluster_index_type clu1Ix, /** i first cluster */
    cluster_index_type clu2Ix, /** i second cluster, might be missing */
    const element_index_set_type& clu2Elements /** i elements of updated clu2 cluster */
){
    unshare();
    cluster_index_type res = wrapped().exchangeProbes( clu1Ix, clu2Ix, clu2Elements );
    wrapped().cleanupClusters();
    return ( res );
}

LLHMetrics FixedProbesPartitionStats::probesLLH(
    const probe_bitset_t&   probes,
    const params_type&      params
) const {
    LLHMetrics llh;

    // quantitative component
    llh.quant = 0;
    ChessboardBiclusteringLLHEval eval = LLHEval( clusFit );
    // enabled cells
    foreach_bit( object_clundex_t, objCluIx, params.blocksMask ) {
        const object_set_t&   objects = clusFit.objectsCluster( objCluIx ).items();
        signal_t signal = params.objectsSignal[ objCluIx ];
        BOOST_ASSERT( !is_unset( signal ) );
        llh.quant += eval.cellsDataLLH( objects, probes, signal );
        BOOST_ASSERT( is_finite( llh ) );
    }

    // disabled cells
    boost::dynamic_bitset<> notMask = ~params.blocksMask;
    const DataSignalNoiseCache& snCache = clusFit.signalNoiseCache();
    foreach_bit( object_clundex_t, objCluIx, notMask ) {
        llh.quant += snCache.noiseLLH( clusFit.objectsCluster( objCluIx ).items(), probes );
        BOOST_ASSERT( !is_unset( llh ) );
    }

    ChessboardBiclusteringStructureLLHEval structEval = StructureLLHEval( clusFit );
    llh.topo = structEval.probesClusterMismatchLLH( probes );
    llh.conf = structEval.objectClustersPerProbeClusterLLH( clusFit.boundObjectsClusters( probes ).size() );

    return ( llh );
}

LLHMetrics FixedProbesPartitionStats::llhDelta(
    const std::vector<ProbesPartition::elements_set_proxy_type>& newClusters,    /** new clusters */
    const std::vector<params_type>&                 newParams,      /** params of new clusters */
    const boost::unordered_set<object_clundex_t>&   oldIndexes      /** indicies of clusters,
                                                                        whose objects are in newClusters */
) const {
    typedef std::set<object_clundex_t> objclu_set_t;

    BOOST_ASSERT( newClusters.size() == newParams.size() );

    // calculate llh
    LLHMetrics llh( 0 );
    for ( size_t i = 0; i < newClusters.size(); i++ ) {
        // new cluster's internal LLH
        llh += probesLLH( newClusters[i],  newParams[i] );
    }
    // subtract llh of old clusters
    objclu_set_t boundObjClus;
    for ( ProbesPartition::cluster_index_set_type::const_iterator oldCluIt = oldIndexes.begin();
          oldCluIt != oldIndexes.end(); ++oldCluIt
    ){
        // subtract block fitting LLH of old clusters
        cluster_index_t oldCluIx = *oldCluIt;
        llh -= clusFit.probesClusterLLH( oldCluIx );
        objclu_set_t objClus = clusFit.boundObjectsClusters( oldCluIx );
        boundObjClus.insert( objClus.begin(), objClus.end() );
    }
    // update delta due to changed probe clusters bound to object cluster
    ChessboardBiclusteringStructureLLHEval structEval = StructureLLHEval( clusFit );
    for ( objclu_set_t::const_iterator ocit = boundObjClus.begin(); ocit != boundObjClus.end(); ++ocit ) {
        std::set<probe_clundex_t> boundProbeClusNew;
        const object_set_t& objs = clusFit.objectsCluster( *ocit ).items();
        for ( object_set_t::const_iterator oit = objs.begin(); oit != objs.end(); ++oit ) {
            std::pair<OPAData::const_bait_to_probe_iterator, OPAData::const_bait_to_probe_iterator> probesRange = clusFit.data().probesOfBait( *oit );
            for ( OPAData::const_bait_to_probe_iterator stit = probesRange.first; stit != probesRange.second; ++stit ) {
                probe_index_t probeIx = stit->second;
                // try to search probe in new clusters
                bool found = false;
                for ( probe_clundex_t i = 0; i < newClusters.size(); i++ ) {
                    if ( newClusters[i].contains( probeIx ) ) {
                        boundProbeClusNew.insert( clusFit.probesClusters().size() + i );
                        found = true;
                        break;
                    }
                }
                // use unmodified clusters
                if ( !found ) boundProbeClusNew.insert( clusFit.clusterOfProbe( probeIx ) );
            }
        }
        llh.conf += structEval.probeClustersPerObjectClusterLLH( boundProbeClusNew.size() )
                 - structEval.probeClustersPerObjectClusterLLH( clusFit.boundProbesClusters( *ocit ).size() );
    }

    return ( llh );
}

log_prob_t FixedProbesPartitionStats::paramsLPP(
    const params_type&  params,
    const probe_bitset_t& clusterProbes,
    const probe_bitset_t& sampledProbes 
) const {
    log_prob_t  lpp = 0;
    ChessboardBiclusteringPriorEval eval( clusFit.data(), clusFit.priors(), clusFit );
    for ( probe_clundex_t i = 0; i < params.blocksMask.size(); i++ ) {
        lpp += params.blocksMask.test( i )
             ? eval.blockEnablementPrior( params.objectsSignal[i] )( true )
             : eval.blockEnablementPrior()( false );
    }
    return ( lpp );
}

bool FixedProbesParamsSampler::operator()(
    params_type&            params,
    const probe_bitset_t&   clusterProbes,
    const probe_bitset_t&   sampledProbes,
    bool                    overwrite,
    bool                    posterior
) const {
    bool randomized = false;
    probe_bitset_t allProbes = clusterProbes;
    allProbes |= sampledProbes;
    if ( allProbes.none() ) {
        THROW_EXCEPTION( std::invalid_argument, "Probes set for sampling is empty" );
    }

    for ( object_clundex_t objCluIx = 0; objCluIx < clusHelper->clustering().objectsClusters().size(); ++objCluIx ) {
        const ObjectsCluster&    objsClu = clusHelper->clustering().objectsCluster( objCluIx );

        bool isEnabled = params.blocksMask.test( objCluIx );
        if ( overwrite && sampleBlockMask ) {
            isEnabled = posterior
                      ? clusHelper->sampleBlockEnablement( objsClu.items(), allProbes, isEnabled ).value
                      : PriorEval( clusHelper->clusteringFit() )
                        .blockEnablementPrior()
                        .generate( clusHelper->rndNumGen() );
            randomized = true;
        }
        params.blocksMask.set( objCluIx, isEnabled );

        if ( isEnabled && sampleSignals ) {
            // sample signals, if required and some signals already present
            // if no signals, generate missing
            randomized = true;
            signal_t signal = params.objectsSignal[ objCluIx ];
            LOG_DEBUG2( "Generating random signal for objects cluster " << objCluIx );
            if ( overwrite && posterior && !is_unset( signal ) ) {
                params.objectsSignal[ objCluIx ] = clusHelper->sampleSignal( objsClu.items(), allProbes, 
                                                                            signal ).value;
            }
            else if ( is_unset( signal ) || overwrite ) {
                params.objectsSignal[ objCluIx ] = clusHelper->initialSignal( objsClu.items(), allProbes,
                                                                              NULL ).value;
            }
        }
    }
    return ( randomized );
}

bool ProbesParamsSampler::operator()(
    const ProbesPartition&  ptn, 
    params_type&            params, 
    const probe_bitset_t&   clusterProbes, 
    const probe_bitset_t&   sampledProbes, 
    bool                    overwrite, 
    bool                    posterior
) const {
    ChessboardBiclusteringGibbsHelper helper = clusSampler->createGibbsHelper( ptn );
    return ( FixedProbesParamsSampler( helper, sampleBlockMask, sampleSignals )
                                     ( params, clusterProbes, sampledProbes, 
                                       overwrite, posterior ) );
}

log_prob_t ProbesParamsSampler::transitionLP(
    const ProbesPartition&  before,
    const ProbesPartition&  after,
    probe_clundex_t         cluIx
) const {
    ChessboardBiclusteringGibbsHelper helperBefore = clusSampler->createGibbsHelper( before );
    ChessboardBiclusteringGibbsHelper helperAfter = clusSampler->createGibbsHelper( after );
    ChessboardBiclusteringPriorEval priorEval = PriorEval( helperBefore.clusteringFit() );
    log_prob_t  lnProb = 0;

    for ( object_clundex_t objCluIx = 0; objCluIx < helperBefore.clustering().objectsClusters().size(); ++objCluIx ) {
        ChessboardBiclustering::block_proxy cluBefore = *((const ChessboardBiclustering&)before).findBlock( objCluIx, cluIx, false );
        ChessboardBiclustering::block_proxy cluAfter = *((const ChessboardBiclustering&)after).findBlock( objCluIx, cluIx, false );
        bool enabledBefore = cluBefore.isEnabled();
        bool enabledAfter = cluAfter.isEnabled();
        if ( enabledBefore != enabledAfter ) {
            // block probe transition probability
            ChessboardBiclusteringGibbsHelper::BlockEnablementDataLLHCached llhBefore( helperBefore.energyEval(), helperBefore.clusteringFit(), objCluIx, cluIx );
            ChessboardBiclusteringGibbsHelper::BlockEnablementDataLLHCached llhAfter( helperAfter.energyEval(), helperAfter.clusteringFit(), objCluIx, cluIx );
            BernoulliDistribution prior = priorEval.blockEnablementPrior();
            log_prob_t  lnOddsRatio = ( prior( enabledAfter ) - prior( enabledBefore ) )
                                    + ( llhAfter( enabledAfter ) - llhBefore( enabledBefore ) );
            lnProb += ln_odds_to_prob( lnOddsRatio );
        }
        else if ( enabledBefore ) {
            const signal_t beforeSignal = cluBefore.signal();
            const signal_t afterSignal = cluAfter.signal();
            // signals transition probability (both clusters enabled)
            const GaussianDistribution& signalPrior = priorEval.signalPrior();
            log_prob_t lnOddsRatio = helperAfter.clusteringFit().blockLLH( objCluIx, cluIx )
                                   - helperBefore.clusteringFit().blockLLH( objCluIx, cluIx );
            lnProb += ln_odds_to_prob( lnOddsRatio )
                    + signalPrior( afterSignal ) - signalPrior( beforeSignal );
        }
        BOOST_ASSERT( !is_unset( lnProb ) );
    }
    return ( lnProb );
}
