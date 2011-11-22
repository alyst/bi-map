#include "math/PitmanYorProcess.h"
#include "mcmc/SplitMergeStep.h"
#include "ObjectsPartition.h"
#include "ProbesPartition.h"
#include "dynamic_bitset_utils.h"

#include "ChessboardBiclusteringGibbsSampler.h"
#include "ChessboardBiclusteringFitInternal.h"

GibbsSamplerParams::GibbsSamplerParams()
    : priorsUpdatePeriod( 50 )
    , blockResamples( 20 )
    , objectsSplitMergeParams( 11, 5 )
    , objectClusterParams( 5, 5 )
    , objectsSplitMergeRate( 0.05 )
    , objectMembershipRate( 0.02 )
    , probesSplitMergeParams( 5, 3 )
    , probeClusterParams( 5, 5 )
    , probesSplitMergeRate( 0.04 )
    , probeMembershipRate( 0.015 )
    , blockFlipRate( 0.15 )
    , objectMultipleRate( 0.07 )
    , signalRate( 0.1 )
    , meanObjectRank( 2 )
    , meanProbeRank( 2 )
{
}

ChessboardBiclusteringGibbsSampler::ChessboardBiclusteringGibbsSampler(
    const gsl_rng* rndNumGen,
    const PrecomputedData& precomputed, 
    const ChessboardBiclusteringHyperPriors& hyperpriors, 
    const GibbsSamplerParams& params, 
    const ChessboardBiclusteringEnergyEval& energyEval,
    const SamplingTransform& samplingTransform
) : _rndNumGen( rndNumGen )
  , _precomputed( precomputed )
  , _hyperpriors( hyperpriors )
  , _params( params )
  , _energyEval( energyEval )
  , _samplingTransform( samplingTransform )
  , _iteration( 0 )
{
}

ChessboardBiclusteringGibbsHelper ChessboardBiclusteringGibbsSampler::createGibbsHelper( const ChessboardBiclusteringFit& clusFit ) const
{
    return ( ChessboardBiclusteringGibbsHelper( rndNumGen(), clusFit, _energyEval, _samplingTransform ) );
}

probability_vector_t ChessboardBiclusteringGibbsSampler::elementsPickupRates(
    const ChessboardBiclusteringFit&    clus,          /** @param[in]  clustering to sample from and to */
    bool                         useObjects,    /** @param[in]  calculate objects or probes rates */
    double MaxRateRatio
) const {
    std::vector<size_t> wrongSCsum( useObjects ? clus.objectsCount() : clus.probesCount(), 0 );
    for ( object_clundex_t objCluIx = 0; objCluIx < clus.objectsClusters().size(); ++objCluIx ) {
        const object_set_t& objs = clus.objectsCluster( objCluIx ).items();
        for ( probe_clundex_t probeCluIx = 0; probeCluIx < clus.probesClusters().size(); ++probeCluIx ) {
            const probe_bitset_t& probes = clus.probesCluster( probeCluIx ).items();
            bool isCCEnabled = clus.isBlockEnabled( objCluIx, probeCluIx );
            for ( object_set_t::const_iterator oit = objs.begin(); oit != objs.end(); ++oit ) {
                foreach_bit( probe_index_t, probeIx, probes ) {
                    const assay_container_t& assays = clus.data().probe( probeIx ).assayIndexes();
                    size_t scSum = 0;
                    for ( int i = 0; i < assays.size(); i++ ) {
                        scSum += clus.data().measurement( *oit, assays[ i ] );
                    }
                    if ( scSum == 0 && isCCEnabled ) {
                        wrongSCsum[ useObjects ? *oit : probeIx ]++;
                    } else if ( scSum > 0 && !isCCEnabled ) {
                        wrongSCsum[ useObjects ? *oit : probeIx ] += scSum;
                    }
                }
            }
        }
    }
    size_t maxWrongSC = *std::max_element( wrongSCsum.begin(), wrongSCsum.end() );
    probability_vector_t res( wrongSCsum.size(), 1.0 );
    if ( maxWrongSC > 0 ) {
        for ( size_t i = 0; i < res.size(); i++ ) {
            res[ i ] = 1.0 + ( MaxRateRatio - 1.0 ) * wrongSCsum[ i ] / maxWrongSC;
        }
    }
    return ( res );
}

object_index_t ChessboardBiclusteringGibbsSampler::randomObjectToModify(
    const ChessboardBiclusteringFit&   clusFit,
    double                      MaxRateRatio
) const {
    return ( GenericDiscreteDistribution( elementsPickupRates( clusFit, true, MaxRateRatio ) )
                .random( rndNumGen() ) );
}

probe_index_t ChessboardBiclusteringGibbsSampler::randomProbeToModify(
    const ChessboardBiclusteringFit&   clusFit,
    double                      MaxRateRatio
) const {
    return ( GenericDiscreteDistribution( elementsPickupRates( clusFit, false, MaxRateRatio ) )
                .random( rndNumGen() ) );
}

/**
    Does MH(Gibbs) step in clusterings space.

    @return @todo true if clusters assignment was modified
 */
log_prob_t ChessboardBiclusteringGibbsSampler::doSamplingStep(
    ChessboardBiclusteringFit&    clus        /** @param[in,out]  clustering to sample from and to */
){
    typedef SplitMergeSamplingStep<ObjectsPartitionEx, ObjectsPartitionStats, ObjectsParamsSampler> obj_split_merge_step;
    typedef SplitMergeSamplingStep<ProbesPartitionEx, ProbesPartitionStats, ProbesParamsSampler> probe_split_merge_step;

    if ( &_precomputed != &clus.precomputed() ) {
        throw std::runtime_error( "Clustering to sample uses incompatible data" );
    }

    LOG_DEBUG2( "Step #" << _iteration << ": totalLnP=" << clus.totalLnP() );
    BOOST_ASSERT( clus.check() );
    // update clusters
    // apply one of two sampling steps
    // apply split/merge step to objects clusters
    if ( gsl_ran_bernoulli( rndNumGen(), 0.5 * _params.objectsSplitMergeRate ) ) {
        LOG_DEBUG2( "Trying objects cluster split" );
        object_clundex_t objCluIx = clus.clusterOfObject( randomObjectToModify( clus ) );
        // find object with farthest distance in neighbouring clusters
        const object_set_t& cluObjs = clus.objectsCluster( objCluIx ).items();
        object_index_t obj1ix = 0;
        PrecomputedData::dist_to_obj_t obj2Dist( OBJECT_NA, unset() );
        for ( object_set_t::const_iterator oit = cluObjs.begin(); oit != cluObjs.end(); ++oit ) {
             PrecomputedData::dist_to_obj_t dist = _precomputed.rankedObject( *oit, cluObjs, true,
                                                    gsl_ran_geometric( rndNumGen(), 1.0 / _params.meanObjectRank )-1, false, true );
             if ( !is_unset( dist.second )
                 && ( is_unset( obj2Dist.second )
                      || obj2Dist.second < dist.second )
             ){
                 obj1ix = *oit;
                 obj2Dist = dist;
             }
        }
        if ( !is_unset( obj2Dist.second ) ) {
            // do the step
            obj_split_merge_step    smStep( rndNumGen(), 
                                            ObjectsPartitionStats( rndNumGen(), clus.data(), clus.priors(),
                                                                   _energyEval.weights.objects ), 
                                            ObjectsParamsSampler( *this, true, true, true ),
                                            _samplingTransform, _params.objectsSplitMergeParams );
            obj_split_merge_step::result_type res = smStep( ObjectsPartitionEx( clus ), obj1ix, obj2Dist.first );
            if ( res.modified ) {
                clus = (const ChessboardBiclusteringFit&)res.ptn;
                clus.setObjectsClusterSamples( res.cluIx1, _params.blockResamples );
                clus.setObjectsClusterSamples( res.cluIx2, _params.blockResamples );
            }
            BOOST_ASSERT( clus.checkObjectsPartition() );
        } else if ( cluObjs.size() > 1 ) {
            THROW_RUNTIME_ERROR( "No farthest objects for cluster " << objCluIx );
        }
    }
    else if ( gsl_ran_bernoulli( rndNumGen(), 0.5 * _params.objectsSplitMergeRate ) ) {
        LOG_DEBUG2( "Trying objects cluster merge" );
        object_index_t obj1ix = randomObjectToModify( clus );
        object_clundex_t obj1CluIx = clus.clusterOfObject( obj1ix );
        PrecomputedData::dist_to_obj_t obj2Dist = _precomputed.rankedObject( obj1ix, clus.objectsCluster( obj1CluIx ).items(), false,
                                                    gsl_ran_geometric( rndNumGen(), 1.0 / _params.meanObjectRank )-1, true, true );
        if ( !is_unset( obj2Dist.second ) ) {
            // do the step
            obj_split_merge_step    smStep( rndNumGen(), 
                                            ObjectsPartitionStats( rndNumGen(), clus.data(), clus.priors(),
                                                                   _energyEval.weights.objects ), 
                                            ObjectsParamsSampler( *this, true, true, true ),
                                            _samplingTransform, _params.objectsSplitMergeParams );
            obj_split_merge_step::result_type res = smStep( ObjectsPartitionEx( clus ), obj1ix, obj2Dist.first );
            if ( res.modified ) {
                clus = (const ChessboardBiclusteringFit&)res.ptn;
                clus.setObjectsClusterSamples( res.cluIx1, _params.blockResamples );
                clus.setObjectsClusterSamples( res.cluIx2, _params.blockResamples );
            }
            BOOST_ASSERT( clus.checkObjectsPartition() );
        } else {
            LOG_DEBUG1( "No clustering partners for object " << obj1ix );
        }
    }

    // apply single element membership step with specified rate to random objects
    if ( gsl_ran_bernoulli( rndNumGen(), _params.objectMembershipRate ) ) {
        object_index_t objIx = randomObjectToModify( clus );
        LOG_DEBUG2( "Sampling object " << objIx << " membership" );
        object_clundex_t objCluIx = clus.clusterOfObject( objIx );
        size_t clusToLookup = std::min( (size_t)clus.objectsClusters().size(),
                                        (size_t)gsl_ran_geometric( rndNumGen(), 1.0 / _params.meanObjectRank ) );

        // single element membership step
        // filter clusters that are very unlikely to embrace objIx
        ObjectsClusterParams  paramsNew( (const ChessboardBiclustering&)clus );
        ObjectsClusterParams  paramsOld( (const ChessboardBiclustering&)clus );
        cluster_set_type obj2clu;
        for ( size_t i = 0; i < _precomputed.objectCoSignalDistances().size()-1 && obj2clu.size() < clusToLookup; i++ ) {
            obj2clu.insert( clus.clusterOfObject( _precomputed.objectCoSignalDistances()( objIx, i ).first ) );
        }
#if defined(_DEBUG)
        log_prob_t prevLLH = clus.llh( _energyEval.weights );
        log_prob_t prevLPP = clus.lpp();
        log_prob_t prevCluLPP = PriorEval( clus ).objectsClusteringLPP();
#endif
        GibbsSample<object_clundex_t> newCluIx = createGibbsHelper( clus ).sampleClusterOfObject( objIx,
                                                                _params.objectClusterParams,
                                                                &paramsNew, &paramsOld,
                                                                obj2clu );
        if ( newCluIx != objCluIx ) {
#if defined(_DEBUG)
            size_t oldCluSize = clus.objectsCluster( objCluIx ).size();
            size_t effCluCnt = clus.objectsClusters().size() - ( oldCluSize == 1 ? 1: 0 );
            log_prob_t cluLPPDeltaOld = log( clus.priors().objectClustering.clusterAssignmentPrior( oldCluSize - 1, effCluCnt, clus.objectsCount() - 1 ) );
 
            log_prob_t cluLPPDeltaNew;
#endif
            // update parameters only if partition really changed
            if ( newCluIx == clus.objectsClusters().size() ) {
                // new cluster
#if defined(_DEBUG)
                cluLPPDeltaNew = log( clus.priors().objectClustering.clusterAssignmentPrior( 0, effCluCnt, clus.objectsCount() - 1 ) );
#endif
                newCluIx.value = clus.addObjectCluster( objIx );
            }
            else {
#if defined(_DEBUG)
                cluLPPDeltaNew = log( clus.priors().objectClustering.clusterAssignmentPrior( clus.objectsCluster( newCluIx ).size(), effCluCnt, clus.objectsCount() - 1 ) );
#endif
                clus.setObjectCluster( objIx, newCluIx );
            }
            clus.objectsClusterParams( newCluIx ) = paramsNew;
            clus.objectsClusterParams( objCluIx ) = paramsOld;
            clus.setObjectsClusterSamples( objCluIx, _params.blockResamples );
            clus.setObjectsClusterSamples( newCluIx, _params.blockResamples );
            clus.cleanupClusters();
#if defined(_DEBUG)
            REPORT_PROB_ERROR_IF( std::abs( prevLLH + newCluIx.oddsRatio.llhRatio - clus.llh( _energyEval.weights ) ) > LN_PROB_ERROR_TOL,
                        "O LLH error: " << prevLLH << "+" << newCluIx.oddsRatio.llhRatio
                        << "=" << prevLLH + newCluIx.oddsRatio.llhRatio
                        << "!=" << clus.llh( _energyEval.weights ) );
            log_prob_t newCluLPP = PriorEval( clus ).objectsClusteringLPP();
            REPORT_PROB_ERROR_IF( std::abs( prevCluLPP + cluLPPDeltaNew - cluLPPDeltaOld - newCluLPP ) > LN_PROB_ERROR_TOL,
                              "OC LPP error: " << prevCluLPP << "+" << cluLPPDeltaNew - cluLPPDeltaOld
                              << "=" << prevCluLPP +  cluLPPDeltaNew - cluLPPDeltaOld
                              << "!=" << newCluLPP );
            log_prob_t newLPP = clus.lpp();
            REPORT_PROB_ERROR_IF( std::abs( prevLPP + newCluIx.oddsRatio.lppRatio - newLPP ) > LN_PROB_ERROR_TOL,
                             "O LPP error: " << prevLPP << "+" << newCluIx.oddsRatio.lppRatio
                             << "=" << prevLPP +  newCluIx.oddsRatio.lppRatio
                             << "!=" << newLPP );
#endif
            LOG_DEBUG2( "Object " << objIx << " put to " << newCluIx << " from " << objCluIx 
                    << " oldLnP=" << (prevLLH + prevLPP) << " newLnP=" << clus.totalLnP() );
        }
        else {
            //clus.objectsClusterParams( newCluIx ) = params;
        }
        BOOST_ASSERT( clus.checkObjectsPartition() );
    }

    if ( gsl_ran_bernoulli( rndNumGen(), 0.5 * _params.probesSplitMergeRate ) ) {
        LOG_DEBUG2( "Trying probes cluster split" );
        // find probe with farthest distance in neighbouring clusters
        probe_clundex_t probeCluIx = clus.clusterOfProbe( randomProbeToModify( clus ) );
        const probe_bitset_t& cluProbes = clus.probesCluster( probeCluIx ).items();
        probe_index_t probe1ix = 0;
        PrecomputedData::dist_to_probe_t probe2Dist( PROBE_NA, unset() );
        foreach_bit( probe_index_t, probeIx, cluProbes ) {
             PrecomputedData::dist_to_probe_t dist = _precomputed.rankedProbe( probeIx, cluProbes, true,
                                                    gsl_ran_geometric( rndNumGen(), 1.0 / _params.meanProbeRank )-1, false, true );
             if ( !is_unset( dist.second )
                 && ( is_unset( probe2Dist.second )
                      || probe2Dist.second < dist.second )
             ){
                 probe1ix = probeIx;
                 probe2Dist = dist;
             }
        }
        if ( !is_unset( probe2Dist.second ) ) {
            // do the step
            probe_split_merge_step    smStep( rndNumGen(), 
                                            ProbesPartitionStats( rndNumGen(), clus.data(), clus.priors(),
                                                                  _energyEval.weights.probes ), 
                                            ProbesParamsSampler( *this, true, true ),
                                            _samplingTransform, _params.probesSplitMergeParams );
            probe_split_merge_step::result_type res = smStep( ProbesPartitionEx( clus ), probe1ix, probe2Dist.first );
            if ( res.modified ) {
                clus = (const ChessboardBiclusteringFit&)res.ptn;
                clus.setProbesClusterSamples( res.cluIx1, _params.blockResamples );
                clus.setProbesClusterSamples( res.cluIx2, _params.blockResamples );
            }
            BOOST_ASSERT( clus.checkProbesPartition() );
        } else if ( cluProbes.count() > 1 ) {
            THROW_RUNTIME_ERROR( "No farthest probes for cluster " << probeCluIx );
        }
    }
    else if ( gsl_ran_bernoulli( rndNumGen(), 0.5 * _params.probesSplitMergeRate ) ) {
        LOG_DEBUG2( "Trying probes cluster merge" );
        probe_index_t probe1ix = randomProbeToModify( clus );
        probe_clundex_t probe1CluIx = clus.clusterOfProbe( probe1ix );
        PrecomputedData::dist_to_probe_t probe2Dist = _precomputed.rankedProbe( probe1ix, clus.probesCluster( probe1CluIx ).items(), false,
                                                    gsl_ran_geometric( rndNumGen(), 1.0 / _params.meanProbeRank )-1, true, true );
        if ( !is_unset( probe2Dist.second ) ) {
            // do the step
            probe_split_merge_step    smStep( rndNumGen(), 
                                            ProbesPartitionStats( rndNumGen(), clus.data(), clus.priors(),
                                                                  _energyEval.weights.probes ), 
                                            ProbesParamsSampler( *this, true, true ),
                                            _samplingTransform, _params.probesSplitMergeParams );
            probe_split_merge_step::result_type res = smStep( ProbesPartitionEx( clus ), probe1ix, probe2Dist.first );
            if ( res.modified ) {
                clus = (const ChessboardBiclusteringFit&)res.ptn;
                clus.setProbesClusterSamples( res.cluIx1, _params.blockResamples );
                clus.setProbesClusterSamples( res.cluIx2, _params.blockResamples );
            }
            BOOST_ASSERT( clus.checkProbesPartition() );
        } else {
            LOG_DEBUG1( "No clustering partners for probe " << probe1ix );
        }
    }

    // apply single element membership step with specified rate to random probes
    if ( gsl_ran_bernoulli( rndNumGen(), _params.probeMembershipRate ) ) {
        probe_index_t probeIx = randomProbeToModify( clus );
        LOG_DEBUG2( "Sampling probe " << probeIx << " membership" );
        probe_clundex_t probeCluIx = clus.clusterOfProbe( probeIx );
        size_t clusToLookup = std::min( (size_t)clus.probesClusters().size(),
                                        (size_t)gsl_ran_geometric( rndNumGen(), 1.0 / _params.meanProbeRank ) );

        // single element membership step
        // filter clusters that are very unlikely to embrace probeIx
        ProbesClusterParams  paramsNew( (const ChessboardBiclustering&)clus );
        ProbesClusterParams  paramsOld( (const ChessboardBiclustering&)clus );
        cluster_set_type probe2clu;
        for ( size_t i = 0; i < _precomputed.probeCoSignalDistances().size()-1 && probe2clu.size() < clusToLookup; i++ ) {
            probe2clu.insert( clus.clusterOfProbe( _precomputed.probeCoSignalDistances()( probeIx, i ).first ) );
        }
#if defined(_DEBUG)
        log_prob_t prevLLH = clus.llh( _energyEval.weights );
        log_prob_t prevLPP = clus.lpp();
        log_prob_t prevCluLPP = PriorEval( clus ).probesClusteringLPP();
#endif
        GibbsSample<probe_clundex_t> newCluIx = createGibbsHelper( clus ).sampleClusterOfProbe( probeIx,
                                                                _params.probeClusterParams,
                                                                &paramsNew, &paramsOld,
                                                                probe2clu );
        if ( newCluIx != probeCluIx ) {
            // update parameters only if partition really changed
#if defined(_DEBUG)
            size_t oldCluSize = clus.probesCluster( probeCluIx ).size();
            size_t effCluCnt = clus.probesClusters().size() - ( oldCluSize == 1 ? 1: 0 );
            log_prob_t cluLPPDeltaOld = log( clus.priors().probeClustering.clusterAssignmentPrior( oldCluSize - 1, effCluCnt, clus.probesCount() - 1 ) );
 
            log_prob_t cluLPPDeltaNew;
#endif
            if ( newCluIx == clus.probesClusters().size() ) {
                // new cluster
#if defined(_DEBUG)
                cluLPPDeltaNew = log( clus.priors().probeClustering.clusterAssignmentPrior( 0, effCluCnt, clus.probesCount() - 1 ) );
#endif
                newCluIx.value = clus.addProbeCluster( probeIx );
            }
            else {
#if defined(_DEBUG)
                cluLPPDeltaNew = log( clus.priors().probeClustering.clusterAssignmentPrior( clus.probesCluster( newCluIx ).size(), effCluCnt, clus.probesCount() - 1 ) );
#endif
                clus.setProbeCluster( probeIx, newCluIx );
            }
            clus.probesClusterParams( newCluIx ) = paramsNew;
            clus.probesClusterParams( probeCluIx ) = paramsOld;
            clus.setProbesClusterSamples( probeCluIx, _params.blockResamples );
            clus.setProbesClusterSamples( newCluIx, _params.blockResamples );
            clus.cleanupClusters();
#if defined(_DEBUG)
            REPORT_PROB_ERROR_IF( std::abs( prevLLH + newCluIx.oddsRatio.llhRatio - clus.llh(_energyEval.weights) ) > LN_PROB_ERROR_TOL,
                        "S LLH error: " << prevLLH << "+" << newCluIx.oddsRatio.llhRatio
                        << "=" << prevLLH + newCluIx.oddsRatio.llhRatio
                        << "!=" << clus.llh(_energyEval.weights) );
            log_prob_t newCluLPP = PriorEval( clus ).probesClusteringLPP();
            REPORT_PROB_ERROR_IF( std::abs( prevCluLPP + cluLPPDeltaNew - cluLPPDeltaOld - newCluLPP ) > LN_PROB_ERROR_TOL,
                              "SC LPP error: " << prevCluLPP << "+" << cluLPPDeltaNew - cluLPPDeltaOld
                              << "=" << prevCluLPP +  cluLPPDeltaNew - cluLPPDeltaOld
                              << "!="<<  newCluLPP );
            log_prob_t newLPP = clus.lpp();
            REPORT_PROB_ERROR_IF( std::abs( prevLPP + newCluIx.oddsRatio.lppRatio - newLPP ) > LN_PROB_ERROR_TOL,
                             "S LPP error: " << prevLPP << "+" << newCluIx.oddsRatio.lppRatio
                             << "=" << prevLPP +  newCluIx.oddsRatio.lppRatio
                             << "!=" << newLPP );
#endif
            LOG_DEBUG1( "Probe " << probeIx << " put to " << newCluIx << " from " << probeCluIx
                    << " oldLnP=" << (prevLLH + prevLPP) << " newLnP=" << clus.totalLnP(_energyEval.weights) );
        }
        else {
            //clus.probesClusterParams( newCluIx ) = params;
        }
        BOOST_ASSERT( clus.checkProbesPartition() );
    }
    BOOST_ASSERT( clus.checkBlocks() );

    typedef std::pair<object_clundex_t, probe_clundex_t> cc_id;

#define MAX_BLOCK_SAMPLES ((size_t)10)
#if 1
    {
        LOG_DEBUG2( "Sampling blocks probes" );
        std::vector<cc_id> ccQueue;
        // sample blocks
        for ( object_clundex_t objCluIx = 0; objCluIx < clus.objectsClusters().size(); objCluIx++ ) {
            for ( probe_clundex_t probeCluIx = 0; probeCluIx < clus.probesClusters().size(); probeCluIx++ ) {
                if ( clus.blockToSample( objCluIx, probeCluIx ) > 0 
                     || gsl_ran_bernoulli( rndNumGen(), _params.blockFlipRate )
                ){
                    ccQueue.push_back( cc_id( objCluIx, probeCluIx ) );
                }
            }
        }
        // limit number of samples -- pick N random, so the iteration doesn't take much time
        std::vector<cc_id> ccToProcess( std::min( MAX_BLOCK_SAMPLES, ccQueue.size() ) );
        gsl_ran_choose( rndNumGen(), ccToProcess.data(), ccToProcess.size(), ccQueue.data(), ccQueue.size(), sizeof( cc_id ) );
        ChessboardBiclusteringGibbsHelper helper = createGibbsHelper( clus );
        for ( size_t i = 0; i < ccToProcess.size(); i++ ) {
            const cc_id& blockId = ccToProcess[ i ];
            LOG_DEBUG2( "Sampling block (" << blockId.first << ", " << blockId.second << ")" );
            bool enableBlock = helper.sampleBlockEnablement( blockId.first, blockId.second ).value;
            bool wasEnabled = clus.isBlockEnabled( blockId.first, blockId.second );
            if ( enableBlock == wasEnabled ) continue;
            ChessboardBiclustering::block_iterator cluIt = clus.setBlock( blockId.first, blockId.second, enableBlock );
            if ( enableBlock ) {
                BOOST_ASSERT( cluIt != clus.blockNotFound() );
                GibbsSample<signal_t> newSignalSample = helper.sampleSignal( *cluIt );
                cluIt->setSignal( newSignalSample.value );
                if ( !wasEnabled ) clus.setBlockSamples( blockId.first, blockId.second, _params.blockResamples );
            }
        }
        BOOST_ASSERT( clus.checkBlocks() );
    }
#endif

    // sample object multiples
    LOG_DEBUG2( "Sampling object multiples" );
    for ( object_index_t objIx = 0; objIx < clus.objectsCount(); objIx++ ){
        if ( gsl_ran_bernoulli( rndNumGen(), _params.objectMultipleRate ) ) {
            clus.setObjectMultiple( objIx, createGibbsHelper( clus ).sampleObjectMultiple( objIx ).value );
        }
    }
    // sample signals
    {
        LOG_DEBUG2( "Sampling signals" );
        std::vector<cc_id> ccQueue;
        ChessboardBiclusteringGibbsHelper helper = createGibbsHelper( clus );
        for ( ChessboardBiclustering::block_iterator cluIt = clus.begin(); cluIt != clus.end(); ++cluIt ) {
            if ( cluIt->isEnabled() && 
                ( clus.blockToSample( cluIt->objectsClusterIndex(), cluIt->probesClusterIndex() ) > 0 
                ||  gsl_ran_bernoulli( rndNumGen(), _params.signalRate ) ) ) {
                ccQueue.push_back( cc_id( cluIt->objectsClusterIndex(), cluIt->probesClusterIndex() ) );
            }
        }
        // limit number of samples -- pick N random, so the iteration doesn't take much time
        std::vector<cc_id> ccToProcess( std::min( MAX_BLOCK_SAMPLES, ccQueue.size() ) );
        gsl_ran_choose( rndNumGen(), ccToProcess.data(), ccToProcess.size(), ccQueue.data(), ccQueue.size(), sizeof( cc_id ) );
        for ( size_t i = 0; i < ccToProcess.size(); i++ ) {
            ChessboardBiclustering::block_iterator cluIt = clus.findBlock( ccToProcess[ i ].first, ccToProcess[ i ].second );
            BOOST_ASSERT( cluIt != clus.blockNotFound() );
            LOG_DEBUG2( "Sampling signals of cluster (" 
                    << cluIt->objectsClusterIndex() << ", " << cluIt->probesClusterIndex() << ")" );
            GibbsSample<signal_t> newSignalSample = helper.sampleSignal( *cluIt );
            cluIt->setSignal( newSignalSample.value );
            clus.decBlockSamples( cluIt->objectsClusterIndex(), cluIt->probesClusterIndex() );
        }
    }

    if ( isPriorsUpdateRequired() ) {
        // update priors
#if defined(_DEBUG)
        log_prob_t oldLnP = clus.totalLnP(_energyEval.weights);
        LOG_DEBUG2( "Iteration " << _iteration << ":" );
        LOG_DEBUG2( "Signal prior: N(" << clus.baselineSignalParams().lnScRate() << "," << clus.derivedPriors().signalPrior.sigma << "^2)" );
        LOG_DEBUG2( "Noise signal: " << clus.noiseParams().successRate );
        LOG_DEBUG2( "OldLnP=" << oldLnP << " newLnP=" << clus.totalLnP(_energyEval.weights) );
#endif
        ChessboardBiclusteringGibbsHelper helper = createGibbsHelper( clus );
        clus.setSignalPrior( helper.sampleSignalPrior( _hyperpriors ) );
        clus.setNoiseParams( helper.sampleNoiseParams() );
    }

    _iteration++;

    return ( clus.totalLnP(_energyEval.weights) );
}
