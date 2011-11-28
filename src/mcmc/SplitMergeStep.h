/**
    Split-merge sampling step template
 */
#pragma once

#include "../BasicTypedefs.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

#include <boost/serialization/serialization.hpp>

#include "../math/GenericDiscreteDistribution.h"
#include "../math/Permutation.h"

#include "MetropolisHastingsStep.h"

/**
 *  Params of split-merge sampling step.
 * 
 *  @see SplitMergeSamplingStep
 */
struct SplitMergeStepParams {
    size_t  splitLaunchSamples; // how many gibbs sampling iterations to do for split launch probe
    size_t  mergeLaunchSamples; // how many gibbs sampling iterations to do for merge launch probe

    SplitMergeStepParams( size_t splitLaunchSamples = 8, size_t mergeLaunchSamples = 4 )
    : splitLaunchSamples( splitLaunchSamples )
    , mergeLaunchSamples( mergeLaunchSamples )
    {}

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( splitLaunchSamples );
        ar & BOOST_SERIALIZATION_NVP( mergeLaunchSamples );
    }
};

/**
    Non-conjugate split-merge sampling step.

    @see
    S.Jain, R.Neal, "Splitting and Merging Components of a Nonconjugate DPMM",
    Bayesian Analysis 2007, N.3
 */
template<class Partition, class PartitionStats, class ParamsSampler>
class SplitMergeSamplingStep
{
public:
    typedef Partition partition_type;

    struct Result {
        typedef typename Partition::cluster_index_type cluster_index_type;
        partition_type      ptn;
        bool                modified;
        cluster_index_type  cluIx1;
        cluster_index_type  cluIx2;

        Result( const partition_type& ptn,
                bool modified,
                cluster_index_type cluIx1,
                cluster_index_type cluIx2 )
        : ptn( ptn )
        , modified( modified )
        , cluIx1( cluIx1 )
        , cluIx2( cluIx2 )
        {}
    };

    typedef Result result_type;

private:
    typedef typename Partition::element_index_type elmix_t;
    typedef typename Partition::cluster_index_type cluix_t;
    typedef typename Partition::cluster_proxy_type clupxy_t;
    typedef typename Partition::elements_set_proxy_type elmix_set_t;
    typedef typename ParamsSampler::params_type cluster_params_type;
    typedef std::pair<elmix_t, log_prob_t> sample_cluster_result;

    struct InternalResult {
        partition_type  ptn;
        log_prob_t      lnTransKernel;  /** log of transition probability */

        InternalResult( const partition_type& launchPtn )
        : ptn( launchPtn ), lnTransKernel( 0.0 )
        {}
    };

    typedef InternalResult internal_result_type;

    const gsl_rng*              rng;
    const PartitionStats&       stats;
    const ParamsSampler&        paramsSampler;
    const SplitMergeStepParams  stepParams;
    const SamplingTransform     samplingTransform;

    log_prob_t clusterAssignmentLOR( const partition_type& ptn, elmix_t elmIx, cluix_t cluIx, cluix_t cluAltIx ) const;
    log_prob_t clusterParamsTransitionLP( const partition_type& ptnLaunch, const partition_type& ptnFinal, cluix_t cluIx ) const;

    sample_cluster_result sampleClusterIndex( const partition_type& ptn, elmix_t elmIx, cluix_t clu1Ix, cluix_t clu2Ix ) const;
    log_prob_t sampleClusterParams( partition_type& ptn, cluix_t cluIx, bool posterior, bool overwrite = true ) const;

    partition_type splitClusterLaunchProbe( const partition_type& ptn, elmix_t elmIx1, elmix_t elmIx2 ) const;
    internal_result_type splitCluster( const partition_type& ptn, elmix_t elmIx1, elmix_t elmIx2 ) const;
    internal_result_type mergeClusters( const partition_type& ptn, elmix_t elmIx1, elmix_t elmIx2 ) const;
    log_prob_t evalTotalLP( const partition_type& ptn ) const;

public:
    SplitMergeSamplingStep( const gsl_rng* rng,
                            const PartitionStats& stats,
                            const ParamsSampler& paramsSampler,
                            const SamplingTransform& samplingTransform = SamplingTransform(),
                            const SplitMergeStepParams& stepParams = SplitMergeStepParams() )
    : rng( rng )
    , stats( stats )
    , paramsSampler( paramsSampler )
    , stepParams( stepParams )
    , samplingTransform( samplingTransform )
    {};

    result_type operator()( const partition_type& ptn ) const;
    result_type operator()( const partition_type& ptn, elmix_t elmIx1, elmix_t elmIx2 ) const;
};

/**
 *  Generate launch (initial) partition for split move.
 *
 *  @return partition (launch probe) with split cluster
 */
template<class Partition, class PartitionStats, class ParamsSampler>
typename SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::partition_type 
SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::splitClusterLaunchProbe(
    const partition_type&   ptn,    /** original partition */
    elmix_t                 elm1Ix, /** element to be present in first split cluster of launch probe */
    elmix_t                 elm2Ix  /** element to be present in second split cluster of launch probe */
) const {
    BOOST_ASSERT( elm1Ix != elm2Ix );
    cluix_t clu1Ix = ptn.clusterIndex( elm1Ix );
    cluix_t clu2Ix = ptn.clusterIndex( elm2Ix );
    elmix_set_t cluElems = ptn.cluster( clu1Ix ).items();
    if ( clu1Ix != clu2Ix ) {
        cluElems += ptn.cluster( clu2Ix ).items();
    }
    BOOST_ASSERT( cluElems.contains( elm1Ix ) );
    BOOST_ASSERT( cluElems.contains( elm2Ix ) );

    cluix_t cluNewIx;
    partition_type launchPtn = ptn;

    // initial random assignment of element between clusters
    {
        LOG_DEBUG2( "splitClusterLaunchProbe(): initial assignment" );
        elmix_set_t newClu2Elems( launchPtn );

        for ( size_t i = 0; i < cluElems.size(); i++ ) {
            elmix_t elmIx = cluElems[ i ];
            if ( elmIx == elm2Ix || ( elmIx != elm1Ix && gsl_rng_uniform( rng ) < 0.5 ) ) {
                newClu2Elems += elmIx;
            }
        }
        BOOST_ASSERT( newClu2Elems.contains( elm2Ix ) );
        BOOST_ASSERT( !newClu2Elems.contains( elm1Ix ) );
        cluNewIx = launchPtn.exchangeElements( clu1Ix, clu1Ix == clu2Ix ? Partition::ClusterNA : clu2Ix, newClu2Elems );
        LOG_DEBUG2( "splitClusterLaunchProbe(): initial params sample for cluster #" << clu1Ix );
        sampleClusterParams( launchPtn, clu1Ix, false );
        LOG_DEBUG2( "splitClusterLaunchProbe(): initial params sample for cluster #" << cluNewIx );
        sampleClusterParams( launchPtn, cluNewIx, false );
    }

    // intermediate restricted Gibbs sampling scans
    for ( size_t i = 0; i < stepParams.splitLaunchSamples; i++ ) {
        LOG_DEBUG2( "splitClusterLaunchProbe(): intermediate scan #" << i );

        Permutation perm = Permutation::Random( rng, cluElems.size() );
        for ( size_t j = 0; j < perm.size(); j++ ) {
            elmix_t elmIx = cluElems[ perm[ j ] ];
            // update element-to-cluster assignment
            if ( elmIx != elm1Ix && elmIx != elm2Ix ) {
                cluix_t prevCluIx = launchPtn.clusterIndex( elmIx );
                cluix_t newCluIx = sampleClusterIndex( launchPtn, elmIx, prevCluIx, prevCluIx == clu1Ix ? cluNewIx : clu1Ix ).first;
                if ( newCluIx != prevCluIx ) {
                    launchPtn.putToCluster( elmIx, newCluIx );
                    LOG_DEBUG2( "splitClusterLaunchProbe():   fill missing params to " << newCluIx );
                    sampleClusterParams( launchPtn, newCluIx, false, false ); // fill-in missing params for posterior samples
                }
            }
        }
        LOG_DEBUG2( "splitClusterLaunchProbe():   Gibbs sampling in #" << clu1Ix );
        sampleClusterParams( launchPtn, clu1Ix, true );
        LOG_DEBUG2( "splitClusterLaunchProbe():   Gibbs sampling in #" << cluNewIx );
        sampleClusterParams( launchPtn, cluNewIx, true );
    }
    LOG_DEBUG3( "splitClusterLaunchProbe(): done" );
    return ( launchPtn );
}

/**
 *  Calculate necessary statistics for splitting cluster in two,
 *  so that 2 selected elements of original cluster
 *  are contained in separate new clusters.
 */
template<class Partition, class PartitionStats, class ParamsSampler>
typename SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::internal_result_type 
SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::splitCluster(
    const partition_type&   ptn,
    elmix_t                 elmIx1,
    elmix_t                 elmIx2
) const {
    BOOST_ASSERT( elmIx1 != elmIx2 );
    cluix_t cluIx = ptn.clusterIndex( elmIx1 );
    BOOST_ASSERT( cluIx == ptn.clusterIndex( elmIx2 ) );

    // generate launch probe
    clupxy_t clu = ptn.cluster( cluIx );
    elmix_set_t cluElems = clu.items();
    LOG_DEBUG2( "splitCluster(): split launch probe for cluster #" << cluIx );
    partition_type launchPtn = splitClusterLaunchProbe( ptn, elmIx1, elmIx2 );
    BOOST_ASSERT( cluIx == launchPtn.clusterIndex( elmIx1 ) );
    cluix_t cluNewIx = launchPtn.clusterIndex( elmIx2 );
    BOOST_ASSERT( cluIx != cluNewIx );

    // final Gibbs sampling scan, each element is scanned once
    LOG_DEBUG2( "splitCluster(): final scan" );
    internal_result_type res = internal_result_type( launchPtn );
    Permutation perm = Permutation::Random( rng, cluElems.size() );
    for ( size_t j = 0; j < perm.size(); j++ ) {
        elmix_t elmIx = cluElems[ perm[ j ] ];
        // update element-to-cluster assignment
        if ( elmIx != elmIx1 && elmIx != elmIx2 ) {
            cluix_t prevCluIx = launchPtn.clusterIndex( elmIx );
            sample_cluster_result elmCluRes = sampleClusterIndex( launchPtn, elmIx, prevCluIx, prevCluIx == cluIx ? cluNewIx : cluIx );
            // update transition ratio from launch probe cluster to new
            res.lnTransKernel -= ln_odds_to_prob( elmCluRes.second );
            res.ptn.putToCluster( elmIx, elmCluRes.first );
            // make sure model has complete valid set of parameters, but don't overwrite the existing ones
            sampleClusterParams( res.ptn, elmCluRes.first, false, false );
        }
    }
    // sample new cluster params once
    launchPtn = res.ptn; // remember launch partition for parameters comparison
    LOG_DEBUG2( "splitCluster(): update params #" << cluIx );
    sampleClusterParams( res.ptn, cluIx, true );
    LOG_DEBUG2( "splitCluster(): update params #" << cluNewIx );
    sampleClusterParams( res.ptn, cluNewIx, true );
    res.lnTransKernel -= clusterParamsTransitionLP( launchPtn, res.ptn, cluIx )
                       + clusterParamsTransitionLP( launchPtn, res.ptn, cluNewIx );

    // intermediate restricted Gibbs sampling scans rollback transition ratio calculation
    LOG_DEBUG2( "splitCluster(): rollback launch probe scans" );
    partition_type rollbackLaunchPtn( ptn );
    for ( size_t i = 1; i < stepParams.mergeLaunchSamples; i++ ) {
        sampleClusterParams( rollbackLaunchPtn, cluIx, true );
    }

    // calculate transition ratio of cluster params for split and rollback-merge partitions
    res.lnTransKernel += clusterParamsTransitionLP( rollbackLaunchPtn, ptn, cluIx );

    return ( res );
}

/**
 *  Calculate statistics for merging two cluster different clusters,
 *  each represented by one element,
 *  into single one.
 */
template<class Partition, class PartitionStats, class ParamsSampler>
typename SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::internal_result_type 
SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::mergeClusters(
    const partition_type&   ptn,
    elmix_t                 elmIx1,
    elmix_t                 elmIx2
) const {
    // elements should be from different clusters
    cluix_t clu1Ix = ptn.clusterIndex( elmIx1 );
    cluix_t clu2Ix = ptn.clusterIndex( elmIx2 );
    BOOST_ASSERT( clu1Ix != clu2Ix );

    // calculate launch partition for merging
    partition_type launchPtn( ptn );
    LOG_DEBUG2( "mergeClusters(): merging " << clu1Ix << " and " << clu2Ix );
    elmix_set_t noClu2Elems = ptn.cluster( clu2Ix ).items();
    noClu2Elems -= noClu2Elems;
    launchPtn.exchangeElements( clu1Ix, clu2Ix, noClu2Elems );
    LOG_DEBUG2( "mergeClusters(): Initial sampling");
    sampleClusterParams( launchPtn, clu1Ix, false );

    LOG_DEBUG2( "mergeClusters(): intermediate scans");
    // intermediate restricted Gibbs sampling scans for launch probe
    for ( size_t i = 1; i < stepParams.mergeLaunchSamples; i++ ) {
        sampleClusterParams( launchPtn, clu1Ix, true );
    }
    // final scan
    internal_result_type res = internal_result_type( launchPtn );
    LOG_DEBUG2( "mergeClusters(): final scan " << clu1Ix << " and " << clu2Ix );
    sampleClusterParams( res.ptn, clu1Ix, true );
    // parameters transition probability
    res.lnTransKernel -= clusterParamsTransitionLP( launchPtn, res.ptn, clu1Ix );

    {
        LOG_DEBUG2( "mergeClusters(): split launch probe" );
        // intermediate restricted Gibbs sampling scans for rollback transition ratio calculation
        partition_type rollbackLaunchPtn = splitClusterLaunchProbe( ptn, elmIx1, elmIx2 );
        // calculate both clusters parameters transition probability
        LOG_DEBUG2( "mergeClusters(): transition kernel for rollback ptn params" );
        res.lnTransKernel += clusterParamsTransitionLP( rollbackLaunchPtn, ptn, clu1Ix )
                           + clusterParamsTransitionLP( rollbackLaunchPtn, ptn, clu2Ix );

        // calculate transition ratio log
        // cluster assignment for elements of first cluster
        LOG_DEBUG2( "mergeClusters(): transition kernel for rollback ptn partition" );
        elmix_set_t clu1Elems = ptn.cluster( clu1Ix ).items();
        Permutation perm1 = Permutation::Random( rng, clu1Elems.size() );
        for ( size_t j = 0; j < perm1.size(); j++ ) {
            elmix_t elmIx = clu1Elems[ perm1[ j ] ];
            cluix_t rbkLCluIx = rollbackLaunchPtn.clusterIndex( elmIx );
            if ( rbkLCluIx != clu1Ix ) {
                res.lnTransKernel += ln_odds_to_prob( clusterAssignmentLOR( rollbackLaunchPtn, elmIx, rbkLCluIx, clu1Ix ) );
                rollbackLaunchPtn.putToCluster( elmIx, clu1Ix );
                // make sure model has complete valid set of parameters, but don't overwrite the existing ones
                sampleClusterParams( rollbackLaunchPtn, clu1Ix, false, false );
            }
        }
        elmix_set_t clu2Elems = ptn.cluster( clu2Ix ).items();
        Permutation perm2 = Permutation::Random( rng, clu2Elems.size() );
        for ( size_t j = 0; j < perm2.size(); j++ ) {
            elmix_t elmIx = clu2Elems[ perm2[ j ] ];
            cluix_t rbkLCluIx = rollbackLaunchPtn.clusterIndex( elmIx );
            if ( rbkLCluIx != clu2Ix ) {
                res.lnTransKernel += ln_odds_to_prob( clusterAssignmentLOR( rollbackLaunchPtn, elmIx, rbkLCluIx, clu2Ix ) );
                rollbackLaunchPtn.putToCluster( elmIx, clu2Ix );
                // make sure model has complete valid set of parameters, but don't overwrite the existing ones
                sampleClusterParams( rollbackLaunchPtn, clu2Ix, false, false );
            }
        }
    }

    return ( res );
}

/**
    Log of odds ratio of assigning element to cluAltIx,
    intsead of cluIx cluster.
 */
template<class Partition, class PartitionStats, class ParamsSampler>
log_prob_t SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::clusterAssignmentLOR( 
    const partition_type&   ptn,
    elmix_t                 elmIx, 
    cluix_t                 cluIx,      /** i index of cluster element is assigned to */
    cluix_t                 cluAltIx    /** i index of alternative cluster */
) const {
    BOOST_ASSERT( ptn.clusterIndex( elmIx ) == cluIx );
    typedef typename Partition::cluster_index_set_type cluix_set_type;
    log_prob_t cluLpp = stats.clusterAssignmentLPP( (int)cluIx >= 0 && cluIx < (cluix_t)ptn.clustersCount() 
                                                    ? ptn.cluster( cluIx ).size() - 1 : 0, 
                                                    ptn.clustersCount(), ptn.elementsCount() );
    log_prob_t cluAltLpp = stats.clusterAssignmentLPP( (int)cluAltIx >= 0 && cluAltIx < (cluix_t)ptn.clustersCount() 
                                                       ? ptn.cluster( cluAltIx ).size(): 0, 
                                                       ptn.clustersCount(), ptn.elementsCount() );
    std::vector<elmix_set_t> newClusters;
    newClusters.push_back( ptn.cluster( cluIx ).items() );
    newClusters[0] -= elmIx;
    newClusters.push_back( ptn.cluster( cluAltIx ).items() );
    newClusters[1] += elmIx;
    std::vector<cluster_params_type> newParams;
    newParams.push_back( ptn.cluster( cluIx ).params() );
    newParams.push_back( ptn.cluster( cluAltIx ).params() );
    cluix_set_type  cluIxs;
    cluIxs.insert( cluIx );
    cluIxs.insert( cluAltIx );
    // make sure parameters for elmIx exist
    paramsSampler( ptn, newParams[0], newClusters[0], elmix_set_t( ptn, elmIx ), false, false );
    paramsSampler( ptn, newParams[1], newClusters[1], elmix_set_t( ptn, elmIx ), false, false );
    log_prob_t llhDelta = stats.llhDelta( ptn, newClusters, newParams, cluIxs );

    return ( cluAltLpp - cluLpp + llhDelta );
}

/**
    Sample element assignment to alternative cluster.
    @return a new cluster index and probability of selected alternative.
 */
template<class Partition, class PartitionStats, class ParamsSampler>
typename SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::sample_cluster_result
SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::sampleClusterIndex(
    const partition_type&   ptn, 
    elmix_t                 elmIx,
    cluix_t                 cluIx,      /** i index of cluster element is currently assigned to */
    cluix_t                 cluAltIx    /** i index of alternative cluster */
) const {
    BOOST_ASSERT( ptn.clusterIndex( elmIx ) == cluIx );
    log_prob_t lnOddRatio = clusterAssignmentLOR( ptn, elmIx, cluIx, cluAltIx );
    log_prob_t adjLogRatio = samplingTransform.clampedDeltaLnP( lnOddRatio, 
                                                                stats.lpp( ptn ) + (log_prob_t)stats.llh( ptn ) )
                           / samplingTransform.temperature;
    sample_cluster_result res;
    // a / ( a + b ) = 1 / ( 1 + b/a )
    double p = 1.0 / ( 1.0 + exp( -adjLogRatio ) );
    res.first = gsl_rng_uniform( rng ) < p ? cluAltIx : cluIx;
    res.second = res.first == cluAltIx ? lnOddRatio : -lnOddRatio;
    return ( res );
}

/**
    Sample parameters of the cluster.
    @return log of conditional probability of the new parameters, 
    given all other parameters of the cluster fixed.
 */
template<class Partition, class PartitionStats, class ParamsSampler>
log_prob_t SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::sampleClusterParams(
    partition_type&         ptn,        /** i  partition */
    cluix_t                 cluIx,      /** i  index of cluster in the partition */
    bool                    posterior,  /** i  if true, posterior sampling, using current values in params, otherwise, generate params using priors */
    bool                    overwrite   /** i  if true, overwrite currently set params, if false, skip their modification */
) const {
    LOG_DEBUG2( "sampleClusterParams( #" << cluIx 
            << ", " << ( posterior ? "posterior" : "prior" )
            << ", " << ( overwrite ? "overwrite" : "no-overwrite" ) << " )" );
    clupxy_t clu = ptn.cluster( cluIx );
    cluster_params_type params = clu.params();
    double logProb = paramsSampler( ptn, params, clu.items(), clu.items(), overwrite, posterior );
    clu.params() = params;
    return ( logProb );
}

/**
    Do split-merge step for two randomly selected elements of the partition.
 */
template<class Partition, class PartitionStats, class ParamsSampler>
typename SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::result_type 
SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::operator()(
    const Partition&    ptn
) const {
    if ( ptn.elementsCount() <= 1 ) {
        return ( result_type( ptn, false, 0, 0 ) );
    }
    elmix_t elmIx1 = gsl_rng_uniform_int( rng, ptn.elementsCount() );
    elmix_t elmIx2 = gsl_rng_uniform_int( rng, ptn.elementsCount() - 1 );
    if ( elmIx2 >= elmIx1 ) elmIx2++;

    return ( operator()( ptn, elmIx1, elmIx2 ) );
}

/**
    Do split-merge step for cluster(s), which contain(s) given elements.
 */
template<class Partition, class PartitionStats, class ParamsSampler>
typename SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::result_type 
SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::operator()(
    const Partition&    ptn,
    elmix_t             elmIx1,
    elmix_t             elmIx2
) const {
    BOOST_ASSERT( elmIx1 != elmIx2 );
    log_prob_t oldTotalLP = evalTotalLP( ptn );
    if ( boost::math::isnan( oldTotalLP ) ) THROW_RUNTIME_ERROR( "Bad partition energy" );
    cluix_t cluIx1 = ptn.clusterIndex( elmIx1 );
    cluix_t cluIx2 = ptn.clusterIndex( elmIx2 );
    if ( cluIx1 > cluIx2 ) {
        std::swap( cluIx1, cluIx2 );
        std::swap( elmIx1, elmIx2 );
    }
    LOG_DEBUG2_IF( cluIx1 == cluIx2, "Proposed splitting of cluster " << cluIx1 << "(" << ptn.cluster( cluIx1 ).label() << ")" );
    LOG_DEBUG2_IF( cluIx1 != cluIx2, "Proposed merging of clusters " << cluIx1 << "(" << ptn.cluster( cluIx1 ).label() << ") and " << cluIx2 << "(" << ptn.cluster( cluIx2 ).label() << ")" );

    internal_result_type res = ( cluIx1 != cluIx2 ) 
                             ? mergeClusters( ptn, elmIx1, elmIx2 ) 
                             : splitCluster( ptn, elmIx1, elmIx2 );
    log_prob_t newTotalLP = evalTotalLP( res.ptn );
    if ( boost::math::isnan( newTotalLP ) ) THROW_RUNTIME_ERROR( "Bad new partition energy" );

    LOG_DEBUG1_IF( cluIx1 == cluIx2, "Proposed splitting of cluster " << ptn.cluster( cluIx1 ).label() <<
                                     " into " << res.ptn.cluster( cluIx1 ).label() << "(" << res.ptn.cluster( cluIx1 ).size() 
        << ") and " << res.ptn.cluster( res.ptn.clustersCount() - 1 ).label() << "(" << res.ptn.cluster( res.ptn.clustersCount() - 1 ).size() <<  ")" );
    LOG_DEBUG1_IF( cluIx1 != cluIx2, "Proposed merging of clusters " << ptn.cluster( cluIx1 ).label()
                                     << " and " << ptn.cluster( cluIx2 ).label()
                                     << " into " << res.ptn.cluster( cluIx1 ).label()
                                     << "(" << res.ptn.cluster( cluIx1 ).size() <<  ")" );
    LOG_DEBUG1( "Old LnP=" << oldTotalLP << " new LnP=" << newTotalLP << " lnTransKernel=" << res.lnTransKernel );
    BOOST_ASSERT( !boost::math::isnan( res.lnTransKernel ) );
    log_prob_t adjLogRatio = samplingTransform.clampedDeltaLnP( newTotalLP - oldTotalLP, oldTotalLP );
    if ( is_accept_step_by_mh( rng, exp( ( adjLogRatio + res.lnTransKernel ) / samplingTransform.temperature ), 1 ) ) {
        LOG_DEBUG1( ( cluIx1 != cluIx2 ? "Merging" : "Splitting" ) << " successfull" );
        return ( result_type( res.ptn, true, 
                              res.ptn.clusterIndex( elmIx1 ), 
                              res.ptn.clusterIndex( elmIx2 ) ) );
    }
    else {
        LOG_DEBUG2( ( cluIx1 != cluIx2 ? "Merging" : "Splitting" ) << " unsuccessfull" );
        return ( result_type( ptn, false, cluIx1, cluIx2 ) );
    }
}

/**
    Calculate cluster params transition probability.
 */
template<class Partition, class PartitionStats, class ParamsSampler>
log_prob_t SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::clusterParamsTransitionLP(
    const partition_type&   ptnLaunch, 
    const partition_type&   ptnFinal,
    cluix_t                 cluIx       /** cluster, that has different parameters in two partitions */
) const {
    return ( paramsSampler.transitionLP( ptnLaunch, ptnFinal, cluIx ) );
}

/**
 *  Evalulate logarithm of total partition probability.
 */
template<class Partition, class PartitionStats, class ParamsSampler>
log_prob_t SplitMergeSamplingStep<Partition, PartitionStats, ParamsSampler>::evalTotalLP(
    const partition_type&   ptn
) const {
    return ( (log_prob_t)stats.llh( ptn ) + stats.lpp( ptn ) );
}

