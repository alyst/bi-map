#pragma once

#include "../BasicTypedefs.h"

#include <vector>
#include <gsl/gsl_rng.h>

#include "../math/GenericDiscreteDistribution.h"
#include "../math/logmath.h"

#include "GibbsSample.h"
#include "SamlpingTransform.h"

/**
 *  Parameters of cluster-of-element gibbs sampling step.
 * 
 *  @see SampleClusterOfElement()
 */
struct ClusterOfElementStepParams
{
    size_t  existingClusterSamples; /** how many parameters samples for existing clusters to generate */
    size_t  newClusterSamples;      /** how many parameters samples for new cluster to generate */

    ClusterOfElementStepParams(
        size_t existingClusterSamples = 3,
        size_t newClusterSamples = 3
    ) : existingClusterSamples( existingClusterSamples )
      , newClusterSamples( newClusterSamples )
    {}

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( existingClusterSamples );
        ar & BOOST_SERIALIZATION_NVP( newClusterSamples );
    }
};

/**
 *  Internal sample of new cluster of the element.
 *  @see SampleClusterOfElement()
 */
template<typename Partition>
struct InternalClusterOfElementSample {
    typedef typename Partition::cluster_params_type params_type;
    typedef typename Partition::cluster_index_type cluix_type;
    typedef typename Partition::element_index_type elmix_type;

    cluix_type      index;          /// index of the new cluster
    OddsRatio       oddsRatio;      /// delta in probability log
    params_type     paramsNew;      /// new parameters for the new cluster
    params_type     paramsOld;      /// new parameters for the old cluster

    InternalClusterOfElementSample(
        cluix_type          index,
        log_prob_t          lppDelta,
        log_prob_t          llhDelta,
        const params_type&  paramsNew,
        const params_type&  paramsOld
    ) : index( index ), oddsRatio( llhDelta, lppDelta )
      , paramsNew( paramsNew ), paramsOld( paramsOld )
    {}
};

/**
    Single element membership step.

    @see Radford M.Neal, "Markov Chain Sampling Methods for Dirichlet Process Mixture Models", JCGS, Vol.9, No.2

    @return index of cluster to put element to,
            if = current number of clusters -- put to new cluster 
 */
template<class Partition, class PartitionStats, class ParamsSampler>
GibbsSample<typename Partition::cluster_index_type> SampleClusterOfElement(
        const gsl_rng*                  rng,                /** i random number generator */
        const Partition&                partition,          /** i partition of elements */
        const PartitionStats&           stats,              /** i partition/cluster/object likelihood and prior probability calculator */
        const typename Partition::element_index_type    elementIx,  /** i element to use */
        const ParamsSampler&            paramsSampler,      /** i cluster and element parameters sampler */
        const ClusterOfElementStepParams&   stepParams,     /** i step parameters */
        const SamplingTransform&        transform = SamplingTransform(),
        log_prob_t                      ptnTotalLnP = unset(), /** i partition total LnP */
        typename ParamsSampler::params_type*    paramsNew = NULL, /** [o] new parameters of element's new cluster */
        typename ParamsSampler::params_type*    paramsOld = NULL, /** [o] new parameters of element's old cluster */
        typename Partition::cluster_index_set_type  cluAlternatives = typename Partition::cluster_index_set_type() /** [i] alternative clusters, if empty, all clusters used */
){
    typedef typename ParamsSampler::params_type params_type;
    typedef typename Partition::cluster_index_type cluix_type;
    typedef typename Partition::element_index_type elmix_type;
    typedef typename Partition::elements_set_proxy_type elmix_set_type;
    typedef typename Partition::cluster_index_set_type cluix_set_type;
    typedef typename Partition::cluster_proxy_type clupxy_type;
    typedef InternalClusterOfElementSample<Partition> clusample_type;

    clupxy_type curClu = partition.cluster( partition.clusterIndex( elementIx ) );
    bool   curClusterSingleton = curClu.size() == 1;
    size_t existingClustersCount = partition.clustersCount();

    // clustering alternatives
    std::vector<clusample_type>      cluSamples;
    if ( curClusterSingleton ) existingClustersCount--;
    cluSamples.reserve( stepParams.existingClusterSamples * existingClustersCount + stepParams.newClusterSamples );

    elmix_set_type singleElm( partition, elementIx );
    elmix_set_type noElm( singleElm );
    elmix_set_type curCluWithoutElem = curClu.items();
    curCluWithoutElem -= singleElm;
    noElm -= singleElm;

    params_type currentParams = curClu.params( singleElm );
    params_type curCluWithoutElemParams = paramsSampler.defaults();

    log_prob_t curCluAssignmentLPP = stats.clusterAssignmentLPP( 
                                        curCluWithoutElem.size(),
                                        existingClustersCount,
                                        partition.elementsCount() )
                                   + stats.paramsLPP( partition, partition.cluster( curClu.index() ).params(),
                                                      curClu.items(), curClu.items() );
    if ( !curClusterSingleton ) {
        // do Gibbs steps to update the current parameters to the cluster without elmIx
        curCluWithoutElemParams = curClu.params( curCluWithoutElem );
        paramsSampler( curCluWithoutElemParams, curCluWithoutElem, noElm, false, false );
        for ( size_t i = 1; i < stepParams.existingClusterSamples; ++i ) {
            paramsSampler( curCluWithoutElemParams, curCluWithoutElem, noElm, true, true );
        }
        curCluAssignmentLPP -= stats.paramsLPP( partition, curCluWithoutElemParams,
                                                curCluWithoutElem, curCluWithoutElem );
    }

    // push current model (newParams already contains currentParams)
    cluSamples.push_back( clusample_type( curClu.index(), 0, 0, currentParams, currentParams ) );

    // generate probabilities for existing clusters
    if ( cluAlternatives.empty() ) {
        // use all clusters by default
        for ( cluix_type cluIx = 0; cluIx < (cluix_type)partition.clustersCount(); cluIx++ ) {
            cluAlternatives.insert( cluIx );
        }
    }
    cluAlternatives.erase( curClu.index() ); // exclude current cluster
    std::vector<params_type>        newParams( curClusterSingleton ? 1 : 2, curCluWithoutElemParams );
    std::vector<elmix_set_type>     newClustersItems( curClusterSingleton ? 1 : 2, curCluWithoutElem );

    for ( typename cluix_set_type::const_iterator cluIxIt = cluAlternatives.begin(); 
         cluIxIt != cluAlternatives.end(); ++cluIxIt
    ){
        cluix_type newCluIx = *cluIxIt;
        clupxy_type newClu = partition.cluster( newCluIx );
        BOOST_ASSERT( !newClu.items().contains( elementIx ) );
        size_t newCluSize = newClu.size();

        newParams[0] = newClu.params( singleElm );
        newClustersItems[0] = newClu.items();
        newClustersItems[0] += elementIx;
        paramsSampler( newParams[0], newClustersItems[0], singleElm, false, false );

        // calculate LPP of modified partition
        log_prob_t clusterPriorDelta =
                stats.clusterAssignmentLPP( newCluSize, existingClustersCount, partition.elementsCount() )
                - stats.paramsLPP( partition, newClu.params(), newClu.items(), newClu.items() )
                - curCluAssignmentLPP;
        BOOST_ASSERT( !is_unset( clusterPriorDelta ) );

        // calculate LLH of modified partition
        cluix_set_type oldClusters;
        oldClusters.insert( newCluIx );
        oldClusters.insert( curClu.index() );

        // generate samples
        //elmix_set_type allElms = curClu.items();
        for ( size_t i = 0; i < stepParams.existingClusterSamples; ++i ) {
            paramsSampler( newParams[0], newClustersItems[0], singleElm, true, true );
            cluSamples.push_back( clusample_type( newCluIx, clusterPriorDelta
                                   + stats.paramsLPP( partition, newParams[0],
                                                      newClustersItems[0], newClustersItems[0] ),
                                   stats.llhDelta( partition, newClustersItems, newParams, oldClusters ),
                                   newParams[0], curCluWithoutElemParams ) );
        }
    }

    // generate newClusterClones(-1) trials for the new cluster parameters
    if ( !curClusterSingleton && stepParams.newClusterSamples > 0 ) {
        newClustersItems[0] = singleElm;
        cluix_set_type oldClusters;
        oldClusters.insert( curClu.index() );

        cluix_type newCluIx = partition.clustersCount();
        log_prob_t newClusterPriorDelta =
                stats.clusterAssignmentLPP( 0, existingClustersCount, partition.elementsCount() )
                - curCluAssignmentLPP;
        newParams[0] = currentParams;
        BOOST_ASSERT( !is_unset( newClusterPriorDelta ) );
        for ( size_t i = 0; i < stepParams.newClusterSamples; ++i ) {
            paramsSampler( newParams[0], noElm, singleElm, true, true );
            cluSamples.push_back( clusample_type( newCluIx, newClusterPriorDelta
                                                  + stats.paramsLPP( partition, newParams[0], noElm, singleElm ),
                                                  stats.llhDelta( partition, newClustersItems, newParams, oldClusters ),
                                                  newParams[0], curCluWithoutElemParams
                                                ) );
        }
    }

    // apply clamping and temperature
    log_prob_t maxCluLnP = -std::numeric_limits<log_prob_t>::infinity();
    if ( is_unset( ptnTotalLnP ) ) ptnTotalLnP = stats.lpp( partition ) + stats.llh( partition );
    probability_vector_t    clusampleProbs( cluSamples.size(), unset() );
    for ( size_t i = 0; i < cluSamples.size(); i++ ) {
        log_prob_t& prob = clusampleProbs[i];
        prob = transform( ptnTotalLnP + cluSamples[i].oddsRatio.lppRatio 
                          + cluSamples[i].oddsRatio.llhRatio );
        // adjust probs due to multiple samples
        if ( cluSamples[i].index != curClu.index() ) {
            prob -= log_int( cluSamples[i].index < (cluix_type)partition.clustersCount()
                             ? stepParams.existingClusterSamples
                             : stepParams.newClusterSamples );
        }
        if ( maxCluLnP < prob ) maxCluLnP = prob;
    }

    // transform logarithms of probabilities to probabilities (unnormalized)
    for ( size_t i = 0; i < cluSamples.size(); i++ ) {
        log_prob_t& prob = clusampleProbs[i];
        BOOST_ASSERT( !is_unset( prob ) );
        prob = is_finite( maxCluLnP )
             ? exp( prob - maxCluLnP )
             : ( maxCluLnP > prob ? 0.0 : 1.0 );
        BOOST_ASSERT( !is_unset( prob ) );
    }
#if DEBUG_LEVEL >= 2
    LOG_DEBUG2( "SEM(" << elementIx << "):" );
    for ( size_t i = 0; i < cluSamples.size(); ++i ) {
        const clusample_type& cluSample = cluSamples[i];
        if ( i > 0 ) std::cerr << " ";
        std::cerr
            << ( cluSample.index < (cluix_type)partition.clustersCount() ? cluSample.index : -1 ) 
            << "=" << ( cluSample.oddsRatio.lppRatio + cluSample.oddsRatio.llhRatio )
            << "=" << clusampleProbs[i];
    }
    std::cerr << " of " << partition.clustersCount();
    std::cerr << std::endl;
#endif
    GenericDiscreteDistribution  clustersDistrib( clusampleProbs );
    const clusample_type& selSample = cluSamples[ clustersDistrib.random( rng ) ];
    if ( selSample.index != curClu.index() ) {
        if ( paramsNew ) *paramsNew = selSample.paramsNew;
        if ( paramsOld ) *paramsOld = selSample.paramsOld;
    }
    LOG_DEBUG2_IF( selSample.index == partition.clustersCount(), "assigned " << elementIx << " to new cluster" ); 

    return ( GibbsSample<cluix_type>( selSample.index, selSample.oddsRatio ) );
}
