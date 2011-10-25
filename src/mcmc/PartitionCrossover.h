/**
    Partitions crossover
 */
#pragma once

#include "../BasicTypedefs.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>

#include "../math/Permutation.h"

/**
 *  Finds a pair of intersecting (but not identical clusters),
 *  one from each partition.
 */
template<class Partition>
bool find_intersecting_pair(
    const gsl_rng*      rng,
    const Partition&    ptn1,   /** i first partition */
    typename Partition::cluster_index_type& clu1Ix,     /** o first cluster */
    typename Partition::cluster_index_type exclClu1Ix,  /** i cluster of ptn1 excluded from consideration */
    const Partition&    ptn2,   /** i second partition */
    typename Partition::cluster_index_type& clu2Ix,     /** o second cluster */
    typename Partition::cluster_index_type exclClu2Ix   /** i cluster of ptn2 excluded from consideration */
){
    typedef typename Partition::cluster_proxy_type clu_t;
    typedef typename Partition::cluster_index_type cluix_t;
    typedef typename Partition::elements_set_proxy_type elmset_t;

    Permutation ptn1Pm = Permutation::Random( rng, ptn1.clustersCount() );
    for ( size_t i = 0; i < ptn1Pm.size(); i++ ) {
        clu1Ix = ptn1Pm[ i ];
        if ( clu1Ix == exclClu1Ix ) continue;
        clu_t clu = ptn1.cluster( clu1Ix );
        Permutation elmPm = Permutation::Random( rng, clu.size() );
        for ( size_t j = 0; j < elmPm.size(); j++ ) {
            size_t elmIx = clu.items()[ elmPm[ j ] ];
            clu2Ix = ptn2.clusterIndex( elmIx );
            if ( clu2Ix == exclClu2Ix ) continue;
            if ( clu.items() != ptn2.cluster( clu2Ix ).items() ) {
                return ( true );
            }
        }
    }
    return ( false );
}

/**
 *  Fills missing parameters of the cluster.
 */
template<class Partition, class ParamsSampler>
void fill_missing_params(
    const gsl_rng*          rng,
    const ParamsSampler&    paramsSampler,
    Partition&              ptn,    /** io  partition to divide */
    typename Partition::cluster_index_type cluIx    /** i cluster to fill the params */
){
    typedef typename Partition::cluster_proxy_type clu_t;
    typedef typename Partition::cluster_params_type params_t;

    clu_t clu = ptn.cluster( cluIx );
    params_t params = clu.params();
    paramsSampler( ptn, params, clu.items(), clu.items(), false, false );
    ptn.cluster( cluIx ).params() = params;
}

/**
    Splits random cluster of the partition by clusters of another.

    Based on 'Evolutionary Monte Carlo Methods for Clustering' by G.Goswami.
 */
template<class Partition, class ParamsSampler>
std::pair<typename Partition::cluster_index_type, typename Partition::cluster_index_type>
split_cluster_by_partition(
    const gsl_rng*      rng,
    const ParamsSampler& paramsSampler,
    Partition&          ptn,    /** io  partition to divide */
    const Partition&    divisor /** i   divisor partition */
){
    typedef typename Partition::cluster_proxy_type clu_t;
    typedef typename Partition::cluster_index_type cluix_t;
    typedef typename Partition::elements_set_proxy_type elmset_t;

    cluix_t dCluIx = Partition::ClusterNA;
    cluix_t clu1Ix = Partition::ClusterNA;
    cluix_t clu2Ix = Partition::ClusterNA;
    // cycle over clusters of divisor, to find one, which is different from each cluster of ptn
    Permutation divPm = Permutation::Random( rng, divisor.clustersCount() );
    for ( size_t i = 0; i < divPm.size(); i++ ) {
        dCluIx = divPm[ i ];
        clu_t dClu = divisor.cluster( dCluIx );
        clu1Ix = Partition::ClusterNA;
        clu2Ix = Partition::ClusterNA;
        // cycle over elements of dClu to find 2 distinct ptn clusters, intersecting with it
        Permutation elmPm = Permutation::Random( rng, dClu.size() );
        for ( size_t j = 0; j < elmPm.size(); j++ ) {
            size_t elmIx = dClu.items()[ elmPm[ j ] ];
            cluix_t cluIx = ptn.clusterIndex( elmIx );
            if ( clu1Ix == Partition::ClusterNA ) {
                clu1Ix = cluIx;
            }
            else if ( clu1Ix != cluIx && clu2Ix == Partition::ClusterNA ) {
                clu2Ix = ptn.clusterIndex( elmIx );
                break; // break, we found both
            }
        }
        if ( clu2Ix != Partition::ClusterNA ) break;
    }
    if ( clu2Ix == Partition::ClusterNA ) {
        // partitions identical, no division
        return ( std::make_pair( Partition::ClusterNA, Partition::ClusterNA ) );
    }
    // create one of new clusters
    elmset_t elms1 = ptn.cluster( clu1Ix ).items();
    elmset_t elms2 = ptn.cluster( clu2Ix ).items();
    elms1 &= divisor.cluster( dCluIx ).items();
    elms2 -= divisor.cluster( dCluIx ).items();
    elms2 += elms1;
    ptn.exchangeElements( clu1Ix, clu2Ix, elms2 );
    // now, clu1' = (clu1 \ dClu) | (clu2 & dClu)
    // clu2' = (clu1 & dClu) | (clu2 \ dClu)

    // fill missing params
    fill_missing_params( rng, paramsSampler, ptn, clu1Ix );
    fill_missing_params( rng, paramsSampler, ptn, clu2Ix );

    return ( std::make_pair( clu1Ix, clu2Ix ) );
}

/**
    Exchanges elements between a pair of clusters from two partitions.
    The clusters in the pairs are intersecting with corresponding cluster of the other pair.
 */
template<class Partition, class ParamsSampler>
std::pair<typename std::pair<typename Partition::cluster_index_type, typename Partition::cluster_index_type>,
          typename std::pair<typename Partition::cluster_index_type, typename Partition::cluster_index_type> >
exchange_cluster_parts(
    const gsl_rng*      rng,
    const ParamsSampler& paramsSampler,
    Partition&          ptn1,    /** i fisrt partition */
    Partition&          ptn2     /** i second partition */
){
    typedef typename Partition::cluster_proxy_type clu_t;
    typedef typename Partition::cluster_index_type cluix_t;
    typedef typename Partition::elements_set_proxy_type elmset_t;
    typedef typename Partition::cluster_params_type params_t;

    cluix_t ptn1Clu1Ix = Partition::ClusterNA;
    cluix_t ptn1Clu2Ix = Partition::ClusterNA;
    cluix_t ptn2Clu1Ix = Partition::ClusterNA;
    cluix_t ptn2Clu2Ix = Partition::ClusterNA;

    // find 2 distinct pairs of intersecting (but not identical) pairs of clusters from both partitions
    if ( !find_intersecting_pair( rng, ptn1, ptn1Clu1Ix, Partition::ClusterNA, ptn2, ptn2Clu1Ix, Partition::ClusterNA ) 
      || !find_intersecting_pair( rng, ptn1, ptn1Clu2Ix, ptn1Clu1Ix, ptn2, ptn2Clu2Ix, ptn2Clu1Ix )
    ){
        return ( std::make_pair( std::make_pair( Partition::ClusterNA, Partition::ClusterNA ), 
                                 std::make_pair( Partition::ClusterNA, Partition::ClusterNA ) ) );
    }

    // isect1 = clu11 & clu21, isect2 = clu12 & clu22
    elmset_t isect1 = ptn1.cluster( ptn1Clu1Ix ).items();
    elmset_t isect2 = ptn1.cluster( ptn1Clu2Ix ).items();
    isect1 &= ptn2.cluster( ptn2Clu1Ix ).items();
    isect2 &= ptn2.cluster( ptn2Clu2Ix ).items();

    // modify first partition
    elmset_t newElms12 = ptn1.cluster( ptn1Clu2Ix ).items();
    newElms12 -= isect2;
    newElms12 += isect1;
    ptn1.exchangeElements( ptn1Clu1Ix, ptn1Clu2Ix, newElms12 );
    // clu12' = (clu11 & clu21) + (clu12 - clu22)
    // clu11' = (clu11 - clu21) + (clu12 & clu22)
    fill_missing_params( rng, paramsSampler, ptn1, ptn1Clu1Ix );
    fill_missing_params( rng, paramsSampler, ptn1, ptn1Clu2Ix );

    // modify second partition
    elmset_t newElms22 = ptn2.cluster( ptn2Clu2Ix ).items();
    newElms22 -= isect2;
    newElms22 += isect1;
    ptn2.exchangeElements( ptn2Clu1Ix, ptn2Clu2Ix, newElms12 );
    // clu22' = (clu21 & clu11) + (clu22 - clu12)
    // clu21' = (clu21 - clu11) + (clu22 & clu12)
    fill_missing_params( rng, paramsSampler, ptn2, ptn2Clu1Ix );
    fill_missing_params( rng, paramsSampler, ptn2, ptn2Clu2Ix );

    return ( std::make_pair( std::make_pair( ptn1Clu1Ix, ptn1Clu2Ix ), 
                             std::make_pair( ptn2Clu1Ix, ptn2Clu2Ix ) ) );
}
