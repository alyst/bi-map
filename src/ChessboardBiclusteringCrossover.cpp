#include "ObjectsPartition.h"
#include "ProbesPartition.h"
#include "mcmc/PartitionCrossover.h"

#include "ChessboardBiclusteringCrossover.h"

#define OBJ_XOVER_RATE 0.7
#define SPLIT_XOVER_RATE 0.5
#define XOVER_CLUSTER_SAMPLES 5

ChessboardBiclusteringCrossoverGenerator::ChessboardBiclusteringCrossoverGenerator(
    const gsl_rng* rndNumGen,
    const PrecomputedData& precomputed,
    const ChessboardBiclusteringPriors& priors,
    const ChessboardBiclusteringHyperPriors& hyperpriors
)   : _precomputed( precomputed )
    , _priors( priors )
    , _sampler( rndNumGen, precomputed, hyperpriors, GibbsSamplerParams() )
{

}

void ChessboardBiclusteringCrossoverGenerator::push_to_result(
    particle_container_type&    res,
    const ChessboardBiclusteringFit&   ptn
) const {
    const_cast<ChessboardBiclusteringFit&>( ptn ).cleanupClusters();
    //const_cast<ChessboardBiclusteringFit&>( ptn ).setAllBlockSamples( 2 );
    BOOST_ASSERT( ptn.check() );
    res.push_back( StaticChessboardBiclustering( ptn ) );
}

ChessboardBiclusteringCrossoverGenerator::particle_container_type 
ChessboardBiclusteringCrossoverGenerator::operator()(
    const gsl_rng*                 rng,
    const particle_cache_type&     cache
) const {
    particle_container_type res;
    if ( cache.size() < 2 )                   return ( res );
    size_t ix1 = gsl_rng_uniform_int( rng, cache.size() );
    size_t ix2 = gsl_rng_uniform_int( rng, cache.size()-1 );
    if ( ix2 == ix1 ) ix2++;
    ChessboardBiclusteringFit ptn1( _precomputed, _priors, *cache.nthParticle( ix1 )->particle.clustering );
    ChessboardBiclusteringFit ptn2( _precomputed, _priors, *cache.nthParticle( ix2 )->particle.clustering );

    if ( gsl_rng_uniform( rng ) < OBJ_XOVER_RATE ) {
        // object's partition crossover
        ObjectsPartitionEx      objPtn1( ptn1 );
        ObjectsPartitionEx      objPtn2( ptn2 );
        ObjectsParamsSampler    paramsSampler( _sampler, true, true, true );

        if ( gsl_rng_uniform( rng ) < SPLIT_XOVER_RATE ) {
            std::pair<object_clundex_t, object_clundex_t> xover =
                split_cluster_by_partition( rng, paramsSampler, objPtn1, objPtn2 );
            if ( xover.first != ObjectsPartitionEx::ClusterNA ) {
                objPtn1.setClusterSamples( xover.first, XOVER_CLUSTER_SAMPLES );
                objPtn1.setClusterSamples( xover.second, XOVER_CLUSTER_SAMPLES );
                push_to_result( res, objPtn1 );
            }
        }
        else {
            std::pair< std::pair<object_clundex_t, object_clundex_t>,
                       std::pair<object_clundex_t, object_clundex_t> > xover =
                       exchange_cluster_parts( rng, paramsSampler, objPtn1, objPtn2 );
            if ( xover.first.first != ObjectsPartitionEx::ClusterNA ) {
                objPtn1.setClusterSamples( xover.first.first, XOVER_CLUSTER_SAMPLES );
                objPtn1.setClusterSamples( xover.first.second, XOVER_CLUSTER_SAMPLES );
                objPtn2.setClusterSamples( xover.second.first, XOVER_CLUSTER_SAMPLES );
                objPtn2.setClusterSamples( xover.second.second, XOVER_CLUSTER_SAMPLES );
                push_to_result( res, objPtn1 );
                push_to_result( res, objPtn2 );
            }
        }
    }
    else {
        // probes's partition crossover
        ProbesPartitionEx probePtn1( ptn1 );
        ProbesPartitionEx probePtn2( ptn2 );
        ProbesParamsSampler    paramsSampler( _sampler, true, true );

        if ( gsl_rng_uniform( rng ) < SPLIT_XOVER_RATE ) {
            std::pair<probe_clundex_t, probe_clundex_t> xover =
                split_cluster_by_partition( rng, paramsSampler, probePtn1, probePtn2 );
            if ( xover.first != ProbesPartitionEx::ClusterNA ) {
                probePtn1.setClusterSamples( xover.first, XOVER_CLUSTER_SAMPLES  );
                probePtn1.setClusterSamples( xover.second, XOVER_CLUSTER_SAMPLES  );
                push_to_result( res, probePtn1 );
            }
        }
        else {
            std::pair< std::pair<probe_clundex_t, probe_clundex_t>,
                       std::pair<probe_clundex_t, probe_clundex_t> > xover =
                       exchange_cluster_parts( rng, paramsSampler, probePtn1, probePtn2 );
            if ( xover.first.first != ProbesPartitionEx::ClusterNA ) {
                probePtn1.setClusterSamples( xover.first.first, XOVER_CLUSTER_SAMPLES );
                probePtn1.setClusterSamples( xover.first.second, XOVER_CLUSTER_SAMPLES );
                probePtn2.setClusterSamples( xover.second.first, XOVER_CLUSTER_SAMPLES );
                probePtn2.setClusterSamples( xover.second.second, XOVER_CLUSTER_SAMPLES );
                push_to_result( res, probePtn1 );
                push_to_result( res, probePtn2 );
            }
        }
    }
    return ( res );
}
