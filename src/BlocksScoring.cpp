#include "cemm/bimap/BlocksScoring.h"

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>

#include "cemm/bimap/IndexedPartitionsCollection.h"

namespace cemm { namespace bimap {

typedef boost::accumulators::accumulator_set<
        signal_t, boost::accumulators::stats<
            boost::accumulators::tag::weighted_mean>
        , prob_t> signal_mean_accum_t;
struct BlockAccums {
    prob_t count_total;     /// weighted total number of object cluster x probe cluster blocks
    prob_t count_on;        /// weighted number of blocks in on-state
    signal_mean_accum_t     signal_mean_accum;
    signal_mean_accum_t     signal_var_accum;
    BlockAccums( prob_t count_total = 0.0 )
    {}
};

BlocksScoring::BlocksScoring(
    const ChessboardBiclusteringsPDFEval& walkStats,
    const IndexedPartitionsCollection<ObjectsCluster>&  objPtnColl,
    const IndexedPartitionsCollection<ProbesCluster>&   probesPtnColl,
    const std::vector<object_clundex_t>&        objectsClusters,
    const std::vector<probe_clundex_t>&         probesClusters
) : _objectsClusters( objectsClusters )
  , _probesClusters( probesClusters )
  , _objClu2Clu( objPtnColl.clusterIntersections( objectsClusters ) )
  , _probeClu2Clu( probesPtnColl.clusterIntersections( probesClusters ) )
{
    typedef boost::unordered_map<block_id, BlockAccums> block_accums_map;
    block_accums_map blockAccums;
    for ( std::size_t i = 0; i < probesClusters.size(); i++ ) {
        for ( std::size_t ii = 0; ii < _probeClu2Clu[i].size(); ii++ ) {
            const ClusterIntersectionStats& probeIsect = _probeClu2Clu[i][ii];
            ChessboardBiclusteringsPDFEval::block_stats_map::const_iterator objMapIt =
                walkStats.blocksStatsMap().find( probeIsect.clusterIx );
            if ( objMapIt == walkStats.blocksStatsMap().end() ) continue;
            const ChessboardBiclusteringsPDFEval::object_block_stats_map& objStats = *objMapIt->second;
            for ( std::size_t j = 0; j < objectsClusters.size(); j++ ) {
                block_id blkId = std::make_pair( objectsClusters[j], probesClusters[i] );
                block_accums_map::iterator blkIt = blockAccums.find( blkId );
                if ( blkIt == blockAccums.end() ) {
                    blkIt = blockAccums.insert( blkIt, std::make_pair( blkId, BlockAccums() ) );
                }
                BlockAccums& accums = blkIt->second;
                for ( std::size_t jj = 0; jj < _objClu2Clu[j].size(); jj++ ) {
                    const ClusterIntersectionStats& objIsect = _objClu2Clu[j][jj];
                    ChessboardBiclusteringsPDFEval::object_block_stats_map::const_iterator oit = objStats.find( objIsect.clusterIx );
                    if ( oit != objStats.end() ) {
                        prob_t blockCoverage = objIsect.coverage * probeIsect.coverage;
                        prob_t weight = oit->second.count_on * blockCoverage;
                        accums.count_total += oit->second.count_total * blockCoverage;
                        if ( weight > 0 ) {
                            accums.signal_mean_accum( oit->second.signal_mean,
                                                      boost::accumulators::weight = weight );
                            accums.signal_var_accum( oit->second.signal_var,
                                                     boost::accumulators::weight = weight );
                        }
                    }
                }
            }
        }
    }
    // extract mean/var from accums
    for ( block_accums_map::const_iterator blkIt = blockAccums.begin();
          blkIt != blockAccums.end(); ++blkIt
    ){
        _blockStats.insert( std::make_pair( blkIt->first,
                            SplitBlockStats(
                                blkIt->second.count_total,
                                boost::accumulators::extract::sum_of_weights( blkIt->second.signal_mean_accum ),
                                boost::accumulators::extract::weighted_mean( blkIt->second.signal_mean_accum ),
                                boost::accumulators::extract::weighted_mean( blkIt->second.signal_var_accum ) ) ) );
    }
}

} }
