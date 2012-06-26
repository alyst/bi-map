#include "cemm/bimap/BlocksScoring.h"

#include "cemm/bimap/IndexedPartitionsCollection.h"

namespace cemm { namespace bimap {

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
    for ( std::size_t i = 0; i < probesClusters.size(); i++ ) {
        for ( std::size_t ii = 0; ii < _probeClu2Clu[i].size(); ii++ ) {
            const ClusterIntersectionStats& probeIsect = _probeClu2Clu[i][ii];
            ChessboardBiclusteringsPDFEval::block_stats_map::const_iterator objMapIt =
                walkStats.blocksStatsMap().find( probeIsect.clusterIx );
            if ( objMapIt == walkStats.blocksStatsMap().end() ) continue;
            const ChessboardBiclusteringsPDFEval::object_block_stats_map& objStats = *objMapIt->second;
            for ( std::size_t j = 0; j < objectsClusters.size(); j++ ) {
                SplitBlockStats stats;
                for ( std::size_t jj = 0; jj < _objClu2Clu[j].size(); jj++ ) {
                    const ClusterIntersectionStats& objIsect = _objClu2Clu[j][jj];
                    ChessboardBiclusteringsPDFEval::object_block_stats_map::const_iterator oit = objStats.find( objIsect.clusterIx );
                    if ( oit != objStats.end() ) {
                        prob_t blockCoverage = objIsect.coverage * probeIsect.coverage;
                        stats.count_total += oit->second.count_total * blockCoverage;
                        stats.count_on += oit->second.count_on * blockCoverage;
                    }
                }
                block_id blkId = std::make_pair( objectsClusters[j], probesClusters[i] );
                block_stats_map::iterator blkIt = _blockStats.find( blkId );
                if ( blkIt == _blockStats.end() ) {
                    _blockStats.insert( blkIt, std::make_pair( blkId, stats ) );
                } else {
                    blkIt->second.count_on += stats.count_on;
                    blkIt->second.count_total += stats.count_total;
                }
            }
        }
    }
}

} }