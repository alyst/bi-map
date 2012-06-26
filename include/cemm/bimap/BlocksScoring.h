#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#include <boost/unordered_map.hpp>

#include "ChessboardBiclusteringsPDFEval.h"
#include "IndexedPartitionsCollection.h"

namespace cemm { namespace bimap {

/**
 *  Scores the blocks resulting from arbitrary objects and probes clusters.
 *  The score is the weighted sum of on-block frequencies for all
 *  blocks that overlap with the arbitrary block,
 *  where the weight is the ratio of the cells of specified block
 *  covered by the existing block.
 */
class BlocksScoring
{
public:
    typedef std::pair<object_clundex_t, probe_clundex_t> block_id;
    struct SplitBlockStats {
        prob_t count_total; /// total number of object cluster x probe cluster blocks
        prob_t count_on;    /// number of blocks in on-state
        SplitBlockStats( prob_t count_total = 0.0, prob_t count_on = 0.0 )
        : count_total( count_total ), count_on( count_on )
        {}
    };
    typedef boost::unordered_map<block_id, SplitBlockStats> block_stats_map;

private:
    const std::vector<object_clundex_t>&    _objectsClusters;
    const std::vector<probe_clundex_t>&     _probesClusters;

    std::vector< std::vector<ClusterIntersectionStats> >   _objClu2Clu;
    std::vector< std::vector<ClusterIntersectionStats> >   _probeClu2Clu;

    block_stats_map _blockStats;

public:
    BlocksScoring( const ChessboardBiclusteringsPDFEval& walkStats,
                   const IndexedPartitionsCollection<ObjectsCluster>& objPtnColl,
                   const IndexedPartitionsCollection<ProbesCluster>& probesPtnColl,
                   const std::vector<object_clundex_t>&     objectsClusters,
                   const std::vector<probe_clundex_t>&      probesClusters
    );

    const std::vector<object_clundex_t>& objectsClusters() const {
        return ( _objectsClusters );
    }

    const std::vector<probe_clundex_t>& probesClusters() const {
        return ( _probesClusters );
    }

    const block_stats_map& blockStats() const {
        return ( _blockStats );
    }
};

} }