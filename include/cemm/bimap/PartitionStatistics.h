#pragma once

#include <boost/unordered/unordered_map.hpp>

namespace cemm { namespace bimap {

/**
 *  Partition part statistics.
 */
struct PartStats {
    size_t  size;               /** part's size */
    size_t  nsteps;             /** N of times (stepwise) this part is seen in the partitions of the collection */
    size_t  nstepsIncluded;     /** N of times (stepwise) this part is seen in the partitions of the collection itself or as a subset of other part */
    prob_t  avgPairCoOccurrence; /** average co-occurrence of element pairs of given part
                                    in all partitions of the collection */
    size_t  componentIndex;     /** index of the partition component the part belongs to */

    PartStats( size_t size = 0, size_t nsteps = 0, size_t componentIndex = (size_t)(-1) )
    : size( size ), nsteps( nsteps ), nstepsIncluded( 0 ), avgPairCoOccurrence( 0.0 )
    , componentIndex( componentIndex )
    {};
};

typedef boost::unordered_map<size_t, PartStats>  part_stats_map_type;
typedef boost::unordered_map<size_t, size_t> counts_map_type;

template<typename Part>
class IndexedPartitionsCollection;

/**
 *  PDF of partition is multiple of frequency of its parts (clusters).
 * 
 *  @see ChessboardBiclusteringsPDFEval
 *  @see IndexedPartitionsCollection
 */
class PartitionPartsPDF
{
public:
    typedef size_t partition_serial;

private:
    friend class boost::serialization::access;

    typedef boost::unordered_map<partition_serial, prob_t> ptn_freq_map;

    ptn_freq_map    _ptnPartsFreq;

   template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( _ptnPartsFreq );
    }

public:
    PartitionPartsPDF() {};

    template<typename Part>
    PartitionPartsPDF( const IndexedPartitionsCollection<Part>& ptnColl );

    PartitionPartsPDF& operator=( const PartitionPartsPDF& t )
    {
        _ptnPartsFreq = t._ptnPartsFreq;
        return ( *this );
    }

    log_prob_t lnPdf( partition_serial ptnSerial ) const {
        return ( _ptnPartsFreq.find( ptnSerial )->second );
    }
};

} }
