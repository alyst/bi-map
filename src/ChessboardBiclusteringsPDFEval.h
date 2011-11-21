#pragma once

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/ptr_container/serialize_ptr_vector.hpp>
#include <boost/unordered_map.hpp>

#include "BIMAPWalk.h"
#include "PartitionStatistics.h"

/**
 *  PDF of partition based on frequency of its independent components:
 *  partition PDF is multiple of marginal distributions of its components.
 * 
 *  @see ChessboardBiclusteringsPDFEval
 *  @see IndexedPartitionsCollection
 */
class PartitionIndependentComponentsPDF
{
public:
    typedef size_t serial_type;
    typedef std::vector<serial_type> serial_vector_type;

private:
    friend class boost::serialization::access;

    typedef boost::unordered_map<serial_type, serial_vector_type> component_map;
    typedef boost::unordered_map<serial_type, prob_t> component_freq_map;
    typedef boost::unordered_map<serial_type, size_t> component_counts_map;

    serial_vector_type      _empty;
    component_map           _ptnComponentsMap;
    boost::ptr_vector<component_freq_map> _ptnComponentsFreq;

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP( _ptnComponentsMap );
        ar & BOOST_SERIALIZATION_NVP( _ptnComponentsFreq );
    }

public:
    PartitionIndependentComponentsPDF() {};

    PartitionIndependentComponentsPDF( size_t normalizer,
                        const component_map& componentsMap, 
                        const boost::ptr_vector<component_counts_map>& countsMap );

    const serial_vector_type& subpartitions( serial_type ptnSerial ) const {
        component_map::const_iterator ptnIt = _ptnComponentsMap.find( ptnSerial );
        return ( ptnIt != _ptnComponentsMap.end() ? ptnIt->second : _empty );
    }
    log_prob_t lnPdf( serial_type ptnSerial ) const;
};

/**
 *  Adjustment of chessboard biclusterings frequency distribution
 *  by adjusting PDFs of objects, probes partitions and
 *  "on"/"off" cell masks.
 * 
 *  @see PartitionPDFEval
 */
class ChessboardBiclusteringsPDFEval
{
public:
    typedef size_t serial_type;
    typedef serial_type object_serial_type;
    typedef serial_type probe_serial_type;
    typedef std::pair<object_serial_type, probe_serial_type> block_id;
    typedef std::pair<size_t, size_t> block_stats;
    typedef boost::unordered_map<block_id, block_stats> block_stats_map;
    typedef boost::unordered_map<object_serial_type, PartStats> obj_clu_stats_map;
    typedef boost::unordered_map<probe_serial_type, PartStats> probe_clu_stats_map;

private:
    friend class boost::serialization::access;

    typedef size_t submask_id_type; /** id (hash) of bitmask of object x probes independent component */

    typedef ChessboardBiclusteringScaffold::cells_mask_t cells_mask_type;
    typedef boost::unordered_map<submask_id_type, log_prob_t> cells_mask_freq_map;
    typedef PartitionIndependentComponentsPDF::serial_vector_type serial_vector_type;

    PartitionIndependentComponentsPDF   _objPtnIcPdf;
    PartitionIndependentComponentsPDF   _probePtnIcPdf;
    PartitionPartsPDF                   _objPtnPartsPdf;
    PartitionPartsPDF                   _probePtnPartsPdf;
    std::vector<object_set_t>           _objComponents;     /** objects independent components */
    std::vector<probe_bitset_t>         _probeComponents;   /** probes independent components */
    std::vector<cells_mask_freq_map>    _cellsSubmaskFreq;  /** map of cell "on"/"off" masks frequencies
                                                                per component pair */
    block_stats_map                     _blockStats;        /** biclustering blocks statistics */
    obj_clu_stats_map                   _objCluStats;       /** average frequency of co-occurrence of pairs of object in a cluster */
    probe_clu_stats_map                 _probeCluStats;     /** average frequency of co-occurrence of pairs of probes in a cluster */

    void encodeMask( const cells_mask_type& mask,
                     std::vector<submask_id_type>& submaskIds,
                     cells_mask_type& submaskBuf,
                     std::string& submaskStrBuf ) const;

    void evalCellsMaskFreqMap( const BIMAPWalk& walk );
    void evalBlocksFreqMap( const BIMAPWalk& walk );

    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        boost::serialization::collection_size_type stepsCount;
        ar & BOOST_SERIALIZATION_NVP( _cellsSubmaskFreq );
        ar & BOOST_SERIALIZATION_NVP( _objPtnIcPdf );
        ar & BOOST_SERIALIZATION_NVP( _probePtnIcPdf );
        ar & BOOST_SERIALIZATION_NVP( _objPtnPartsPdf );
        ar & BOOST_SERIALIZATION_NVP( _probePtnPartsPdf );
        ar & BOOST_SERIALIZATION_NVP( _objComponents );
        ar & BOOST_SERIALIZATION_NVP( _probeComponents );
    }

    ChessboardBiclusteringsPDFEval() {};

public:
    ChessboardBiclusteringsPDFEval( BIMAPWalk& walk, 
                               prob_t objects_threshold,
                               prob_t probes_threshold );

    const std::vector<object_set_t>& objectComponents() const {
        return ( _objComponents );
    }

    const std::vector<probe_bitset_t>& probeComponents() const {
        return ( _probeComponents );
    }

    const PartitionIndependentComponentsPDF& objectsIcPdf() const {
        return ( _objPtnIcPdf );
    }
    const PartitionIndependentComponentsPDF& probesIcPdf() const {
        return ( _probePtnIcPdf );
    }

    const PartitionPartsPDF& objectsPartsPdf() const {
        return ( _objPtnPartsPdf );
    }
    const PartitionPartsPDF& probesPartsPdf() const {
        return ( _probePtnPartsPdf );
    }

    log_prob_t lnObjectsPartitionIcPdf( const ChessboardBiclusteringScaffold& clustering ) const
    {
        return ( _objPtnIcPdf.lnPdf( clustering.pObjectsPartition->serial() ) );
    }
    log_prob_t lnProbesPartitionIcPdf( const ChessboardBiclusteringScaffold& clustering ) const
    {
        return ( _probePtnIcPdf.lnPdf( clustering.pProbesPartition->serial() ) );
    }
    log_prob_t lnObjectsPartitionPartsPdf( const ChessboardBiclusteringScaffold& clustering ) const
    {
        return ( _objPtnPartsPdf.lnPdf( clustering.pObjectsPartition->serial() ) );
    }
    log_prob_t lnProbesPartitionPartsPdf( const ChessboardBiclusteringScaffold& clustering ) const
    {
        return ( _probePtnPartsPdf.lnPdf( clustering.pProbesPartition->serial() ) );
    }
    log_prob_t lnCellsMaskPdf( const ChessboardBiclusteringScaffold& clustering ) const;

    log_prob_t lnIcPdf( const ChessboardBiclusteringScaffold& clustering ) const
    {
        return (   lnObjectsPartitionIcPdf( clustering )
                 + lnProbesPartitionIcPdf( clustering )
                 + lnCellsMaskPdf( clustering ) );
    }
    log_prob_t lnPartsPdf( const ChessboardBiclusteringScaffold& clustering ) const
    {
        return (   lnObjectsPartitionPartsPdf( clustering )
                 + lnProbesPartitionPartsPdf( clustering )
                 + lnCellsMaskPdf( clustering ) );
    }
    const obj_clu_stats_map& objectsClustersStatsMap() const {
        return ( _objCluStats );
    }
    const probe_clu_stats_map& probesClustersStatsMap() const {
        return ( _probeCluStats );
    }
    block_stats blockStats( const block_id& blockId ) const {
        return ( _blockStats.find( blockId )->second );
    }
    const block_stats_map& blocksStatsMap() const {
        return ( _blockStats );
    }
};
