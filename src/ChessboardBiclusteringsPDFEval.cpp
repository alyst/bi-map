#include "cemm/bimap/ChessboardBiclusteringsPDFEval.h"

#include <cemm/containers/dynamic_bitset_foreach.h>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/variance.hpp>

#include "cemm/bimap/IndexedPartitionsCollection.h"
#include "PartitionStatisticsImpl.h"

namespace cemm { namespace bimap {

PartitionIndependentComponentsPDF::PartitionIndependentComponentsPDF(
    size_t normalizer,
    const component_map& componentsMap,
    const boost::ptr_vector<component_counts_map>& countsMap
) : _ptnComponentsMap( componentsMap )
  , _ptnComponentsFreq( countsMap.size() )
{
    _ptnComponentsFreq.resize( countsMap.size() );
    for ( size_t compIx = 0; compIx < countsMap.size(); compIx++ ) {
        for ( component_counts_map::const_iterator it = countsMap[ compIx ].begin();
              it != countsMap[ compIx ].end(); ++it )
        {
            _ptnComponentsFreq[ compIx ][ it->first ] = (prob_t)it->second / normalizer;
        }
    }
}

log_prob_t PartitionIndependentComponentsPDF::lnPdf(
    serial_type ptnSerial
) const {
    component_map::const_iterator compsIt = _ptnComponentsMap.find( ptnSerial );
    if ( compsIt == _ptnComponentsMap.end() ) {
        //THROW_EXCEPTION( std::invalid_argument, "No components for partition #" << ptnSerial );
        return ( std::numeric_limits<log_prob_t>::signaling_NaN() );
    }
    serial_vector_type  comps = compsIt->second;
    double  res = 0;
    for ( size_t compIx = 0; compIx < comps.size(); compIx++ ) {
        component_freq_map::const_iterator freqIt = _ptnComponentsFreq[ compIx ].find( comps[ compIx ] );
        if ( freqIt == _ptnComponentsFreq[ compIx ].end() ) {
            THROW_RUNTIME_ERROR( "No frequency data for partition #" << ptnSerial << " of component " << compIx );
        }
        res += log( freqIt->second );
    }
    return ( res );
}

ChessboardBiclusteringsPDFEval::ChessboardBiclusteringsPDFEval(
    const BIMAPWalk&  walk,    /** i random walk, the hidden effect is
                                     addition new objects/probes clusters
                                     to the their indices -- due to splitting of existing
                                     clusters by the independent components
                                */
    const IndexedPartitionsCollection<ObjectsCluster>& objsPtnColl,
    const IndexedPartitionsCollection<ProbesCluster>& probesPtnColl
) : _objComponents( objsPtnColl.components() )
  , _objPtnIcPdf( walk.stepsCount(),
                  objsPtnColl.partitionCompositionMap(),
                  objsPtnColl.subpartitionCountsMaps() )
  , _objPtnPartsPdf( objsPtnColl )
  , _objCluStats( objsPtnColl.partStats() )

  , _probeComponents( probesPtnColl.components() )
  , _probePtnIcPdf( walk.stepsCount(),
                    probesPtnColl.partitionCompositionMap(),
                    probesPtnColl.subpartitionCountsMaps() )
  , _probePtnPartsPdf( PartitionPartsPDF( probesPtnColl ) )
  , _probeCluStats( probesPtnColl.partStats() )
{
    evalCellsMaskFreqMap( walk );
    evalBlocksFreqMap( walk );
}

/**
 *  Adjusts the probability density function of the models
 *  by identifying independent components in the partitions.
 */
ChessboardBiclusteringsPDFEval* ChessboardBiclusteringsPDFEval::Create(
    BIMAPWalk&     walk,       /** i random walk, the hidden effect is
                                      addition new objects/probes clusters
                                      to the their indices -- due to splitting of existing
                                      clusters by the independent components
                                  */
    gsl_rng*        rng,
    prob_t          objects_threshold, /** i frequency threshold for objects to be
                                             considered co-clustered */
    prob_t          probes_threshold,  /** i frequency threshold for probes to be
                                             considered co-clustered */
    prob_t          objects_clot_threshold /** i frequency threshold for objects forming 'clots' */
){
    IndexedPartitionsCollection<ObjectsCluster> objPtnColl( walk );
    objPtnColl.init( rng, walk.stepsCount() * objects_threshold, walk.stepsCount() * objects_clot_threshold );

    IndexedPartitionsCollection<ProbesCluster> probePtnColl( walk );
    probePtnColl.init( NULL, walk.stepsCount() * probes_threshold, std::numeric_limits<prob_t>::quiet_NaN() );
    return ( new ChessboardBiclusteringsPDFEval( walk, objPtnColl, probePtnColl ) );
}

/**
 *  Encodes given mask into matrix of submask Ids.
 */
void ChessboardBiclusteringsPDFEval::encodeMask(
    const cells_mask_type&          mask,           /** i total mask of cells probes */
    std::vector<submask_id_type>&   submaskIds,     /** o object x probes matrix of submask IDs (per independent component pair */
    cells_mask_type&                submaskBuf,     /** o buffer for intermediate submask storage */
    std::string&                    submaskStrBuf   /** o buffer for ID generation */
) const {
    BOOST_ASSERT( submaskIds.size() == _objComponents.size() * _probeComponents.size() );
    // extract submask of cells for each object x probe independent component pair
    // and register it in array
    for ( size_t objCompIx = 0; objCompIx < _objComponents.size(); objCompIx++ ) {
        const object_set_t& objs = _objComponents[ objCompIx ];
        for ( size_t probeCompIx = 0; probeCompIx < _probeComponents.size(); probeCompIx++ ) {
            // extract submask
            size_t ccompIx = objCompIx * _probeComponents.size() + probeCompIx;
            const probe_bitset_t& probes = _probeComponents[ probeCompIx ];
            submaskBuf.resize( objs.size() * probes.count() );
            submaskBuf.reset();
            size_t offset = 0;
            for ( object_set_t::const_iterator oit = objs.begin(); oit != objs.end(); ++oit ) {
                size_t maskOffset = (*oit) * probes.size();
                foreach_bit( probe_index_t, probeIx, probes ) {
                    if ( mask.test( maskOffset + probeIx ) ) submaskBuf.set( offset );
                    offset++;
                }
            }
            // update counts of this submask
            boost::to_string( submaskBuf, submaskStrBuf );
            submaskIds[ ccompIx ] = boost::hash_value( submaskStrBuf );
        }
    }
}

/**
 *  Calculate frequency of objects x probes cells "on"/"off" pattern.
 */
void ChessboardBiclusteringsPDFEval::evalCellsMaskFreqMap(
    const BIMAPWalk& walk
){
    typedef boost::unordered_map<submask_id_type, size_t> submask_id_counts_map;

    // collection of maps from submask id to its counts,
    // per independent component pair
    std::vector<submask_id_counts_map> submasksCounts( _objComponents.size() * _probeComponents.size() );

    cells_mask_type mask; // buffer for complete mask
    cells_mask_type submask; // buffer for component's submask
    std::string maskString;
    std::vector<submask_id_type> submaskIds( submasksCounts.size() );

    for ( BIMAPWalk::const_step_iterator stepIt = walk.stepsBegin(); stepIt != walk.stepsEnd(); ++stepIt )
    {
        const ChessboardBiclusteringScaffold& clustering = stepIt->clustering.scaffold();
        clustering.getCellsMask( mask );
        encodeMask( mask, submaskIds, submask, maskString );
        // extract submask of cells for each object x probe independent component pair
        // and register it in array
        for ( size_t submaskIx = 0; submaskIx < submaskIds.size(); submaskIx++ ) {
            submask_id_counts_map& submaskCounts = submasksCounts[ submaskIx ];
            submask_id_type submaskId = submaskIds[ submaskIx ];
            submask_id_counts_map::iterator cit = submaskCounts.find( submaskId );
            if ( cit != submaskCounts.end() ) {
                cit->second++;
            } else {
                submaskCounts[ submaskId ] = 1;
            }
        }
    }
    // fill submask frequency maps
    log_prob_t ln_normalizer = log( walk.stepsCount() );
    _cellsSubmaskFreq.resize( submasksCounts.size() );
    for ( size_t submaskIx = 0; submaskIx < submaskIds.size(); submaskIx++ ) {
        submask_id_counts_map& countsMap = submasksCounts[ submaskIx ];

        cells_mask_freq_map& freqMap = _cellsSubmaskFreq[ submaskIx ];
        freqMap.clear();

        for ( submask_id_counts_map::const_iterator cit = countsMap.begin(); cit != countsMap.end(); ++cit )
        {
            freqMap[ cit->first ] = log( (prob_t)cit->second ) - ln_normalizer;
        }
    }
}

/**
 *  Calculate frequency of blocks in the samples, and frequency of their on/off probes.
 */
void ChessboardBiclusteringsPDFEval::evalBlocksFreqMap(
    const BIMAPWalk& walk
){
    typedef boost::accumulators::accumulator_set<
            signal_t, boost::accumulators::stats<
                boost::accumulators::tag::mean,
                boost::accumulators::tag::variance> > signal_accum_t;

    typedef ChessboardBiclusteringIndexed::block_key_type blk_id_t;
    typedef boost::unordered_map<blk_id_t, signal_accum_t> signal_accum_map;

    signal_accum_map blockAccums; // internal blocks signal statistics

    // collection of maps from submask id to its counts,
    // per independent component pair
    _blockStats.clear();
    for ( BIMAPWalk::const_step_iterator stepIt = walk.stepsBegin(); stepIt != walk.stepsEnd(); ++stepIt )
    {
        const ChessboardBiclusteringScaffold& clustering = stepIt->clustering.scaffold();
        for ( ChessboardBiclusteringScaffold::const_probe_cluster_iterator prcluIt = clustering.probeClusterBegin();
              prcluIt != clustering.probeClusterEnd(); ++prcluIt
        ){
            probe_clundex_t prcluIx = (*prcluIt)->serial();
            block_stats_map::iterator objBlockMapIt = _blockStats.find( prcluIx );
            if ( objBlockMapIt == _blockStats.end() ) {
                objBlockMapIt = _blockStats.insert( objBlockMapIt, prcluIx, new object_block_stats_map() );
            }
            object_block_stats_map& objBlockMap = *objBlockMapIt->second;

            for ( ChessboardBiclusteringScaffold::const_object_cluster_iterator ocluIt = clustering.objectClusterBegin();
                ocluIt != clustering.objectClusterEnd(); ++ocluIt
            ){
                object_clundex_t objcluIx = (*ocluIt)->serial();
                bool isEnabled = clustering.isBlockEnabled( objcluIx, prcluIx );
                object_block_stats_map::iterator blkIt = objBlockMap.find( objcluIx );
                if ( blkIt == objBlockMap.end() ) {
                    objBlockMap.insert( blkIt, std::make_pair( objcluIx,
                                        block_stats( 1, isEnabled ? 1 : 0 ) ) );
                } else {
                    blkIt->second.count_total++;
                    if ( isEnabled ) blkIt->second.count_on++;
                }
            }
        }
        // store signals to accums
        ChessboardBiclusteringIndexed::block_data_map_type::const_iterator blkItEnd = stepIt->clustering.blocksData().end();
        for ( ChessboardBiclusteringIndexed::block_data_map_type::const_iterator blkIt = stepIt->clustering.blocksData().begin();
              blkIt != blkItEnd; ++blkIt
        ){
            signal_accum_map::iterator blkAccIt = blockAccums.find( blkIt->first );
            if ( blkAccIt == blockAccums.end() ) {
                blkAccIt = blockAccums.insert( blkAccIt, std::make_pair( blkIt->first, signal_accum_t() ) );
            }
            blkAccIt->second( blkIt->second );
        }
    }
    // copy the signal statistics from the accumulators
    for ( block_stats_map::iterator it = _blockStats.begin(); it != _blockStats.end(); ++it )
    {
        object_block_stats_map& objMap = *it->second;
        for ( object_block_stats_map::iterator jt = objMap.begin(); jt != objMap.end(); ++jt ) {
            blk_id_t blkId = blk_id_t( jt->first, it->first );
            signal_accum_map::iterator accIt = blockAccums.find( blkId );
            if ( accIt != blockAccums.end() ) {
                jt->second.signal_mean = boost::accumulators::extract::mean( accIt->second );
                jt->second.signal_var = boost::accumulators::extract::variance( accIt->second );
            } else {
                LOG_WARN( "Pair (" << jt->first << ", " << it->first << ")"
                          << " was not found in accumulators map" );
            }
        }
    }
}

/**
 *  Calculates the adjusted frequency of enabled cells mask.
 */
log_prob_t ChessboardBiclusteringsPDFEval::lnCellsMaskPdf(
    const ChessboardBiclusteringScaffold& clustering
) const {
    cells_mask_type mask;
    clustering.getCellsMask( mask );
    cells_mask_type submask; // buffer for component's submask
    std::string maskString;
    std::vector<submask_id_type> submaskIds( _objComponents.size() * _probeComponents.size() );
    encodeMask( mask, submaskIds, submask, maskString );

    log_prob_t res = 0;
    for ( size_t submaskIx = 0; submaskIx < submaskIds.size(); submaskIx++ ) {
        const cells_mask_freq_map& freqMap = _cellsSubmaskFreq[ submaskIx ];
        cells_mask_freq_map::const_iterator cfmit = freqMap.find( submaskIds[ submaskIx ] );
        if ( cfmit != freqMap.end() ) {
            res += cfmit->second;
        } else {
            // no frequency for given submask
            return ( std::numeric_limits<log_prob_t>::quiet_NaN() );
        }
    }
    return ( res );
}

} }
