#include "ChessboardBiclusteringsIndexing.h"
#include "dynamic_bitset_utils.h"

object_clundex_t ChessboardBiclusteringScaffold::objectsClusterIndex(
    objects_cluster_serial_type objectsCluSerial
) const {
    object_clundex_t cluIx = 0;
    for ( const_object_cluster_iterator cluIt = objectClusterBegin(); cluIt != objectClusterEnd(); ++cluIt ){
        if ( (*cluIt)->serial() == objectsCluSerial ) return ( cluIx );
        cluIx++;
    }
    return ( CLUSTER_NA );
}

probe_clundex_t ChessboardBiclusteringScaffold::probesClusterIndex(
    probes_cluster_serial_type probesCluSerial
) const {
    probe_clundex_t cluIx = 0;
    for ( const_probe_cluster_iterator cluIt = probeClusterBegin(); cluIt != probeClusterEnd(); ++cluIt ){
        if ( (*cluIt)->serial() == probesCluSerial ) return ( cluIx );
        cluIx++;
    }
    return ( CLUSTER_NA );
}

bool ChessboardBiclusteringScaffold::isCrossClusterEnabled(
    objects_cluster_serial_type     objCluSerial,
    probes_cluster_serial_type      probeCluSerial
) const {
    object_clundex_t objCluIx = objectsClusterIndex( objCluSerial );
    probe_clundex_t probeCluIx = probesClusterIndex( probeCluSerial );
    return ( objCluIx != CLUSTER_NA && probeCluIx != CLUSTER_NA
             && crossClusterMask.test( objCluIx * pProbesPartition->value().size() + probeCluIx ) );
}

/**
 *  Gets mask of enabled/disabled probes of cells of object x probes
 *  array, induced by crossClustersMask().
 */
void ChessboardBiclusteringScaffold::getCellsMask(
    cells_mask_t&   mask /** o  objects x probes mask of enabled/disabled probes */
) const {
    // get the number of objects and probes
    size_t  nProbeClus = pProbesPartition->value().size();
    size_t  nProbes = (*probeClusterBegin())->value().size();
    size_t  nObjects = 0;
    for ( const_object_cluster_iterator oit = objectClusterBegin(); 
          oit != objectClusterEnd(); ++oit
    ){
        nObjects += (*oit)->value().size();
    }

    // reset mask
    mask.resize( nProbes * nObjects );
    mask.reset();

    foreach_bit( size_t, ccIx, crossClusterMask ) {
        // enable cells of each enabled cross-cluster
        object_clundex_t objCluIx = ccIx / nProbeClus;
        probe_clundex_t  probeCluIx = ccIx % nProbeClus;

        const_object_cluster_iterator oCluIt = objectClusterBegin();
        std::advance( oCluIt, objCluIx );
        const_probe_cluster_iterator sCluIt = probeClusterBegin();
        std::advance( sCluIt, probeCluIx );

        const object_set_t&     objs = (*oCluIt)->value();
        const probe_bitset_t&   probes = (*sCluIt)->value();
        for ( object_set_t::const_iterator oit = objs.begin(); oit != objs.end(); ++oit ) {
            foreach_bit( probe_index_t, probeIx, probes ) {
                mask.set( *oit * nProbes + probeIx );
            }
        }
    }
}

bool ChessboardBiclusteringScaffold::check() const
{
    if ( pProbesPartition->value().empty() ) {
        THROW_RUNTIME_ERROR( "Probes partition #" << pProbesPartition->serial() << " is empty" );
    }
    if ( pObjectsPartition->value().empty() ) {
        std::ostringstream out;
        THROW_RUNTIME_ERROR( "Objects partition #" << pObjectsPartition->serial() << " is empty" );
    }
    size_t expCCMaskSize = pProbesPartition->value().size() * pObjectsPartition->value().size();
    if ( crossClusterMask.size() != expCCMaskSize ) {
        THROW_RUNTIME_ERROR( "Incorrect cross clusters mask: size " << crossClusterMask.size() << ", should be " << expCCMaskSize );
    }
    return ( true );
}

bool ChessboardBiclusteringIndexed::check() const
{
    scaffold().check();

    // check the signals
    for ( ChessboardBiclusteringScaffold::const_object_cluster_iterator oit = scaffold().objectClusterBegin();
          oit != scaffold().objectClusterEnd(); ++oit
    ){
        objects_cluster_serial_type objCluSerial = (*oit)->serial();
        object_clundex_t objCluIx = scaffold().objectsClusterIndex( objCluSerial );
        for ( ChessboardBiclusteringScaffold::const_probe_cluster_iterator sit = scaffold().probeClusterBegin();
            sit != scaffold().probeClusterEnd(); ++sit
        ){
            probes_cluster_serial_type probeCluSerial = (*sit)->serial();
            probe_clundex_t probeCluIx = scaffold().probesClusterIndex( probeCluSerial );
            cross_cluster_data_map_type::const_iterator ccDataIt = crossClustersData().find( cross_cluster_key_type( objCluSerial, probeCluSerial ) );
            if ( scaffold().isCrossClusterEnabled( objCluSerial, probeCluSerial ) ) {
                if ( ccDataIt != crossClustersData().end() ) {
                    if ( is_unset( ccDataIt->second ) ) {
                        THROW_RUNTIME_ERROR( "Unset (NaN) signal for object cluster #" << objCluIx
                                                << ", probe cluster #" << probeCluIx );
                    }
                } else {
                        THROW_RUNTIME_ERROR( "No signal for object cluster #" << objCluIx
                                                << ", probe cluster #" << probeCluIx );
                }
            } else if ( ccDataIt != crossClustersData().end() ) {
                THROW_RUNTIME_ERROR( "Signal exists for disabled cross cluster (" << objCluIx
                                        << ", " << probeCluIx << ")" );
            }
        }
    }
    return ( true );
}

ChessboardBiclusteringsIndexing::ChessboardBiclusteringsIndexing(
    size_t      min_size,
    size_t      max_size
) : EntityIndexing< ChessboardBiclusteringScaffold, default_serial_t >( min_size, max_size)
 , _objectPartitionIndexing( min_size, max_size )
 , _probePartitionIndexing( min_size, max_size )
{

}

ChessboardBiclusteringIndexed ChessboardBiclusteringsIndexing::index(
    const ChessboardBiclustering& clustering
){
    clustering.check();

    ChessboardBiclusteringScaffold    scaffold;
    ChessboardBiclusteringIndexed::cross_cluster_data_map_type     crossClusterData;
    scaffold.pObjectsPartition = _objectPartitionIndexing.index( clustering.objectsClusters().begin(), clustering.objectsClusters().end() );
    scaffold.pProbesPartition = _probePartitionIndexing.index( clustering.probesClusters().begin(), clustering.probesClusters().end() );
    // create a new cross clusters map, where object clusters and probe clusters would be reindexed with respect to their order in indexed partition
    scaffold.crossClusterMask = ChessboardBiclusteringScaffold::cross_cluster_mask_t( clustering.objectsClusters().size() * clustering.probesClusters().size() );

    for ( ChessboardBiclustering::const_cross_cluster_iterator ccIt = clustering.begin(); ccIt != clustering.end(); ++ccIt ) {
        ChessboardBiclusteringIndexed::objects_cluster_serial_type objCluSerial = _objectPartitionIndexing.elementIndexing()
                    .index( ccIt->objectsCluster().items() )->serial();
        // get the index of the cluster in the indexed partition
        object_clundex_t newObjCluIx = scaffold.objectsClusterIndex( objCluSerial );
        ChessboardBiclusteringIndexed::probes_cluster_serial_type probeCluSerial = _probePartitionIndexing.elementIndexing()
                    .index( ccIt->probesCluster().items() )->serial();
        // get the index of the cluster in the indexed partition
        probe_clundex_t newProbeCluIx = scaffold.probesClusterIndex( probeCluSerial );
        scaffold.crossClusterMask.set( newObjCluIx * clustering.probesClusters().size() + newProbeCluIx );
        crossClusterData[ ChessboardBiclusteringIndexed::cross_cluster_key_type( objCluSerial, probeCluSerial ) ] = ccIt->signal();
    }
    BOOST_ASSERT( scaffold.crossClusterMask.count() == clustering.crossClustersMask().count() );
    base_type::entity_pointer_type pScaffold = base_type::index( scaffold );

    return ( ChessboardBiclusteringIndexed( pScaffold, crossClusterData, clustering.objectMultiples(), clustering.clusteringData() ) );
}
