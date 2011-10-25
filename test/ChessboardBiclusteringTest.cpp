#include <set>

#define TEST_CHECKPOINT( a )

#include "TestCommon.h"
#include "OPADataTest.h"

#include <fstream>
#include <boost/archive/tmpdir.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/filesystem.hpp>

#include <ChessboardBiclustering.h>
#include <ChessboardBiclusteringsIndexing.h>
#include <BIMAPSampler.h>
#include <OPAData.h>
#include <OPADataImportCSV.h>

/*************************************************/

TEST_F2( GslRngTestF, ChessboardBiclusteringOld, DISABLED_edit )
{
    OPAData                     data = generateTestOPAData();
    ChessboardBiclusteringPriors       priors;
    priors.objectClustering.concentration = 0.8;
    priors.probeClustering.concentration = 0.8;
    ObjectsClusterSignal        dummy;
    PitmanYorSample             objectClusters = priors.objectClustering.random( rndGen, data.objectsCount() );
    PitmanYorSample             probeClusters = priors.objectClustering.random( rndGen, data.probesCount() );
    ChessboardBiclustering             clus( ChessboardBiclusteringDerivedPriors(), dummy, GeometricDistribution::BySuccessRate( 0.9 ),
                                      objectClusters, probeClusters );

    EXPECT_EQ( ((const OPAData)data).object( "A" )->sequenceLength(), 2 );
    TEST_CHECKPOINT( "Checking initial object clusters assignment" );
    EXPECT_EQ( clus.objectsCount(), objectClusters.samplesCount() );
    EXPECT_EQ( clus.objectsClusters().size(), objectClusters.clustersCount() );
    for ( size_t i = 0; i < objectClusters.samplesCount(); i++ ) {
        EXPECT_EQ( clus.clusterOfObject( i ), objectClusters.clusterOfSample( i ) );
        EXPECT_TRUE( clus.objectsCluster( objectClusters.clusterOfSample( i ) ).items().find( i )
                    != clus.objectsCluster( objectClusters.clusterOfSample( i ) ).items().end() );
    }

    TEST_CHECKPOINT( "Checking initial probe clusters assignment" );
    EXPECT_EQ( clus.probesCount(), probeClusters.samplesCount() );
    EXPECT_EQ( clus.probesClusters().size(), probeClusters.clustersCount() );
    for ( size_t i = 0; i < probeClusters.samplesCount(); i++ ) {
        EXPECT_EQ( clus.clusterOfProbe( i ), probeClusters.clusterOfSample( i ) );
        EXPECT_TRUE( clus.probesCluster( probeClusters.clusterOfSample( i ) ).items().test( i ) );
    }

    TEST_CHECKPOINT( "Object cluster addition" );
    object_clundex_t oCluIx = clus.addObjectCluster( 1 );
    EXPECT_EQ( clus.clusterOfObject( 1 ), oCluIx );
    EXPECT_TRUE( clus.objectsCluster( oCluIx ).items().find( 1 ) != clus.objectsCluster( oCluIx ).items().end() );

    TEST_CHECKPOINT( "Object-to-cluster assignment" );
    clus.setObjectCluster( 1, 0 );
    EXPECT_EQ( clus.clusterOfObject( 1 ), 0 );
    EXPECT_TRUE( clus.objectsCluster( 0 ).items().find( 1 ) != clus.objectsCluster( 1 ).items().end() );
    EXPECT_TRUE( clus.objectsCluster( oCluIx ).items().find( 1 ) == clus.objectsCluster( oCluIx ).items().end() );
    
    TEST_CHECKPOINT( "Empty objects cluster removal" );
    EXPECT_TRUE( clus.objectsCluster( oCluIx ).items().empty() );
    clus.cleanupClusters();
    EXPECT_GE( oCluIx, clus.objectsClusters().size() );

    TEST_CHECKPOINT( "Probe cluster addition" );
    probe_clundex_t sCluIx = clus.addProbeCluster( 2 );
    EXPECT_EQ( clus.clusterOfProbe( 2 ), sCluIx );
    EXPECT_TRUE( clus.probesCluster( sCluIx ).items().test( 2 ) );

    TEST_CHECKPOINT( "Probe-to-cluster assignment" );
    clus.setProbeCluster( 2, 0 );
    EXPECT_EQ( clus.clusterOfProbe( 2 ), 0 );
    EXPECT_TRUE( !clus.probesCluster( sCluIx ).items().test( 2 ) );
    EXPECT_TRUE( clus.probesCluster( 0 ).items().test( 2 ) );
    
    TEST_CHECKPOINT( "Empty probes cluster removal" );
    EXPECT_TRUE( clus.probesCluster( sCluIx ).items().none() );
    clus.cleanupClusters();
    EXPECT_GE( sCluIx, clus.probesClusters().size() );

#if 0
    TEST_CHECKPOINT( "Cluster indices correctness" );
    EXPECT_EQ( clus[0].index(), 0  );
    EXPECT_EQ( clus[1].index(), 1  );
    EXPECT_EQ( clus[2].index(), 2  );

    TEST_CHECKPOINT( "Object inserting/multiples setting" );
    EXPECT_TRUE( !clus.objectClusters( 0 ).test( 0 ) );
    EXPECT_TRUE( clus[ 0 ].isEmpty() );
    EXPECT_TRUE( !clus[ 0 ].containsObject( 0 ) );
    EXPECT_EQ( clus[ 0 ].objectMultiple( 0 ), 0 );
    EXPECT_TRUE_THROW( clus[ 0 ].setObject( 0 ), std::runtime_error );
    clus[ 0 ].unlockObjectSet();
    clus[ 0 ].setObject( 0 );
    EXPECT_TRUE( clus[ 0 ].isEmpty() );
    EXPECT_TRUE( clus[ 0 ].containsObject( 0 ) );
    EXPECT_EQ( clus[ 0 ].objectMultiple( 0 ), 0 );
    clus[ 0 ].objectMultiples()[ 0 ] = 2;
    EXPECT_TRUE( clus[ 0 ].isEmpty() );
    EXPECT_EQ( clus[ 0 ].objectMultiple( 0 ), 2 );
    EXPECT_TRUE( clus.objectClusters( 0 ).test( 0 ) );
    EXPECT_TRUE( !clus.objectClusters( 0 ).test( 1 ) );
    clus[ 0 ].lockObjectSet();
    
    clus[ 0 ].pushObjectMultiples();
    
    TEST_CHECKPOINT( "Probe inserting/signals setting" );
    EXPECT_TRUE_THROW( clus[ 0 ].setProbe( 0 ), std::runtime_error );
    clus[ 1 ].unlockProbeSet();
    EXPECT_TRUE( !clus.probeClusters( 0 ).test( 0 ) );
    EXPECT_TRUE( !clus.probeClusters( 0 ).test( 1 ) );
    EXPECT_TRUE( clus[ 1 ].isEmpty() );
    EXPECT_TRUE( !clus[ 1 ].containsProbe( 0 ) );
    clus[ 1 ].probeSignals()[ 0 ] = 1.5;
    EXPECT_TRUE( clus[ 1 ].isEmpty() );
    EXPECT_EQ( clus[ 1 ].probeSignal( 0 ), 1.5 );
    EXPECT_TRUE( !clus[ 1 ].containsProbe( 0 ) );
    clus[ 1 ].setProbe( 0 );
    EXPECT_TRUE( clus[ 1 ].containsProbe( 0 ) );
    EXPECT_TRUE( clus.probeClusters( 0 ).test( 1 ) );
    clus[ 1 ].lockProbeSet();
    
    clus[ 1 ].pushProbeSignals();

    TEST_CHECKPOINT( "Simultaneous object+probe inserting/multiples+signals setting" );
    EXPECT_TRUE( clus[ 2 ].isEmpty() );
    clus[2].unlockObjectSet();
    clus[2].unlockProbeSet();
    clus[2].setObject( 1 );
    clus[2].setProbe( 1 );
    EXPECT_TRUE( !clus[ 2 ].isEmpty() );
    EXPECT_TRUE( !clus[ 2 ].containsProbe( 0 ) );
    EXPECT_TRUE( clus[ 2 ].containsProbe( 1 ) );
    EXPECT_TRUE( clus[ 2 ].containsObject( 1 ) );
    EXPECT_TRUE( clus.probeClusters( 1 ).test( 2 ) );
    EXPECT_TRUE( clus.objectClusters( 1 ).test( 2 ) );
    // checking the accessability of multiples/signals
    EXPECT_TRUE_THROW( clus[2].objectMultiples(), std::runtime_error );
    EXPECT_TRUE_THROW( clus[2].objectMultiple( 0 ), std::runtime_error );
    EXPECT_TRUE_THROW( clus[2].probeSignals(), std::runtime_error );
    EXPECT_TRUE_THROW( clus[2].probeSignal( 0 ), std::runtime_error );
    clus[2].lockObjectSet();
    EXPECT_TRUE_NO_THROW( clus[2].probeSignals() );
    EXPECT_TRUE_NO_THROW( clus[2].probeSignal( 0 ) );
    clus[2].lockProbeSet();
    EXPECT_TRUE_NO_THROW( clus[2].objectMultiple( 0 ) );
    EXPECT_TRUE_NO_THROW( clus[2].objectMultiples() );
    
    TEST_CHECKPOINT( "Simultaneous object+probe deletion/multiples+signals setting" );
    /// @warning - it's good to delete the object automatically on setting its count to zero
    clus[2].unlockObjectSet();
    clus[2].setObject( 0 );
    EXPECT_TRUE_NO_THROW( clus[2].objectMultiples()[ 0 ] = 0 );
    EXPECT_TRUE( clus[2].containsObject( 0 ) );
    EXPECT_TRUE( clus.objectClusters( 0 ).test( 2 ) );
    clus[2].setObject( 0, false );
    EXPECT_TRUE( !clus[2].containsObject( 0 ) );
    EXPECT_TRUE( !clus.objectClusters( 0 ).test( 2 ) );

    clus[2].setObject( 1, false );
    EXPECT_TRUE( !clus[2].containsObject( 1 ) );
    EXPECT_TRUE( !clus.objectClusters( 1 ).test( 2 ) );
    EXPECT_TRUE( clus[2].isEmpty() );
    clus[2].lockObjectSet();

    clus[1].unlockProbeSet();
    clus[1].setProbe( 1, false );
    EXPECT_TRUE( !clus[1].containsProbe( 2 ) );
    EXPECT_TRUE( !clus.probeClusters( 2 ).test( 1 ) );

    EXPECT_TRUE( clus[1].containsProbe( 0 ) );
    clus[1].setProbe( 0, false );
    EXPECT_TRUE( !clus[1].containsProbe( 0 ) );
    EXPECT_TRUE( !clus.probeClusters( 0 ).test( 1 ) );
    clus[1].lockProbeSet();
#endif
}

TEST( ChessboardBiclustering, generate_random )
{
    OPAData                     data = generateTestOPAData();
    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    ChessboardBiclusteringPriors       priors;
    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    GibbsSamplerParams          params;
    priors.objectClustering.concentration = 0.8;
    priors.probeClustering.concentration = 0.8;
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    BIMAPSamplerHelper         helper( precomputed, hyperpriors, priors, params );
    ChessboardBiclustering             clus = helper.randomClustering();
    //EXPECT_TRUE( !clusIndexed.objectClusterData().begin()->second.empty() );

    EXPECT_TRUE( clus.checkObjectsPartition() );
    EXPECT_TRUE( clus.checkProbesPartition() );
    EXPECT_TRUE( clus.checkCrossClusters() );
    
    TEST_CHECKPOINT( "Cross cluster enabling/disabling" );
    clus.setCrossCluster( 0, 0, false );
    EXPECT_EQ( clus.isCrossClusterEnabled( 0, 0 ), false );
    clus.setCrossCluster( 0, 0, true );
    EXPECT_EQ( clus.isCrossClusterEnabled( 0, 0 ), true );
}

TEST( ChessboardBiclustering, indexing )
{
    OPAData                     data = generateTestOPAData();
    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    ChessboardBiclusteringPriors       priors;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    GibbsSamplerParams          params;
    priors.objectClustering.concentration = 0.8;
    priors.probeClustering.concentration = 0.8;
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    BIMAPSamplerHelper         helper( precomputed, hyperpriors, priors, params );
    ChessboardBiclustering             clus = helper.randomClustering();
    EXPECT_TRUE( clus.checkObjectsPartition() );
    EXPECT_TRUE( clus.checkProbesPartition() );
    EXPECT_TRUE( clus.checkCrossClusters() );

    ChessboardBiclusteringIndexed clusIndexed = crossClusteringsIndexing.index( clus );
    EXPECT_EQ( clusIndexed.objectsClusters().size(), clus.objectsClusters().size() );
    EXPECT_EQ( clusIndexed.probesClusters().size(), clus.probesClusters().size() );
    EXPECT_EQ( clusIndexed.objectsData().size(), data.objectsCount() );
    EXPECT_NO_THROW( clusIndexed.check() );
    //EXPECT_EQ( clusIndexed.objectClusterData().size(), clus.objectsClusters().size() );
    //EXPECT_TRUE( !clusIndexed.objectClusterData().begin()->second.empty() );
}

TEST( ChessboardBiclustering, csv_in )
{
    boost::filesystem::path test_data_path( BIMAPtest_data_path );
    if ( !boost::filesystem::exists( test_data_path ) || !boost::filesystem::is_directory( test_data_path ) ) {
        THROW_RUNTIME_ERROR( "OPAData test data folder not found: " << test_data_path );
    }

    OPAData                     data = OPADataImportCSV(
        ( test_data_path / "test_proteins.csv" ).string().c_str(),
        ( test_data_path / "test_exp_design.csv" ).string().c_str(),
        ( test_data_path / "test_measurements.csv" ).string().c_str()
    );
}

TEST( ChessboardBiclustering, DISABLED_serialization_out )
{
    OPAData                     data = generateTestOPAData();
    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    ChessboardBiclusteringPriors       priors;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    GibbsSamplerParams          params;
    priors.objectClustering.concentration = 0.8;
    priors.probeClustering.concentration = 0.8;
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    BIMAPSamplerHelper         helper( precomputed, hyperpriors, priors, params );
    ChessboardBiclustering             clusOut = helper.randomClustering();
    EXPECT_TRUE( clusOut.checkObjectsPartition() );
    EXPECT_TRUE( clusOut.checkProbesPartition() );
    EXPECT_TRUE( clusOut.checkCrossClusters() );
    
    // make an archive
    {
        std::ofstream ofs("ChessboardBiclusteringSerializationTest.xml");
        ASSERT_TRUE( ofs.good() );
        boost::archive::xml_oarchive oa(ofs);
        oa << BOOST_SERIALIZATION_NVP( clusOut );
    }

    // read an archive
    {
        std::ifstream ifs("ChessboardBiclusteringSerializationTest.xml");
        ASSERT_TRUE( ifs.good() );
        boost::archive::xml_iarchive ia(ifs);
        ChessboardBiclustering  clusIn;
        ia >> BOOST_SERIALIZATION_NVP( clusIn );
        
        EXPECT_TRUE( clusIn.checkObjectsPartition() );
        EXPECT_TRUE( clusIn.checkProbesPartition() );
        EXPECT_TRUE( clusIn.checkCrossClusters() );

        EXPECT_TRUE( clusIn.objectsCount() == clusOut.objectsCount() );
        EXPECT_TRUE( clusIn.probesCount() == clusOut.probesCount() );
        EXPECT_TRUE( clusIn.objectsClusters().size() == clusOut.objectsClusters().size() );
        EXPECT_TRUE( clusIn.probesClusters().size() == clusOut.probesClusters().size() );
        EXPECT_TRUE( clusIn.enabledCrossClustersCount() == clusOut.enabledCrossClustersCount() );
        EXPECT_TRUE( clusIn.objectMultiples() == clusOut.objectMultiples() );
    }
}
