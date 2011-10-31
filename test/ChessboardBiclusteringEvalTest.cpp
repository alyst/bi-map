#include <set>
#include <map>

#include <dynamic_bitset_utils.h>

#include <ChessboardBiclusteringFit.h>

#include "TestCommon.h"
#include "OPADataTest.h"

#include "ObjectsPartition.h"
#include "ProbesPartition.h"

/*************************************************/

TEST( ChessboardBiclusteringEval, objects_partition )
{
//    msrun prey_ac    sc
// 1    1-1       A    59
// 2    1-1       B    57
// 3    1-1       C 20938
// 4    1-2       A    51
// 5    1-2       B    63
// 6    1-2       C 21097
// 7    1-3       A    45
// 8    1-3       B    55
// 9    1-3       C 21095
// 10   1-4       A    56
// 11   1-4       B    52
// 12   1-4       C 21051
// 13   1-5       A    43
// 14   1-5       B    53
// 15   1-5       C 21061
// 16   1-6       A    51
// 17   1-6       B    58
// 18   1-6       C 20997
// objects A+B(2), C(8)

    OPAData    data = loadTestOPAData( "osadata_3x1_2_clusters.xml" );

    PrecomputedDataParams   precomputedDataParams;
    CellSignalParams signalParams;
    signalParams.sequenceLengthFactor = 0.5;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    hyperpriors.signalHyperprior.varDistrib.scale = 1;
    ChessboardBiclusteringPriors       priors;
    priors.objectClustering.concentration = 0.1;
    priors.objectClustering.discount = 0.0;
    priors.cellEnablementProb = 0.1;

    ChessboardBiclustering ccActual( 3, 1 );
    ccActual.setSignalPrior( GaussianDistribution( 0, 0.5 ) );
    //ccActual.setBaselineShape( priors.signalShape );
    ccActual.setNoiseParams( GeometricDistribution::ByFailureRate( 0.1, 0 ) );
    {
        ccActual.addProbeCluster( 0 );
        ccActual.addObjectCluster( 0 );
        ccActual.setObjectCluster( 1, 0 );
        ccActual.addObjectCluster( 2 );
        ObjectsClusterParams oc0Params = ccActual.objectsClusterParams( 0 );
        oc0Params.blocksMask.set( 0 );
        oc0Params.objectMultiple[ 0 ] = 1;
        oc0Params.objectMultiple[ 1 ] = 1;
        oc0Params.probeSignal[ 0 ] = 2;
        ccActual.objectsClusterParams( 0 ) = oc0Params;
        ObjectsClusterParams oc1Params = ccActual.objectsClusterParams( 1 );
        oc1Params.blocksMask.set( 0 );
        oc1Params.objectMultiple[ 2 ] = 1;
        oc1Params.probeSignal[ 0 ] = 8;
        ccActual.objectsClusterParams( 1 ) = oc1Params;
        ccActual.cleanupClusters();
    }
    BOOST_ASSERT( ccActual.checkObjectsPartition() );
    BOOST_ASSERT( ccActual.checkProbesPartition() );
    BOOST_ASSERT( ccActual.checkBlocks() );
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    ChessboardBiclusteringFit fitCcActual( precomputed, priors, ccActual );
    LOG_INFO( "Actual partition llh=" << fitCcActual.llh() );

    ChessboardBiclustering ccWrongSignal( ccActual );
    {
        ObjectsClusterParams oc0Params = ccWrongSignal.objectsClusterParams( 0 );
        oc0Params.probeSignal[ 0 ] = 4;
        ccWrongSignal.objectsClusterParams( 0 ) = oc0Params;
        ObjectsClusterParams oc1Params = ccWrongSignal.objectsClusterParams( 1 );
        oc1Params.probeSignal[ 0 ] = 7;
        ccWrongSignal.objectsClusterParams( 1 ) = oc1Params;
    }
    ChessboardBiclusteringFit fitCcWrongSignal( precomputed, priors, ccWrongSignal );
    LOG_INFO( "Wrong signal partition llh=" << fitCcWrongSignal.llh() );
    EXPECT_GT( fitCcActual.llh(), fitCcWrongSignal.llh() );

    ChessboardBiclustering ccWrongMultiple( ccActual );
    {
        ObjectsClusterParams oc1Params = ccWrongMultiple.objectsClusterParams( 1 );
        oc1Params.objectMultiple[ 2 ] = 2;
        oc1Params.probeSignal[ 0 ] = 7;
        ccWrongMultiple.objectsClusterParams( 1 ) = oc1Params;
    }
    ChessboardBiclusteringFit fitCcWrongMultiple( precomputed, priors, ccWrongMultiple );
    LOG_INFO( "Wrong multiple partition llh=" << fitCcWrongMultiple.llh() );
    EXPECT_GT( fitCcActual.llh(), fitCcWrongMultiple.llh() );

    ChessboardBiclustering ccOneCluster( ccActual );
    {
        ccOneCluster.setObjectCluster( 2, 0 );
        ccOneCluster.cleanupClusters();
        ObjectsClusterParams oc0Params = ccOneCluster.objectsClusterParams( 0 );
        oc0Params.probeSignal[ 0 ] = 4;
        ccOneCluster.objectsClusterParams( 0 ) = oc0Params;
    }
    ChessboardBiclusteringFit fitCcOneCluster( precomputed, priors, ccOneCluster );
    LOG_INFO( "One cluster partition llh=" << fitCcOneCluster.llh() );
    EXPECT_GT( fitCcActual.llh(), fitCcOneCluster.llh() );

    ChessboardBiclustering ccThreeClusters( ccActual );
    {
        ccThreeClusters.addObjectCluster( 1 );
        ObjectsClusterParams oc2Params = ccThreeClusters.objectsClusterParams( 2 );
        oc2Params.blocksMask.set( 0 );
        oc2Params.probeSignal[ 0 ] = 2;
        ccThreeClusters.objectsClusterParams( 2 ) = oc2Params;
        ccThreeClusters.cleanupClusters();
    }
    ChessboardBiclusteringFit fitCcThreeClusters( precomputed, priors, ccThreeClusters );
    LOG_INFO( "3 clusters partition llh=" << fitCcThreeClusters.llh() );
    EXPECT_GE( fitCcActual.llh(), fitCcThreeClusters.llh() );
    EXPECT_GT( fitCcActual.lpp(), fitCcThreeClusters.lpp() );

    ChessboardBiclustering ccTwoWrongClusters( ccActual );
    {
        ccTwoWrongClusters.setObjectCluster( 1, 1 );
        ccTwoWrongClusters.cleanupClusters();
    }
    ChessboardBiclusteringFit fitCcTwoWrongClusters( precomputed, priors, ccTwoWrongClusters );
    LOG_INFO( "2 wrong clusters partition llh=" << fitCcTwoWrongClusters.llh() );
    EXPECT_GT( fitCcActual.llh(), fitCcTwoWrongClusters.llh() );
}

TEST( ChessboardBiclusteringEval, objects_partition_llh_delta )
{
//    msrun prey_ac    sc
// 1    1-1       A    59
// 2    1-1       B    57
// 3    1-1       C 20938
// 4    1-2       A    51
// 5    1-2       B    63
// 6    1-2       C 21097
// 7    1-3       A    45
// 8    1-3       B    55
// 9    1-3       C 21095
// 10   1-4       A    56
// 11   1-4       B    52
// 12   1-4       C 21051
// 13   1-5       A    43
// 14   1-5       B    53
// 15   1-5       C 21061
// 16   1-6       A    51
// 17   1-6       B    58
// 18   1-6       C 20997
// objects A+B(2), C(8)

    OPAData    data = loadTestOPAData( "osadata_3x1_2_clusters.xml" );

    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    signalParams.sequenceLengthFactor = 0.5;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    hyperpriors.signalHyperprior.varDistrib.scale = 1;
    ChessboardBiclusteringPriors       priors;
    priors.objectClustering.concentration = 0.1;
    priors.objectClustering.discount = 0.0;
    priors.cellEnablementProb = 0.1;

    ChessboardBiclustering ccActual( 3, 1 );
    ccActual.setSignalPrior( GaussianDistribution( 0, 0.5 ) );
    ccActual.setNoiseParams( GeometricDistribution::ByFailureRate( 0.1, 0 ) );
     {
        ccActual.addProbeCluster( 0 );
        ccActual.addObjectCluster( 0 );
        ccActual.setObjectCluster( 1, 0 );
        ccActual.addObjectCluster( 2 );
        ObjectsClusterParams oc0Params = ccActual.objectsClusterParams( 0 );
        oc0Params.blocksMask.set( 0 );
        oc0Params.objectMultiple[ 0 ] = 1;
        oc0Params.objectMultiple[ 1 ] = 1;
        oc0Params.probeSignal[ 0 ] = 2;
        ccActual.objectsClusterParams( 0 ) = oc0Params;
        ObjectsClusterParams oc1Params = ccActual.objectsClusterParams( 1 );
        oc1Params.blocksMask.set( 0 );
        oc1Params.objectMultiple[ 2 ] = 1;
        oc1Params.probeSignal[ 0 ] = 8;
        ccActual.objectsClusterParams( 1 ) = oc1Params;
        ccActual.cleanupClusters();
    }
    BOOST_ASSERT( ccActual.checkObjectsPartition() );
    BOOST_ASSERT( ccActual.checkProbesPartition() );
    BOOST_ASSERT( ccActual.checkBlocks() );
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    ChessboardBiclusteringFit fitCcActual( precomputed, priors, ccActual );
    LOG_INFO( "Actual partition llh=" << fitCcActual.llh() );

    ChessboardBiclustering ccOneCluster( ccActual );
    {
        ccOneCluster.setObjectCluster( 2, 0 );
        ccOneCluster.cleanupClusters();
        ObjectsClusterParams oc0Params = ccOneCluster.objectsClusterParams( 0 );
        oc0Params.probeSignal[ 0 ] = 4;
        ccOneCluster.objectsClusterParams( 0 ) = oc0Params;
    }
    ChessboardBiclusteringFit fitCcOneCluster( precomputed, priors, ccOneCluster );
    LOG_INFO( "One cluster partition llh=" << fitCcOneCluster.llh() );

    ChessboardBiclustering ccThreeClusters( ccActual );
    {
        ccThreeClusters.addObjectCluster( 1 );
        ObjectsClusterParams oc2Params = ccThreeClusters.objectsClusterParams( 2 );
        oc2Params.blocksMask.set( 0 );
        oc2Params.probeSignal[ 0 ] = 2;
        ccThreeClusters.objectsClusterParams( 2 ) = oc2Params;
        ccThreeClusters.cleanupClusters();
    }
    ChessboardBiclusteringFit fitCcThreeClusters( precomputed, priors, ccThreeClusters );
    LOG_INFO( "3 clusters partition llh=" << fitCcThreeClusters.llh() );

    ObjectsPartition opActual( fitCcActual );
    FixedObjectsPartitionStats opsActual( fitCcActual );

    std::vector<ObjectsPartition::elements_set_proxy_type> actualToOneItems;
    actualToOneItems.push_back( opActual.cluster( 0 ).items() );
    actualToOneItems[ 0 ] += 2;
    std::vector<ObjectsPartition::cluster_params_type> actualToOneParams;
    actualToOneParams.push_back( opActual.cluster( 0 ).params() );
    actualToOneParams[ 0 ].objectMultiple[ 2 ] = 1;
    actualToOneParams[ 0 ].probeSignal[ 0 ] = 4;
    ObjectsPartition::cluster_index_set_type actualToOneOldClusters;
    actualToOneOldClusters.insert( 0 );
    actualToOneOldClusters.insert( 1 );
    double actualToOneLLHDelta = opsActual.llhDelta( actualToOneItems, actualToOneParams, actualToOneOldClusters );
    LOG_INFO( "Delta between Actual and One cluster partition llhDelta=" << actualToOneLLHDelta );
    EXPECT_DOUBLE_EQ( fitCcOneCluster.llh() - opsActual.llh(), actualToOneLLHDelta );

    std::vector<ObjectsPartition::elements_set_proxy_type> actualToThreeItems;
    actualToThreeItems.push_back( opActual.cluster( 0 ).items() );
    actualToThreeItems.push_back( opActual.cluster( 0 ).items() );
    actualToThreeItems[ 0 ] -= 1;
    actualToThreeItems[ 1 ] -= 0;
    std::vector<ObjectsPartition::cluster_params_type> actualToThreeParams;
    actualToThreeParams.push_back( opActual.cluster( 0 ).params() );
    actualToThreeParams.push_back( opActual.cluster( 0 ).params() );
    actualToThreeParams[ 1 ].probeSignal[ 0 ] = 2;
    ObjectsPartition::cluster_index_set_type actualToThreeOldClusters;
    actualToThreeOldClusters.insert( 0 );
    double actualToThreeDelta = opsActual.llhDelta( actualToThreeItems, actualToThreeParams, actualToThreeOldClusters );
    LOG_INFO( "Delta between Actual and One cluster partition llhDelta=" << actualToThreeDelta );
    EXPECT_DOUBLE_EQ( fitCcThreeClusters.llh() - opsActual.llh(), actualToThreeDelta );
}

TEST( ChessboardBiclusteringEval, objects2_probes2_partition )
{
/*       msrun
protein 1-1 1-2 1-3 1-4   2-1   2-2   2-3   2-4 3-1 3-2 3-3 3-4   4-1   4-2   4-3   4-4
      A  59  51  45  56 21070 21068 21167 21083   0   0   0   0     0     0     0     0
      B  52  49  52  52 20982 20951 21060 21158   0   0   0   0     0     0     1     0
      C   0   0   0   0     0     0     0     0  39  53  42  45 21015 21036 21116 21015
      D   1   0   0   0     0     0     0     0  51  48  46  53 21056 21088 21082 20977*/
// objects A(2)+B(8), C(2)+D(8)

    OPAData    data = loadTestOPAData( "osadata_4x4_diag.xml" );

    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    signalParams.sequenceLengthFactor = 0.5;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    hyperpriors.signalHyperprior.varDistrib.scale = 1;
    ChessboardBiclusteringPriors       priors;
    priors.objectClustering.concentration = 0.1;
    priors.objectClustering.discount = 0.0;
    priors.cellEnablementProb = 0.1;

    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );

    ChessboardBiclustering ccActual( 4, 4 );
    ccActual.setSignalPrior( GaussianDistribution( 0, 0.5 ) );
    ccActual.setNoiseParams( GeometricDistribution::ByFailureRate( 0.1, 0 ) );
     {
        ccActual.addProbeCluster( 0 );
        ccActual.setProbeCluster( 1, 0 );
        ccActual.addProbeCluster( 2 );
        ccActual.setProbeCluster( 3, 1 );
        ccActual.addObjectCluster( 0 );
        ccActual.setObjectCluster( 1, 0 );
        ccActual.addObjectCluster( 2 );
        ccActual.setObjectCluster( 3, 1 );
        ObjectsClusterParams oc0Params = ccActual.objectsClusterParams( 0 );
        oc0Params.blocksMask.set( 0 );
        oc0Params.blocksMask.set( 1, false );
        oc0Params.objectMultiple[ 0 ] = 1;
        oc0Params.objectMultiple[ 1 ] = 1;
        oc0Params.probeSignal[ 0 ] = 2;
        ccActual.objectsClusterParams( 0 ) = oc0Params;
        ObjectsClusterParams oc1Params = ccActual.objectsClusterParams( 1 );
        oc1Params.blocksMask.set( 0, false );
        oc1Params.blocksMask.set( 1 );
        oc1Params.objectMultiple[ 2 ] = 1;
        oc1Params.objectMultiple[ 3 ] = 1;
        oc1Params.probeSignal[ 1 ] = 2;
        ccActual.objectsClusterParams( 1 ) = oc1Params;
    }
    BOOST_ASSERT( ccActual.checkObjectsPartition() );
    BOOST_ASSERT( ccActual.checkProbesPartition() );
    BOOST_ASSERT( ccActual.checkBlocks() );
    ChessboardBiclusteringFit fitCcActual( precomputed, priors, ccActual );
    LOG_INFO( "Actual partition llh=" << fitCcActual.llh() );

    ChessboardBiclustering ccWrongSignal( ccActual );
    {
        ObjectsClusterParams oc0Params = ccWrongSignal.objectsClusterParams( 0 );
        oc0Params.probeSignal[ 0 ] = 7;
        ccWrongSignal.objectsClusterParams( 0 ) = oc0Params;
    }
    ChessboardBiclusteringFit fitCcWrongSignal( precomputed, priors, ccWrongSignal );
    LOG_INFO( "Wrong signal partition llh=" << fitCcWrongSignal.llh() );
    EXPECT_GT( fitCcActual.llh(), fitCcWrongSignal.llh() );

    ChessboardBiclustering ccWrongMultiple( ccActual );
    {
        ObjectsClusterParams oc1Params = ccWrongMultiple.objectsClusterParams( 1 );
        oc1Params.objectMultiple[ 0 ] = 2;
        oc1Params.objectMultiple[ 1 ] = 2;
        oc1Params.probeSignal[ 0 ] = 1;
        oc1Params.probeSignal[ 1 ] = 7;
        ccWrongMultiple.objectsClusterParams( 1 ) = oc1Params;
    }
    ChessboardBiclusteringFit fitCcWrongMultiple( precomputed, priors, ccWrongMultiple );
    LOG_INFO( "Wrong multiple partition llh=" << fitCcWrongMultiple.llh() );
    EXPECT_GT( fitCcActual.llh(), fitCcWrongMultiple.llh() );

    ChessboardBiclustering ccThreeObjectClusters( ccActual );
    {
        ccThreeObjectClusters.addObjectCluster( 3 );
        ObjectsClusterParams oc2Params = ccThreeObjectClusters.objectsClusterParams( 2 );
        oc2Params.blocksMask.set( 0, false );
        oc2Params.blocksMask.set( 1, true );
        oc2Params.probeSignal[ 1 ] = 8;
        ccThreeObjectClusters.objectsClusterParams( 2 ) = oc2Params;
    }
    ChessboardBiclusteringFit fitCcThreeObjectClusters( precomputed, priors, ccThreeObjectClusters );
    LOG_INFO( "Wrong 3 objects clusters partition llh=" << fitCcThreeObjectClusters.llh() );
    EXPECT_GT( fitCcActual.llh(), fitCcThreeObjectClusters.llh() );

    ChessboardBiclustering ccThreeProbeClusters( ccActual );
    {
        ccThreeProbeClusters.addProbeCluster( 3 );
        ProbesClusterParams oc2Params = ccThreeProbeClusters.probesClusterParams( 2 );
        oc2Params.blocksMask.set( 0, false );
        oc2Params.blocksMask.set( 1, true );
        oc2Params.objectsSignal[ 1 ] = 8;
        ccThreeProbeClusters.probesClusterParams( 2 ) = oc2Params;
    }
    ChessboardBiclusteringFit fitCcThreeProbeClusters( precomputed, priors, ccThreeProbeClusters );
    LOG_INFO( "Wrong 3 probes clusters partition llh=" << fitCcThreeProbeClusters.llh() );
    EXPECT_GT( fitCcActual.llh(), fitCcThreeProbeClusters.llh() );

    ChessboardBiclustering ccThreeProbeClusters2( ccThreeProbeClusters );
    {
        ProbesClusterParams oc2Params = ccThreeProbeClusters2.probesClusterParams( 2 );
        oc2Params.blocksMask.set( 0, true );
        oc2Params.blocksMask.set( 1, true );
        oc2Params.objectsSignal[ 0 ] = -3;
        ccThreeProbeClusters2.probesClusterParams( 2 ) = oc2Params;
    }
    ChessboardBiclusteringFit fitCcThreeProbeClusters2( precomputed, priors, ccThreeProbeClusters2 );
    LOG_INFO( "Wrong 3 probes clusters partition-2 llh=" << fitCcThreeProbeClusters2.llh() );
    EXPECT_GT( fitCcActual.llh(), fitCcThreeProbeClusters2.llh() );
}
