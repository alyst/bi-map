#include <set>
#include <map>

#include <dynamic_bitset_utils.h>

#include <BIMAPSampler.h>
#include <ChessboardBiclusteringsPDFEval.h>

#include <ConsolePTCExecutionMonitor.h>

#include "TestCommon.h"
#include "OPADataTest.h"

/*************************************************/

TEST( BIMAPSampler, one_by_one_parameter_estimation )
{
//      msrun
//protein  1-1  1-2  1-3  1-4  1-5  1-6
//      A 1003 1046 1029 1022 1045 1068
//      probe signal = 5, object.multiple = 1
    OPAData    data = loadTestOPAData( "osadata_1x1.xml" );

    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    ChessboardBiclusteringPriors       priors;
    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    GibbsSamplerParams          gibbsParams;
    TurbineCascadeParams        cascadeParams;
    cascadeParams.levelsCount = 5;
    cascadeParams.turbinesCount = 5;
    cascadeParams.turbineParams.particleSamplingPeriod = 20;
    gibbsParams.objectMultipleRate = 0;
    gibbsParams.objectMembershipRate = gibbsParams.objectsSplitMergeRate = 0.0;
    gibbsParams.probeMembershipRate = gibbsParams.probesSplitMergeRate = 0.0;
    gibbsParams.crossClusterFlipRate = gibbsParams.objectMembershipRate = 0.0;
    signalParams.sequenceLengthFactor = 0.5;
    BIMAPSampleCollectorParams collectorParams;
    collectorParams.walkSamples = 200;
    collectorParams.priorsStoragePeriod = 5;

    ChessboardBiclustering clus( data.objectsCount(), data.probesCount() );
    clus.addObjectCluster( 0 );
    clus.addProbeCluster( 0 );
    clus.setObjectMultiple( 0, 1 );
    clus.setCrossCluster( 0, 0, true );
    StdOutPTCExecutionMonitor   mon( 1 );
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    BIMAPSamplerHelper         helper( precomputed,
                                        hyperpriors, priors, gibbsParams );
    BIMAPWalk walk = BIMAPSampler_run( helper, crossClusteringsIndexing, gibbsParams, cascadeParams, 
                                         collectorParams, clus, &mon );
    EXPECT_NO_THROW( walk.check() );
    LOG_INFO( "Clustering steps recorded: " << walk.stepsCount() );
    LOG_INFO( "Prior parameters steps recorded: " << walk.priorParamsStepsCount() );
    RecordProperty( "Clustering steps recorded", walk.stepsCount() );
    RecordProperty( "Prior parameters steps recorded", walk.priorParamsStepsCount() );

    double signalSum = 0.0;
    double signalSqrSum = 0.0;

    for ( BIMAPWalk::const_step_iterator it = walk.stepsBegin(); it != walk.stepsEnd(); ++it ) {
        const ChessboardBiclusteringIndexed& clus = it->clustering;
        ASSERT_EQ( clus.crossClustersData().size(), 1 );
        signal_t signal = clus.crossClustersData().begin()->second;
        signalSum += signal;
        signalSqrSum += signal * signal;
    }
    double signalMean = signalSum / walk.stepsCount();
    double signalVar = signalSqrSum / ( walk.stepsCount() - 1 ) - signalMean * signalMean;
    double signalRelVar = signalVar / signalMean;
    EXPECT_NEAR( signalMean, 5.0, 1 );
    LOG_INFO( "Signal variance: " << signalVar << " (" << signalRelVar * 100 << "%)" );
    RecordProperty( "Signal variance", signalVar );
    RecordProperty( "Signal variance percent", signalRelVar * 100 );
}

struct CluStat {
    typedef boost::dynamic_bitset<> elems_t;
    char   zeroChar;
    elems_t elems;
    size_t cnt;
    size_t enabledCnt;

    static elems_t Elems( bool hasA, bool hasB, bool hasC )
    {
        elems_t elems( 3 );
        elems.set( 0, hasA );
        elems.set( 1, hasB );
        elems.set( 2, hasC );
        return ( elems );
    }

    static elems_t Elems( bool hasA, bool hasB, bool hasC, bool hasD )
    {
        elems_t elems( 4 );
        elems.set( 0, hasA );
        elems.set( 1, hasB );
        elems.set( 2, hasC );
        elems.set( 3, hasD );
        return ( elems );
    }

    CluStat( char zeroChar, bool hasA, bool hasB, bool hasC )
    : zeroChar( zeroChar ), elems( Elems( hasA, hasB, hasC ) )
    , cnt( 0 ), enabledCnt( 0 )
    {
    }

    CluStat( char zeroChar, bool hasA, bool hasB, bool hasC, bool hasD )
    : zeroChar( zeroChar ), elems( Elems( hasA, hasB, hasC, hasD ) )
    , cnt( 0 ), enabledCnt( 0 )
    {
    }

    CluStat( char zeroChar, const elems_t& elems )
    : zeroChar( zeroChar ), elems( elems )
    , cnt( 0 ), enabledCnt( 0 )
    {}

    bool match( const elems_t& elems ) const 
    {
        return ( CluStat::elems == elems );
    }

    void update( const elems_t& elems, bool enabled )
    {
        if ( match( elems ) ) {
            cnt++;
            if ( enabled ) enabledCnt++;
        }
    }

    void update( bool hasA, bool hasB, bool hasC, bool enabled )
    {
        update( Elems( hasA, hasB, hasC ), enabled );
    }

    void update( bool hasA, bool hasB, bool hasC, bool hasD, bool enabled )
    {
        update( Elems( hasA, hasB, hasC, hasD ), enabled );
    }

    friend std::ostream& operator<<( std::ostream& out, const CluStat& stat ) {
        for ( size_t i = 0; i < stat.elems.size(); i++ ) {
            if ( stat.elems.test( i ) ) out << (char)(stat.zeroChar + i);
        }
        out << ":" << stat.enabledCnt << "/" << stat.cnt;
        return ( out );
    }
};

TEST( BIMAPSampler, objects_partition )
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

    for ( object_index_t objIx = 0; objIx < data.objectsCount(); objIx++ ) {
        for ( assay_index_t assayIx = 0; assayIx < data.assaysCount(); assayIx++ ) {
            EXPECT_GT( data.measurement( objIx, assayIx ), 0 );
        }
    }
    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    signalParams.sequenceLengthFactor = 0.5;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    hyperpriors.signalHyperprior.varDistrib.scale = 1;
    ChessboardBiclusteringPriors       priors;
    priors.objectClustering.concentration = 0.1;
    priors.objectClustering.discount = 0.0;
    priors.cellEnablementProb = 0.1;
    GibbsSamplerParams          gibbsParams;
    gibbsParams.objectMembershipRate = 0;
    gibbsParams.objectsSplitMergeRate = 0.1; 
    gibbsParams.crossClusterFlipRate = 0.5;
    gibbsParams.probeMembershipRate = gibbsParams.probesSplitMergeRate = 0;
    gibbsParams.objectMultipleRate = 0;
    TurbineCascadeParams        cascadeParams;
    cascadeParams.levelsCount = 1;
    cascadeParams.turbinesCount = 1;
    cascadeParams.turbineParams.particleSamplingPeriod = 20;
    BIMAPSampleCollectorParams collectorParams;
    collectorParams.walkSamples = 300;
    collectorParams.priorsStoragePeriod = 5;

    StdOutPTCExecutionMonitor   mon( 1 );
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    BIMAPSamplerHelper         helper( precomputed,
                                        hyperpriors, priors, gibbsParams );
    BIMAPWalk walk = BIMAPSampler_run( helper, crossClusteringsIndexing, gibbsParams, cascadeParams,
                                         collectorParams, helper.randomClustering(), &mon );
    RecordProperty( "Clustering steps recorded", walk.stepsCount() );
    RecordProperty( "Prior parameters steps recorded", walk.priorParamsStepsCount() );
    EXPECT_NO_THROW( walk.check() );

    size_t  proprerClusterings = 0;

    std::vector<CluStat> stats;
    for ( int i = 0; i < 8; i++ ) {
        stats.push_back( CluStat( 'A', i & 1, i & 2, i & 4 ) );
    }

    for ( BIMAPWalk::const_step_iterator it = walk.stepsBegin(); it != walk.stepsEnd(); ++it ) {
        const ChessboardBiclusteringIndexed& clus = it->clustering;
        bool  proper = true;
        for ( ChessboardBiclusteringIndexed::objects_cluster_collection_type::const_iterator cluIt = clus.objectsClusters().begin(); cluIt != clus.objectsClusters().end(); ++cluIt ) {
            const object_set_t& objs = (*cluIt)->value();
            bool  hasA = objs.find( 0 ) != objs.end();
            bool  hasB = objs.find( 1 ) != objs.end();
            bool  hasC = objs.find( 2 ) != objs.end();
            bool  isEnabled = clus.isCrossClusterEnabled( (*cluIt)->serial(), (*clus.probesClusters().begin())->serial() );

            for ( size_t i = 0; i < stats.size(); i++ ) {
                stats[i].update( hasA, hasB, hasC, isEnabled );
            }

            if ( ( hasA != hasB ) || ( hasA && hasC ) || ( hasC && hasB ) || !isEnabled ) {
                proper = false;
                //break;
            }
        }
        if ( proper ) {
            proprerClusterings++;
        }
    }
    RecordProperty( "Total clusterings", walk.stepsCount() );
    RecordProperty( "Proper clusterings", proprerClusterings );
    LOG_INFO( "Total clusterings: " << walk.stepsCount() );
    LOG_INFO( "Proper clusterings: " << proprerClusterings );
    for ( size_t i = 0; i < stats.size(); i ++ ) {
        LOG_INFO( stats[i] );
    }
    EXPECT_GT( (double)proprerClusterings / walk.stepsCount(), 0.75 );

    ChessboardBiclusteringsPDFEval   pdfAdj( walk, 0.25, 0.25 );
    EXPECT_EQ( pdfAdj.objectComponents().size(), 2 );
    EXPECT_EQ( pdfAdj.probeComponents().size(), 1 );
}

TEST_F2( GslRngTestF, BIMAPSamplerRng, noise_prediction )
{
    OPAData    data = loadTestOPAData( "osadata_4x2_signal_noise.xml" );

    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    signalParams.sequenceLengthFactor = 0.5;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    ChessboardBiclusteringPriors       priors;
    GibbsSamplerParams          gibbsParams;

    ChessboardBiclustering crossClu( data.objectsCount(), data.probesCount() );
    //crossClu.derivedPriors().clusterProbeSignalSigma = 0.5;
    crossClu.setSignalPrior( GaussianDistribution( 0, 0.5 ) );
    crossClu.setNoiseParams( GeometricDistribution::BySuccessRate( 0.9, 0 ) );
    crossClu.addObjectCluster( 0 );
    crossClu.setObjectCluster( 1, 0 );
    crossClu.addObjectCluster( 2 );
    crossClu.setObjectCluster( 3, 1 );

    crossClu.addProbeCluster( 0 );
    crossClu.setProbeCluster( 1, 0 );

    crossClu.setObjectMultiple( 0, 1 );
    crossClu.setObjectMultiple( 1, 1 );
    crossClu.setObjectMultiple( 2, 1 );
    crossClu.setObjectMultiple( 3, 1 );

    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    ChessboardBiclusteringGibbsSampler   sampler( rndGen, precomputed, hyperpriors, gibbsParams );
    ChessboardBiclusteringFit clusFit( precomputed, priors, crossClu );
    ChessboardBiclusteringGibbsHelper gibbsHelper = sampler.createGibbsHelper( clusFit );

    LOG_INFO( "With signal shape = 0.0" );
    signalParams.scShape = 0.0;
    const size_t trials = 1000;
    size_t  noiseProb = 0;
    size_t  signalProb = 0;

    for ( size_t i = 0; i < trials; i++ ) {
        crossClu.setCrossCluster( 0, 0, gibbsHelper.sampleCrossClusterEnablement( 0, 0 ).value );
        if ( crossClu.isCrossClusterEnabled( 0, 0 ) ) {
            ObjectsClusterParams params = crossClu.objectsClusterParams( 0 );
            if ( params.probeSignal.empty() ) {
                params.probeSignal[0] = gibbsHelper.initialSignal( crossClu.objectsCluster( 0 ).items(),
                                                                   crossClu.probesCluster( 0 ).items(), NULL ).value;
                crossClu.objectsClusterParams( 0 ) = params;
            }
            signalProb++;
        }
        crossClu.setCrossCluster( 1, 0, gibbsHelper.sampleCrossClusterEnablement( 1, 0 ).value );
        if ( crossClu.isCrossClusterEnabled( 1, 0 ) ) {
            ObjectsClusterParams params = crossClu.objectsClusterParams( 1 );
            if ( params.probeSignal.empty() ) {
                params.probeSignal[0] = gibbsHelper.initialSignal( crossClu.objectsCluster( 1 ).items(),
                                                         crossClu.probesCluster( 0 ).items(), NULL ).value;
                crossClu.objectsClusterParams( 1 ) = params;
            }
        }
        else {
            noiseProb++;
        }
    }

    LOG_INFO( "Noise prediction: " << noiseProb << "/" << trials );
    LOG_INFO( "Signal prediction: " << signalProb << "/" << trials );
    EXPECT_GT( (double)noiseProb/trials, 0.9 );
    EXPECT_GT( (double)signalProb/trials, 0.9 );

#if 1 /// @todo: parametrize test
    LOG_INFO( "With signal shape = 0.9" );
    signalParams.scShape = 0.9;
    noiseProb = 0;
    signalProb = 0;

    for ( size_t i = 0; i < trials; i++ ) {
        if ( gibbsHelper.sampleCrossClusterEnablement( 0, 0 ).value ) signalProb++;
        if ( !gibbsHelper.sampleCrossClusterEnablement( 1, 0 ).value ) noiseProb++;
    }

    LOG_INFO( "Noise prediction: " << noiseProb << "/" << trials );
    LOG_INFO( "Signal prediction: " << signalProb << "/" << trials );
    EXPECT_GT( (double)noiseProb/trials, 0.9 );
    EXPECT_GT( (double)signalProb/trials, 0.9 );
#endif
}

TEST( BIMAPSampler, objects_n_probes_partition )
{
/*       msrun
protein 1-1 1-2 1-3 1-4   2-1   2-2   2-3   2-4 3-1 3-2 3-3 3-4   4-1   4-2   4-3   4-4
      A  59  51  45  56 21070 21068 21167 21083   0   0   0   0     0     0     0     0
      B  52  49  52  52 20982 20951 21060 21158   0   0   0   0     0     0     1     0
      C   0   0   0   0     0     0     0     0  39  53  42  45 21015 21036 21116 21015
      D   1   0   0   0     0     0     0     0  51  48  46  53 21056 21088 21082 20977*/
// objects A(2)+B(8), C(2)+D(8)

    OPAData    data = loadTestOPAData( "osadata_4x4_diag.xml" );

    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    //hyperpriors.signalHyperprior.meanVarScale = 20;
    ChessboardBiclusteringPriors       priors;
    //priors.objectClustering.concentration = 0.5;
    //priors.probeClustering.concentration = 0.02;
    //priors.objectClustering.discount = 0.5;
    //priors.probeClustering.discount= 0.001;
    //priors.enabledCC = BetaDistribution( 1000, 1 );
    //priors.noise = BetaDistribution( 1000, 1 );
    priors.cellEnablementProb = 0.01;
    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    signalParams.sequenceLengthFactor = 0.5;
    signalParams.scShape = 0.5;
    GibbsSamplerParams          gibbsParams;
    gibbsParams.objectMultipleRate = 0;
    TurbineCascadeParams        cascadeParams;
    //cascadeParams.maxTemperature = 5;
    cascadeParams.burnInIterations = 12000;
    cascadeParams.ladderAdjustPeriod = 2000;
    cascadeParams.levelsCount = 3;
    cascadeParams.turbinesCount = 5;
    cascadeParams.turbineParams.particleSamplingPeriod = 81;
    BIMAPSampleCollectorParams collectorParams;
    collectorParams.walkSamples = 1000;
    collectorParams.priorsStoragePeriod = 25;

    StdOutPTCExecutionMonitor   mon( 1 );
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    BIMAPSamplerHelper         helper( precomputed,
                                        hyperpriors, priors, gibbsParams );
    BIMAPWalk walk = BIMAPSampler_run( helper, crossClusteringsIndexing, gibbsParams, cascadeParams,
                                         collectorParams, helper.randomClustering(), &mon );
    EXPECT_NO_THROW( walk.check() );
    RecordProperty( "Clustering steps recorded", walk.stepsCount() );
    RecordProperty( "Prior parameters steps recorded", walk.priorParamsStepsCount() );

    std::vector<CluStat> objCluStats;
    for ( int i = 0; i < 16; i++ ) {
        objCluStats.push_back( CluStat( 'A', i & 1, i & 2, i & 4, i & 8 ) );
    }
    std::vector<CluStat> probeCluStats;
    for ( int i = 0; i < 16; i++ ) {
        probeCluStats.push_back( CluStat( '1', i & 1, i & 2, i & 4, i & 8 ) );
    }

    size_t  proprerClusterings = 0;

    double meanNoise = 0;
    double meanSignal = 0;

    typedef std::map<std::pair<object_index_t,probe_index_t>, int>  cc_enabled_counts_map;
    cc_enabled_counts_map   enabledCounts;

    for ( BIMAPWalk::const_step_iterator it = walk.stepsBegin(); it != walk.stepsEnd(); ++it ) {
        const ChessboardBiclusteringIndexed& clus = it->clustering;
        bool  proper = true;
        for ( ChessboardBiclusteringIndexed::objects_cluster_collection_type::const_iterator cluIt = clus.objectsClusters().begin(); cluIt != clus.objectsClusters().end(); ++cluIt ) {
            const object_set_t& objs = (*cluIt)->value();
            bool  hasA = objs.find( 0 ) != objs.end();
            bool  hasB = objs.find( 1 ) != objs.end();
            bool  hasC = objs.find( 2 ) != objs.end();
            bool  hasD = objs.find( 3 ) != objs.end();

            if ( ( hasA != hasB ) || ( hasC != hasD ) 
                || ( hasA == hasC ) || ( hasA == hasD ) 
                || ( hasB == hasC ) || ( hasB == hasD )
            ){
                proper = false;
                //break;
            }
            for ( size_t i = 0; i < probeCluStats.size(); i++ ) {
                objCluStats[ i ].update( hasA, hasB, hasC, hasD, true );
            }
        }
        for ( ChessboardBiclusteringIndexed::probes_cluster_collection_type::const_iterator cluIt = clus.probesClusters().begin(); cluIt != clus.probesClusters().end(); ++cluIt ) {
            const probe_bitset_t& probes = (*cluIt)->value();
            bool  has1 = probes.test( 0 );
            bool  has2 = probes.test( 1 );
            bool  has3 = probes.test( 2 );
            bool  has4 = probes.test( 3 );

            if ( ( has1 != has2 ) || ( has3 != has4 ) 
                || ( has1 == has3 ) || ( has1 == has4 ) 
                || ( has2 == has3 ) || ( has2 == has4 )
            ){
                proper = false;
            }
            for ( size_t i = 0; i < probeCluStats.size(); i++ ) {
                probeCluStats[ i ].update( has1, has2, has3, has4, true );
            }
        }
        for ( ChessboardBiclusteringIndexed::objects_cluster_collection_type::const_iterator oCluIt = clus.objectsClusters().begin(); oCluIt != clus.objectsClusters().end(); ++oCluIt ) {
            const object_set_t& objs = (*oCluIt)->value();
            for ( ChessboardBiclusteringIndexed::probes_cluster_collection_type::const_iterator sCluIt = clus.probesClusters().begin(); sCluIt != clus.probesClusters().end(); ++sCluIt ) {
                const probe_bitset_t& probes = (*sCluIt)->value();
                bool  hasA = objs.find( 0 ) != objs.end();
                bool  hasB = objs.find( 1 ) != objs.end();
                bool  hasC = objs.find( 2 ) != objs.end();
                bool  hasD = objs.find( 3 ) != objs.end();
                bool  has1 = probes.test( 0 );
                bool  has2 = probes.test( 1 );
                bool  has3 = probes.test( 2 );
                bool  has4 = probes.test( 3 );

                if ( clus.isCrossClusterEnabled( (*oCluIt)->serial(), (*sCluIt)->serial() ) ) {
                    if ( ( ( hasA || hasB ) && ( has3 || has4 ) )
                        || ( ( hasC || hasD ) && ( has1 || has2 ) )
                    ){
                        proper = false;
                    }

                    for ( object_set_t::const_iterator objIt = objs.begin(); objIt != objs.end(); ++objIt ) {
                        foreach_bit( probe_index_t, probeIx, probes ) {
                            cc_enabled_counts_map::iterator it = enabledCounts.find( cc_enabled_counts_map::key_type( *objIt, probeIx ) );
                            if ( it != enabledCounts.end() ) {
                                it->second++;
                            }
                            else {
                                enabledCounts.insert( it, cc_enabled_counts_map::value_type( cc_enabled_counts_map::key_type( *objIt, probeIx ), 1 ) );
                            }
                        }
                    }
                }
                else {
                    if ( ( ( hasA || hasB ) && ( has1 || has2 ) )
                        || ( ( hasC || hasD ) && ( has3 || has4 ) )
                    ){
                        proper = false;
                    }
                }
            }
        }
        if ( proper ) proprerClusterings++;
        meanSignal += clus.baselineSignalParams().lnScRate();
        meanNoise += clus.noiseParams().successRate;
    }

    double meanSignalSigma = 0;
    for ( BIMAPWalk::const_priors_step_iterator it = walk.priorParamsStepsBegin(); it != walk.priorParamsStepsEnd(); ++it ) {
        meanSignalSigma += it->priors.signalPrior.sigma;
    }

    for ( cc_enabled_counts_map::const_iterator it = enabledCounts.begin(); it != enabledCounts.end(); ++it ) {
        LOG_INFO( "[" << data.object( it->first.first ).label() << "," << data.probe( it->first.second ).label() << "]=" << it->second << "/" << walk.stepsCount() );
    }

    LOG_INFO( "Objects:" );
    for ( size_t i = 0; i < objCluStats.size(); i++ ) {
        LOG_INFO_IF( objCluStats[ i ].cnt > 0, objCluStats[ i ] );
    }
    LOG_INFO( "Probes:" );
    for ( size_t i = 0; i < probeCluStats.size(); i++ ) {
        LOG_INFO_IF( probeCluStats[ i ].cnt > 0, probeCluStats[ i ] );
    }
    EXPECT_GT( (double)proprerClusterings / walk.stepsCount(), 0.75 );
    LOG_INFO( "Noise: " << meanNoise / walk.stepsCount() );
    LOG_INFO( "Signal: " << meanSignal / walk.stepsCount() );
    LOG_INFO( "Signal sigma: " << meanSignalSigma / walk.priorParamsStepsCount() );

    ChessboardBiclusteringsPDFEval   pdfAdj( walk, 0.25, 0.25 );
    EXPECT_EQ( 2, pdfAdj.objectComponents().size() );
    EXPECT_EQ( 2, pdfAdj.objectComponents()[0].size() );
    if ( pdfAdj.objectComponents().size() >= 2 ) EXPECT_EQ( 2, pdfAdj.objectComponents()[1].size() );
    EXPECT_EQ( 2, pdfAdj.probeComponents().size() );
    EXPECT_EQ( 2, pdfAdj.probeComponents()[0].count() );
    if ( pdfAdj.probeComponents().size() >= 2 ) EXPECT_EQ( 2, pdfAdj.probeComponents()[1].count() );
}

TEST( BIMAPSampler, DISABLED_objects_n_probes_partition_2 )
{
// objects A+B(2,8), C(8)
    OPAData    data = loadTestOPAData( "osadata_4x4_triag_lower.xml" );
//        msrun
// protein 1-1 1-2 1-3 1-4 2-1 2-2 2-3 2-4
//       A  58  45  34  46  10   7   8  13
//       B  47  65  41  47  16   9  14   5
//       C  21  16  14  17  23  20  43  29
//       D  23  19  17  24  24  32  26  28

//        msrun
// protein 3-1 3-2 3-3 3-4 4-1 4-2 4-3 4-4
//       A   0   0   0   0   0   0   0   0
//       B   0   0   0   0   0   0   1   0
//       C  45  48  59  53  37  35  37  36
//       D  72  50  65  54  35  39  34  28

    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    ChessboardBiclusteringPriors       priors;
    priors.objectClustering.concentration = 0.1;
    priors.objectClustering.discount = 0.01;
    priors.probeClustering.concentration = 0.1;
    priors.probeClustering.discount = 0.01;
    priors.cellEnablementProb = 0.1;
    GibbsSamplerParams          gibbsParams;
    //gibbsParams.objectMembershipUpdatePeriod = 3; 
    TurbineCascadeParams        cascadeParams;
    cascadeParams.ladderAdjustPeriod = 1000;
    cascadeParams.turbineParams.particleSamplingPeriod = 40;
    cascadeParams.turbineParams.eeJumpRate = 0.1;
    cascadeParams.turbineParams.generateRate = 0.005;
    cascadeParams.turbinesCount = 1;
    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    signalParams.sequenceLengthFactor = 0.5;
    signalParams.scShape = 0.8;
    BIMAPSampleCollectorParams collectorParams;
    collectorParams.walkSamples = 800;
    collectorParams.priorsStoragePeriod = 5;

    StdOutPTCExecutionMonitor   mon( 1 );
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    BIMAPSamplerHelper         helper( precomputed,
                                        hyperpriors, priors, gibbsParams );
    BIMAPWalk walk = BIMAPSampler_run( helper, crossClusteringsIndexing, gibbsParams, cascadeParams,
                                         collectorParams, helper.randomClustering(), &mon );
    EXPECT_NO_THROW( walk.check() );
    RecordProperty( "Clustering steps recorded", walk.stepsCount() );
    RecordProperty( "Prior parameters steps recorded", walk.priorParamsStepsCount() );

    std::vector<CluStat> probeCluStats;
    for ( int i = 0; i < 16; i++ ) {
        probeCluStats.push_back( CluStat( '1', i & 1, i & 2, i & 4, i & 8 ) );
    }

    size_t  proprerClusterings = 0;

    size_t  abCount = 0;
    size_t  cdCount = 0;

    double meanNoise = 0;
    double meanSignal = 0;

    typedef std::map<std::pair<object_index_t,probe_index_t>, int>  cc_enabled_counts_map;
    cc_enabled_counts_map   enabledCounts;

    for ( BIMAPWalk::const_step_iterator it = walk.stepsBegin(); it != walk.stepsEnd(); ++it ) {
        const ChessboardBiclusteringIndexed& clus = it->clustering;
        bool  proper = true;
        for ( ChessboardBiclusteringIndexed::objects_cluster_collection_type::const_iterator cluIt = clus.objectsClusters().begin(); cluIt != clus.objectsClusters().end(); ++cluIt ) {
            const object_set_t& objs = (*cluIt)->value();
            bool  hasA = objs.find( 0 ) != objs.end();
            bool  hasB = objs.find( 1 ) != objs.end();
            bool  hasC = objs.find( 2 ) != objs.end();
            bool  hasD = objs.find( 3 ) != objs.end();

            if ( ( hasA != hasB ) || ( hasC != hasD ) 
                || ( hasA == hasC ) || ( hasA == hasD ) 
                || ( hasB == hasC ) || ( hasB == hasD )
            ){
                proper = false;
                //break;
            }
            if ( hasA && hasB && !hasC && !hasD ) {
                abCount++;
            }
            if ( !hasA && !hasB && hasC && hasD ) {
                cdCount++;
            }
        }
        for ( ChessboardBiclusteringIndexed::probes_cluster_collection_type::const_iterator cluIt = clus.probesClusters().begin(); cluIt != clus.probesClusters().end(); ++cluIt ) {
            const probe_bitset_t& probes = (*cluIt)->value();
            bool  has1 = probes.test( 0 );
            bool  has2 = probes.test( 1 );
            bool  has3 = probes.test( 2 );
            bool  has4 = probes.test( 3 );
            
            if ( ( has1 != has2 ) || ( has3 != has4 ) 
                || ( has1 == has3 ) || ( has1 == has4 ) 
                || ( has2 == has3 ) || ( has2 == has4 )
            ){
                proper = false;
            }
            for ( size_t i = 0; i < probeCluStats.size(); i++ ) {
                probeCluStats[ i ].update( has1, has2, has3, has4, true );
            }
        }
        if ( proper ) proprerClusterings++;
        for ( ChessboardBiclusteringIndexed::objects_cluster_collection_type::const_iterator oCluIt = clus.objectsClusters().begin(); oCluIt != clus.objectsClusters().end(); ++oCluIt ) {
            const object_set_t& objs = (*oCluIt)->value();
            for ( ChessboardBiclusteringIndexed::probes_cluster_collection_type::const_iterator sCluIt = clus.probesClusters().begin(); sCluIt != clus.probesClusters().end(); ++sCluIt ) {
                const probe_bitset_t& probes = (*sCluIt)->value();
                bool  hasA = objs.find( 0 ) != objs.end();
                bool  hasB = objs.find( 1 ) != objs.end();
                bool  hasC = objs.find( 2 ) != objs.end();
                bool  hasD = objs.find( 3 ) != objs.end();
                bool  has1 = probes.test( 0 );
                bool  has2 = probes.test( 1 );
                bool  has3 = probes.test( 2 );
                bool  has4 = probes.test( 3 );

                if ( clus.isCrossClusterEnabled( (*oCluIt)->serial(), (*sCluIt)->serial() ) ) {
                    if ( ( ( hasA || hasB ) && ( has3 || has4 ) )
                        && ( ( hasC || hasD ) && ( has1 || has2 ) )
                    ){
                        proper = false;
                    }

                    for ( object_set_t::const_iterator objIt = objs.begin(); objIt != objs.end(); ++objIt ) {
                        foreach_bit( probe_index_t, probeIx, probes ) {
                            cc_enabled_counts_map::iterator it = enabledCounts.find( cc_enabled_counts_map::key_type( *objIt, probeIx ) );
                            if ( it != enabledCounts.end() ) {
                                it->second++;
                            }
                            else {
                                enabledCounts.insert( it, cc_enabled_counts_map::value_type( cc_enabled_counts_map::key_type( *objIt, probeIx ), 1 ) );
                            }
                        }
                    }
                }
                else {
                    if ( ( ( hasA || hasB ) && ( has1 || has2 ) )
                        && ( ( hasC || hasD ) && ( has3 || has4 ) )
                    ){
                        proper = false;
                    }
                }
            }
        }
        meanSignal += clus.baselineSignalParams().lnScRate();
        meanNoise += clus.noiseParams().successRate;
    }

    double meanSignalSigma = 0;
    for ( BIMAPWalk::const_priors_step_iterator it = walk.priorParamsStepsBegin(); it != walk.priorParamsStepsEnd(); ++it ) {
        meanSignalSigma += it->priors.signalPrior.sigma;
    }

    for ( cc_enabled_counts_map::const_iterator it = enabledCounts.begin(); it != enabledCounts.end(); ++it ) {
        LOG_INFO( "[" << data.object( it->first.first ).label() << "," << data.probe( it->first.second ).label() << "]=" << it->second << "/" << walk.stepsCount() );
    }

    LOG_INFO( "AB count: " << abCount );
    LOG_INFO( "CD count: " << cdCount );
    LOG_INFO( "Probes:" );
    for ( size_t i = 0; i < probeCluStats.size(); i++ ) {
        LOG_INFO( probeCluStats[ i ] );
    }
    EXPECT_GT( (double)proprerClusterings / walk.stepsCount(), 0.75 );
    LOG_INFO( "Noise: " << meanNoise / walk.stepsCount() );
    LOG_INFO( "Signal: " << meanSignal / walk.stepsCount() );
    LOG_INFO( "Signal sigma: " << meanSignalSigma / walk.priorParamsStepsCount() );
}

TEST( BIMAPSampler, DISABLED_sampling_test_run )
{
    OPAData                     data = generateTestOPAData();
    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    ChessboardBiclusteringPriors       priors;
    GibbsSamplerParams          gibbsParams;
    TurbineCascadeParams        cascadeParams;
    cascadeParams.turbineParams.particleSamplingPeriod = 20;
    BIMAPSampleCollectorParams collectorParams;
    collectorParams.walkSamples = 300;
    collectorParams.priorsStoragePeriod = 5;

    StdOutPTCExecutionMonitor   mon( 1 );
    PrecomputedData    precomputed( data, precomputedDataParams, signalParams );
    BIMAPSamplerHelper         helper( precomputed,
                                        hyperpriors, priors, gibbsParams );
    BIMAPWalk walk = BIMAPSampler_run( helper, crossClusteringsIndexing, gibbsParams, cascadeParams,
                                         collectorParams, helper.randomClustering(), &mon );
    EXPECT_NO_THROW( walk.check() );
    RecordProperty( "Clustering steps recorded", walk.stepsCount() );
    RecordProperty( "Prior parameters steps recorded", walk.priorParamsStepsCount() );
}
