#include <gtest/gtest.h>

#include <fstream>
#include <boost/archive/tmpdir.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <OPAData.h>
#include <BIMAPWalk.h>
#include <BIMAPSampler.h>

#include <ConsolePTCExecutionMonitor.h>

TEST( BIMAPWalk, serialize_write )
{
/*       msrun
protein 1-1 1-2 1-3 1-4   2-1   2-2   2-3   2-4 3-1 3-2 3-3 3-4   4-1   4-2   4-3   4-4
      A  59  51  45  56 21070 21068 21167 21083   0   0   0   0     0     0     0     0
      B  52  49  52  52 20982 20951 21060 21158   0   0   0   0     0     0     1     0
      C   0   0   0   0     0     0     0     0  39  53  42  45 21015 21036 21116 21015
      D   1   0   0   0     0     0     0     0  51  48  46  53 21056 21088 21082 20977*/
// objects A+B(2,8), C(8)

    OPAData    data;
    data.addObject( "A", 50 );
    data.addObject( "B", 50 );
    data.addObject( "C", 50 );
    data.addObject( "D", 50 );

    data.addProbe( "1", "A" );
    data.addProbe( "2", "B" );
    data.addProbe( "3", "C" );
    data.addProbe( "4", "D" );

    data.addAssay( "1-1", "1" );
    data.addAssay( "1-2", "1" );
    data.addAssay( "1-3", "1" );
    data.addAssay( "1-4", "1" );
    data.addAssay( "2-1", "2" );
    data.addAssay( "2-2", "2" );
    data.addAssay( "2-3", "2" );
    data.addAssay( "2-4", "2" );
    data.addAssay( "3-1", "3" );
    data.addAssay( "3-2", "3" );
    data.addAssay( "3-3", "3" );
    data.addAssay( "3-4", "3" );
    data.addAssay( "4-1", "4" );
    data.addAssay( "4-2", "4" );
    data.addAssay( "4-3", "4" );
    data.addAssay( "4-4", "4" );

    data.addMeasurement( "A", "1-1", 59 );
    data.addMeasurement( "B", "1-1", 52 );
    data.addMeasurement( "C", "1-1", 0 );
    data.addMeasurement( "D", "1-1", 1 );
    data.addMeasurement( "A", "1-2", 51 );
    data.addMeasurement( "B", "1-2", 49 );
    data.addMeasurement( "C", "1-2", 0 );
    data.addMeasurement( "D", "1-2", 0 );
    data.addMeasurement( "A", "1-3", 45 );
    data.addMeasurement( "B", "1-3", 52 );
    data.addMeasurement( "C", "1-3", 0 );
    data.addMeasurement( "D", "1-3", 0 );
    data.addMeasurement( "A", "1-4", 56 );
    data.addMeasurement( "B", "1-4", 52 );
    data.addMeasurement( "C", "1-4", 0 );
    data.addMeasurement( "D", "1-4", 0 );

    data.addMeasurement( "A", "2-1", 21070 );
    data.addMeasurement( "B", "2-1", 20982 );
    data.addMeasurement( "C", "2-1", 0 );
    data.addMeasurement( "D", "2-1", 0 );
    data.addMeasurement( "A", "2-2", 21068 );
    data.addMeasurement( "B", "2-2", 20951 );
    data.addMeasurement( "C", "2-2", 0 );
    data.addMeasurement( "D", "2-2", 0 );
    data.addMeasurement( "A", "2-3", 21167 );
    data.addMeasurement( "B", "2-3", 21060 );
    data.addMeasurement( "C", "2-3", 0 );
    data.addMeasurement( "D", "2-3", 0 );
    data.addMeasurement( "A", "2-4", 21083 );
    data.addMeasurement( "B", "2-4", 21158 );
    data.addMeasurement( "C", "2-4", 0 );
    data.addMeasurement( "D", "2-4", 0 );

    data.addMeasurement( "A", "3-1", 0 );
    data.addMeasurement( "B", "3-1", 0 );
    data.addMeasurement( "C", "3-1", 39 );
    data.addMeasurement( "D", "3-1", 51 );
    data.addMeasurement( "A", "3-2", 0 );
    data.addMeasurement( "B", "3-2", 0 );
    data.addMeasurement( "C", "3-2", 53 );
    data.addMeasurement( "D", "3-2", 48 );
    data.addMeasurement( "A", "3-3", 0 );
    data.addMeasurement( "B", "3-3", 0 );
    data.addMeasurement( "C", "3-3", 42 );
    data.addMeasurement( "D", "3-3", 46 );
    data.addMeasurement( "A", "3-4", 0 );
    data.addMeasurement( "B", "3-4", 0 );
    data.addMeasurement( "C", "3-4", 45 );
    data.addMeasurement( "D", "3-4", 53 );

    data.addMeasurement( "A", "4-1", 0 );
    data.addMeasurement( "B", "4-1", 0 );
    data.addMeasurement( "C", "4-1", 21015 );
    data.addMeasurement( "D", "4-1", 21056 );
    data.addMeasurement( "A", "4-2", 0 );
    data.addMeasurement( "B", "4-2", 0 );
    data.addMeasurement( "C", "4-2", 21036 );
    data.addMeasurement( "D", "4-2", 21088 );
    data.addMeasurement( "A", "4-3", 0 );
    data.addMeasurement( "B", "4-3", 0 );
    data.addMeasurement( "C", "4-3", 21116 );
    data.addMeasurement( "D", "4-3", 21082 );
    data.addMeasurement( "A", "4-4", 0 );
    data.addMeasurement( "B", "4-4", 0 );
    data.addMeasurement( "C", "4-4", 21015 );
    data.addMeasurement( "D", "4-4", 20977 );

    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    ChessboardBiclusteringHyperPriors  hyperpriors;
    hyperpriors.signalHyperprior.meanVarScale = 2;
    ChessboardBiclusteringPriors       priors;
    //params.eeSamplerParams.maxTemperature = 5;
    //priors.probeClustering.concentration = 0.01;
    priors.cellEnablementProb = 0.1;
    PrecomputedDataParams       precomputedDataParams;
    CellSignalParams            signalParams;
    signalParams.sequenceLengthFactor = 0.5;
    GibbsSamplerParams          gibbsParams;
    TurbineCascadeParams        cascadeParams;
    cascadeParams.ladderAdjustPeriod = 40;
    cascadeParams.turbinesCount = 2;
    BIMAPSampleCollectorParams collectorParams;
    collectorParams.walkSamples = 50;
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

    // make an archive
    std::ofstream ofs("BIMAPWalkTest.xml");
    ASSERT_TRUE( ofs.good() );
    boost::archive::xml_oarchive oa(ofs);
    boost::unordered_map<object_index_t, object_label_t> objects;
    boost::unordered_map<probe_index_t, probe_label_t> probes;
    for ( object_index_t i = 0; i < data.objectsCount(); i++ ) {
        objects[ i ] = data.object( i ).label();
    }
    for ( probe_index_t i = 0; i < data.probesCount(); i++ ) {
        probes[ i ] = data.probe( i ).label();
    }
    oa << BOOST_SERIALIZATION_NVP( objects );
    oa << BOOST_SERIALIZATION_NVP( probes );
    oa << BOOST_SERIALIZATION_NVP( crossClusteringsIndexing );
    oa << BOOST_SERIALIZATION_NVP( walk );
}

TEST( BIMAPWalk, serialize_read )
{
    boost::unordered_map<object_index_t, object_label_t> objects;
    boost::unordered_map<probe_index_t, probe_label_t> probes;
    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    BIMAPWalk walk( crossClusteringsIndexing );
    std::ifstream ifs( "BIMAPWalkTest.xml" );//"/home/astukalov/test2walk.xml" );//);
    ASSERT_TRUE( ifs.good() );
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP( objects );
    LOG_INFO( "Objects count: " << objects.size() );
    ia >> BOOST_SERIALIZATION_NVP( probes );
    LOG_INFO( "Probes count: " << probes.size() );
    ia >> BOOST_SERIALIZATION_NVP( crossClusteringsIndexing );
    LOG_INFO( "Indexing size: " << crossClusteringsIndexing.size() );
    ia >> BOOST_SERIALIZATION_NVP( walk );
    EXPECT_NO_THROW( walk.check() );
    LOG_INFO( "Steps read: " << walk.stepsCount() );
    LOG_INFO( "Prior param steps read: " << walk.priorParamsStepsCount() );
}

TEST( BIMAPWalk, DISABLED_serialize_read_compressed )
{
    boost::unordered_map<object_index_t, object_label_t> objects;
    boost::unordered_map<probe_index_t, probe_label_t> probes;
    ChessboardBiclusteringsIndexing    crossClusteringsIndexing;
    BIMAPWalk walk( crossClusteringsIndexing );
    std::ifstream compressed_ifs( "/home/astukalov/projects/pcp/results/egfr/egfr_cc_walk.xml.gz" );//"/home/astukalov/test2walk.xml" );//);
    boost::iostreams::filtering_stream<boost::iostreams::input> ifs;
    ifs.push( boost::iostreams::gzip_decompressor() );
    ifs.push( compressed_ifs );
    ASSERT_TRUE( ifs.good() );
    boost::archive::xml_iarchive ia(ifs);
    ia >> BOOST_SERIALIZATION_NVP( objects );
    LOG_INFO( "Objects count: " << objects.size() );
    ia >> BOOST_SERIALIZATION_NVP( probes );
    LOG_INFO( "Probes count: " << probes.size() );
    ia >> BOOST_SERIALIZATION_NVP( crossClusteringsIndexing );
    LOG_INFO( "Indexing size: " << crossClusteringsIndexing.size() );
    ia >> BOOST_SERIALIZATION_NVP( walk );
    EXPECT_NO_THROW( walk.check() );
    LOG_INFO( "Steps read: " << walk.stepsCount() );
    LOG_INFO( "Prior param steps read: " << walk.priorParamsStepsCount() );
}
