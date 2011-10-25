#include <set>
#include <fstream>

#include "TestCommon.h"

#include <gtest/gtest.h>
#include <boost/archive/xml_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/filesystem.hpp>

#include <OPAData.h>
#include "OPADataTest.h"

OPAData generateTestOPAData()
{
    OPAData    res;
    res.addObject( "A", 2 );
    res.addObject( "B", 3 );
    res.addObject( "C", 1 );
    res.addObject( "D", 2 );

    res.addProbe( "X", "A" );
    res.addProbe( "Y", "B" );
    res.addProbe( "Z", "C" );

    res.addAssay( "X1", "X" );
    res.addAssay( "Y1", "Y" );
    res.addAssay( "Y2", "Y" );
    res.addAssay( "Z1", "Z" );
    res.addAssay( "Z2", "Z" );

    res.addMeasurement( "A", "X1", 100 );
    res.addMeasurement( "A", "Z1", 50 );
    res.addMeasurement( "A", "Z2", 40 );
    res.addMeasurement( "B", "Y1", 20 );
    res.addMeasurement( "B", "Y2", 18 );
    res.addMeasurement( "C", "X1", 80 );
    res.addMeasurement( "C", "Z1", 18 );
    res.addMeasurement( "C", "Z2", 25 );
    res.addMeasurement( "D", "Y1", 10 );
    res.addMeasurement( "D", "Y2", 11 );

    return ( res );
}

OPAData loadTestOPAData( const char* filename )
{
    boost::filesystem::path test_data_path( BIMAPtest_data_path );
    if ( !boost::filesystem::exists( test_data_path ) || !boost::filesystem::is_directory( test_data_path ) ) {
        THROW_RUNTIME_ERROR( "OPAData test data folder not found: " << test_data_path );
    }
    boost::filesystem::path test_data_file = test_data_path / filename;
    if ( !boost::filesystem::exists( test_data_file ) ) {
        THROW_RUNTIME_ERROR( "OPAData test data file not found: " << filename );
    }
    std::ifstream dataStream( test_data_file.string().c_str(), std::ios_base::in );
    boost::archive::xml_iarchive dataArchive( dataStream );
    OPAData data;
    dataArchive >> BOOST_SERIALIZATION_NVP( data );
    return ( data );
}

void saveTestOPAData( const OPAData& data, const char* filename )
{
    boost::filesystem::path test_data_path( BIMAPtest_data_path );
    if ( !boost::filesystem::exists( test_data_path ) || !boost::filesystem::is_directory( test_data_path ) ) {
        THROW_RUNTIME_ERROR( "OPAData test data folder not found: " << test_data_path );
    }
    boost::filesystem::path test_data_file = test_data_path / filename;
    if ( boost::filesystem::exists( test_data_file ) ) {
        LOG_INFO( "Overwriting " << filename );
    }
    std::ofstream dataStream( test_data_file.string().c_str(), std::ios_base::out );
    boost::archive::xml_oarchive dataArchive( dataStream );
    dataArchive << BOOST_SERIALIZATION_NVP( data );
}

/*************************************************/

TEST( OPAData, DISABLED_edit )
{
    OPAData    res;
    EXPECT_EQ( res.objectsCount(), 0 );
    EXPECT_EQ( res.probesCount(), 0 );
    EXPECT_EQ( res.assaysCount(), 0 );

    OPAData::const_object_ptr_t pA;
    ASSERT_TRUE( ( pA = res.addObject( "A", 2 ) ) );
    EXPECT_EQ( ((const OPAData&)res).object( "A" ), pA );
    EXPECT_EQ( pA->index(), 0 );
    EXPECT_EQ( pA->label(), "A" );
    EXPECT_EQ( pA->sequenceLength(), 2 );
    EXPECT_EQ( res.objectsCount(), 1 );
    EXPECT_THROW( res.addObject( "A", 3 ), entity_already_exists );

    OPAData::const_probe_ptr_t pX;
    ASSERT_TRUE( ( pX = res.addProbe( "X", "A" ) ) );
    EXPECT_EQ( ((const OPAData&)res).probe( "X" ), pX );
    EXPECT_EQ( pX->index(), 0 );
    EXPECT_EQ( pX->baitIndex(), pA->index() );
    EXPECT_EQ( pX->label(), "X" );
    EXPECT_EQ( res.probesCount(), 1 );
    EXPECT_EQ( pX->assayIndexes().size(), 0 );

    EXPECT_THROW( res.addProbe( "X", "C" ), entity_not_found );
    EXPECT_THROW( res.addProbe( "X", "A" ), entity_already_exists );

    assay_index_t pX1 = res.addAssay( "X1", "X" );
    EXPECT_TRUE( ((const OPAData&)res).assayIndex( "X1" ) == pX1 );
    EXPECT_EQ( pX1, 0 );
    EXPECT_EQ( res.assay( pX1 ).label(), "X1" );
    EXPECT_EQ( pX->assayIndexes().size(), 1 );
    EXPECT_EQ( res.assaysCount(), 1 );

    OPAData::const_object_ptr_t pB;
    ASSERT_TRUE( ( pB = res.addObject( "B", 2 ) ) );
    EXPECT_EQ( res.objectsCount(), 2 );

    OPAData::const_probe_ptr_t pY;
    ASSERT_TRUE( ( pY = res.addProbe( "Y", "B" ) ) );
    EXPECT_EQ( res.probesCount(), 2 );

    assay_index_t pX2 = res.addAssay( "X2", "X" );
    assay_index_t pY1 = res.addAssay( "Y1", "Y" );
    EXPECT_EQ( res.assaysCount(), 3 );

    EXPECT_NO_THROW( res.addMeasurement( "A", "X1", 5 ) );
    EXPECT_THROW( res.addMeasurement( "C", "X1", 5 ), entity_not_found );
    EXPECT_THROW( res.addMeasurement( "A", "X3", 5 ), entity_not_found );
    EXPECT_EQ( res.measurement( pA->index(), pX1 ), OPAData::celldata_t( 5u ) );
    EXPECT_EQ( res.measurement( pA->index(), pX2 ), OPAData::celldata_t( 0u ) );
    EXPECT_EQ( res.measurement( 0, 0 ), OPAData::celldata_t( 5u ) );
    EXPECT_NO_THROW( res.addMeasurement( "B", "X1", 2 ) );
    EXPECT_NO_THROW( res.addMeasurement( "B", "Y1", 4 ) );
}

TEST( OPAData, DISABLED_save_BIMAP_tests_input )
{
    {
//      msrun
//protein  1-1  1-2  1-3  1-4  1-5  1-6
//      A 1003 1046 1029 1022 1045 1068
//      probe signal = 5, object.multiple = 1
        OPAData    data;
        data.addObject( "A", 50 );

        data.addProbe( "1", "A" );

        data.addAssay( "1-1", "1" );
        data.addAssay( "1-2", "1" );
        data.addAssay( "1-3", "1" );
        data.addAssay( "1-4", "1" );
        data.addAssay( "1-5", "1" );
        data.addAssay( "1-6", "1" );

        data.addMeasurement( "A", "1-1", 1003 );
        data.addMeasurement( "A", "1-2", 1046 );
        data.addMeasurement( "A", "1-3", 1029 );
        data.addMeasurement( "A", "1-4", 1022 );
        data.addMeasurement( "A", "1-5", 1045 );
        data.addMeasurement( "A", "1-5", 1068 );
        saveTestOPAData( data, "osadata_1x1.xml" );
    }
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

        OPAData    data;
        data.addObject( "A", 50 );
        data.addObject( "B", 50 );
        data.addObject( "C", 50 );
        
        data.addProbe( "1", "A" );
        
        data.addAssay( "1-1", "1" );
        data.addAssay( "1-2", "1" );
        data.addAssay( "1-3", "1" );
        data.addAssay( "1-4", "1" );
        
        data.addMeasurement( "A", "1-1", 59 );
        data.addMeasurement( "B", "1-1", 57 );
        data.addMeasurement( "C", "1-1", 20938 );
        data.addMeasurement( "A", "1-2", 51 );
        data.addMeasurement( "B", "1-2", 63 );
        data.addMeasurement( "C", "1-2", 21097 );
        data.addMeasurement( "A", "1-3", 45 );
        data.addMeasurement( "B", "1-3", 55 );
        data.addMeasurement( "C", "1-3", 21095 );
        data.addMeasurement( "A", "1-4", 51 );
        data.addMeasurement( "B", "1-4", 58 );
        data.addMeasurement( "C", "1-4", 20997 );
        saveTestOPAData( data, "osadata_3x1_2_clusters.xml" );
    }
    {
        OPAData    data;
        data.addObject( "A", 50 );
        data.addObject( "B", 50 );
        data.addObject( "C", 50 );
        data.addObject( "D", 50 );

        data.addProbe( "1", "A" );
        data.addProbe( "2", "B" );

        data.addAssay( "1-1", "1" );
        data.addAssay( "1-2", "1" );
        data.addAssay( "1-3", "1" );
        data.addAssay( "1-4", "1" );
        data.addAssay( "2-1", "2" );
        data.addAssay( "2-2", "2" );
        data.addAssay( "2-3", "2" );
        data.addAssay( "2-4", "2" );

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

        saveTestOPAData( data, "osadata_4x2_signal_noise.xml" );
    }
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

        saveTestOPAData( data, "osadata_4x4_diag.xml" );
    }
    {
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

//        msrun
// protein 1-1 1-2 1-3 1-4 2-1 2-2 2-3 2-4
//       A  58  45  34  46  10   7   8  13
//       B  47  65  41  47  16   9  14   5
//       C  21  16  14  17  23  20  43  29
//       D  23  19  17  24  24  32  26  28
        data.addMeasurement( "A", "1-1", 58 );
        data.addMeasurement( "B", "1-1", 47 );
        data.addMeasurement( "C", "1-1", 21 );
        data.addMeasurement( "D", "1-1", 23 );
        data.addMeasurement( "A", "1-2", 45 );
        data.addMeasurement( "B", "1-2", 65 );
        data.addMeasurement( "C", "1-2", 16 );
        data.addMeasurement( "D", "1-2", 19 );
        data.addMeasurement( "A", "1-3", 34 );
        data.addMeasurement( "B", "1-3", 41 );
        data.addMeasurement( "C", "1-3", 14 );
        data.addMeasurement( "D", "1-3", 17 );
        data.addMeasurement( "A", "1-4", 46 );
        data.addMeasurement( "B", "1-4", 47 );
        data.addMeasurement( "C", "1-4", 17 );
        data.addMeasurement( "D", "1-4", 24 );

        data.addMeasurement( "A", "2-1", 10 );
        data.addMeasurement( "B", "2-1", 16 );
        data.addMeasurement( "C", "2-1", 23 );
        data.addMeasurement( "D", "2-1", 24 );
        data.addMeasurement( "A", "2-2", 7 );
        data.addMeasurement( "B", "2-2", 9 );
        data.addMeasurement( "C", "2-2", 20 );
        data.addMeasurement( "D", "2-2", 32 );
        data.addMeasurement( "A", "2-3", 8 );
        data.addMeasurement( "B", "2-3", 14 );
        data.addMeasurement( "C", "2-3", 43 );
        data.addMeasurement( "D", "2-3", 26 );
        data.addMeasurement( "A", "2-4", 13 );
        data.addMeasurement( "B", "2-4", 5 );
        data.addMeasurement( "C", "2-4", 29 );
        data.addMeasurement( "D", "2-4", 28 );

//        msrun
// protein 3-1 3-2 3-3 3-4 4-1 4-2 4-3 4-4
//       A   0   0   0   0   0   0   0   0
//       B   0   0   0   0   0   0   1   0
//       C  45  48  59  53  37  35  37  36
//       D  72  50  65  54  35  39  34  28
        data.addMeasurement( "A", "3-1", 0 );
        data.addMeasurement( "B", "3-1", 0 );
        data.addMeasurement( "C", "3-1", 45 );
        data.addMeasurement( "D", "3-1", 72 );
        data.addMeasurement( "A", "3-2", 0 );
        data.addMeasurement( "B", "3-2", 0 );
        data.addMeasurement( "C", "3-2", 48 );
        data.addMeasurement( "D", "3-2", 50 );
        data.addMeasurement( "A", "3-3", 0 );
        data.addMeasurement( "B", "3-3", 0 );
        data.addMeasurement( "C", "3-3", 59 );
        data.addMeasurement( "D", "3-3", 65 );
        data.addMeasurement( "A", "3-4", 0 );
        data.addMeasurement( "B", "3-4", 0 );
        data.addMeasurement( "C", "3-4", 53 );
        data.addMeasurement( "D", "3-4", 54 );

        data.addMeasurement( "A", "4-1", 0 );
        data.addMeasurement( "B", "4-1", 0 );
        data.addMeasurement( "C", "4-1", 37 );
        data.addMeasurement( "D", "4-1", 35 );
        data.addMeasurement( "A", "4-2", 0 );
        data.addMeasurement( "B", "4-2", 0 );
        data.addMeasurement( "C", "4-2", 35 );
        data.addMeasurement( "D", "4-2", 39 );
        data.addMeasurement( "A", "4-3", 0 );
        data.addMeasurement( "B", "4-3", 0 );
        data.addMeasurement( "C", "4-3", 37 );
        data.addMeasurement( "D", "4-3", 34 );
        data.addMeasurement( "A", "4-4", 0 );
        data.addMeasurement( "B", "4-4", 0 );
        data.addMeasurement( "C", "4-4", 36 );
        data.addMeasurement( "D", "4-4", 28 );
        saveTestOPAData( data, "osadata_4x4_triag_lower.xml" );
    }
}
