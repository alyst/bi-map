#include "TestCommon.h"

#include <boost/program_options.hpp>

namespace po = boost::program_options;

std::string BIMAPtest_data_path;

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    // Declare the supported options.
    po::options_description psibbtest_options_desc( "Additional options" );
    psibbtest_options_desc.add_options()
        ( "data_path", po::value< std::string >( &BIMAPtest_data_path )->default_value( "../../../test/data" ), "path to test data" );

    po::variables_map psibbtest_options;
    po::store( po::command_line_parser(argc, argv).options(psibbtest_options_desc).allow_unregistered().run(), psibbtest_options );
    po::notify( psibbtest_options );

    return RUN_ALL_TESTS();
}
