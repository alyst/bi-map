#pragma once

#include <cemm/bimap/OPAData.h>

namespace cemm { namespace test {

using namespace cemm::bimap;

OPAData generateTestOPAData();
OPAData loadTestOPAData( const char* filename );
void saveTestOPAData( const OPAData& data, const char* filename );

} }
