#pragma once

#include <OPAData.h>

OPAData generateTestOPAData();
OPAData loadTestOPAData( const char* filename );
void saveTestOPAData( const OPAData& data, const char* filename );
