#pragma once

#include "BasicTypedefs.h"

#include "OPAData.h"

OPAData OPADataImportCSV( const char* proteinsFilename,
                          const char* expDesignFilename,
                          const char* measurementsFilename,
                          bool        mapBaitsToObjects = true,
                          char        sep = '\t' );