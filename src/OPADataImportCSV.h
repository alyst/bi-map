#pragma once

#include "BasicTypedefs.h"

#include "OPAData.h"

OPAData OPADataImportCSV( const char* proteinsFilename,
                          const char* expDesignFilename,
                          const char* measurementsFilename,
                          char        sep = '\t' );