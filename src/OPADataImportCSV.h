#pragma once

#include "BasicTypedefs.h"

#include "OPAData.h"

struct BIMAPIOParams {
    size_t          minCrossClusRefCount;
    size_t          minObjectsPtnRefCount;
    size_t          minProbesPtnRefCount;
    std::string     outputFilename;
    std::string     dataFilename;
    std::string     proteinsFilename;
    std::string     expDesignFilename;
    std::string     measurementsFilename;
    char            csvColumnSeparator;
    bool            mapBaitsToObjects;

    BIMAPIOParams();

    void checkFilenames();
};

OPAData OPADataImportCSV( BIMAPIOParams& ioParams );