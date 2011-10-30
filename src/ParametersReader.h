#pragma once

#include "BasicTypedefs.h"
#include <string>

#include "BIMAPSampler.h"

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

    BIMAPIOParams();
};

bool BIMAPParamsRead(
    int argc, char* argv[],
    ChessboardBiclusteringHyperPriors&     hyperpriors,
    GibbsSamplerParams&             gibbsParams,
    TurbineCascadeParams&           cascadeParams,
    CellSignalParams&               signalParams,
    PrecomputedDataParams&          precomputedDataParams,
    ChessboardBiclusteringPriors&   priors,
    BIMAPSampleCollectorParams&     collectorParams,
    BIMAPIOParams&                  ioParams
);

