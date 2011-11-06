#pragma once

#include "BasicTypedefs.h"
#include <string>

#include "BIMAPSampler.h"
#include "OPADataImportCSV.h"

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

