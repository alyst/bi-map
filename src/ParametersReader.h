#pragma once

#include "cemm/bimap/BasicTypedefs.h"
#include <string>

#include "cemm/bimap/BIMAPSampler.h"
#include "cemm/bimap/OPADataImportCSV.h"

namespace cemm { namespace bimap {

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

} }