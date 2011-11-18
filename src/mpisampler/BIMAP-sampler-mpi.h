#pragma once

#include "../BasicTypedefs.h"

#include "../BIMAPSampler.h"
#include "MPIUnitCommunicator.h"

BOOST_IS_MPI_DATATYPE( ChessboardBiclusteringEnergyEval )

boost::optional<BIMAPWalk> MPI_BIMAPSampler_run(
    const BIMAPSamplerHelper&      helper,
    boost::mpi::communicator&       communicator,
    ChessboardBiclusteringsIndexing*       ccIndexing,
    const TurbineCascadeParams&     eeCascadeParams,
    const ChessboardBiclustering&          iniClus,
    const TurbineCascadeExecutionMonitor*   pMon = NULL,
    size_t      walkSamples = 100,
    size_t      priorsStoragePeriod = 5,
    size_t      samplesReportingPeriod = 100,
    double      maxSamplesReportingDelay = 60
);
