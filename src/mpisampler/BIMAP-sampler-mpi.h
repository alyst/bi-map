#pragma once

#include "cemm/bimap/BasicTypedefs.h"

#include <cemm/eesampler/MPIUnitCommunicator.h>

#include "cemm/bimap/BIMAPSampler.h"

BOOST_IS_MPI_DATATYPE( cemm::bimap::LLHPartitionWeights )
BOOST_IS_MPI_DATATYPE( cemm::bimap::LLHWeights )
BOOST_IS_MPI_DATATYPE( cemm::bimap::ChessboardBiclusteringEnergyEval )

namespace cemm { namespace bimap {

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

} }