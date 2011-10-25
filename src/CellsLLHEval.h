/******************************************************************************
    This file is a part of BI-MAP project.

    Copyright (C) 2009-2011 Alexey Stukalov

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any
    later version, which shall act as a proxy defined in
    Section 6 of version 3 of the license.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public 
    License along with this library.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/

#pragma once

#include "BasicTypedefs.h"

#include "PrecomputedData.h"

#include "math/Distributions.h"
#include "ChessboardBiclustering.h"

/**
    Likelihood evaluations for OPAData.
    Evaluates likelihoods of types for cell of objects x probes matrix.
 */
class CellsLLHEval {
public:
    typedef ChessboardBiclusteringData::noise_params_type noise_params_type;
    typedef ChessboardBiclusteringData::signal_params_type signal_params_type;
    typedef array2d<log_prob_t> lnprob_matrix_type;

protected:
    const PrecomputedData&      _precomputed;
    const ChessboardBiclusteringData&  _clusData;
    noise_params_type           _noiseModel;

public:

    CellsLLHEval( const PrecomputedData&        precomputed,
                  const ChessboardBiclusteringData&    clusData );

    const OPAData& data() const {
        return ( _precomputed.data() );
    }

    const ChessboardBiclusteringData& crossClusteringParams() const {
        return ( _clusData );
    }

    const CellSignalParams& measurementParams() const {
        return ( _precomputed.signalParams() );
    }

    log_prob_t cellNoiseLnPdf( object_index_t objectIx, probe_index_t probeIx ) const;
    log_prob_t cellSignalLnPdf( object_index_t objectIx, probe_index_t probeIx ) const;

    static log_prob_t CrossClusterBlockLLH(
        const object_set_t&         objects,
        const probe_bitset_t&       probes,
        const lnprob_matrix_type&   lnPdf
    );
};
