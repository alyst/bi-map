#pragma once

#include "BasicTypedefs.h"

#include "BIMAPWalk.h"
#include "ChessboardBiclusteringsPDFEval.h"

void BIMAPResultsLoad(
    const char* filename, ChessboardBiclusteringsIndexing& indexing,
    objects_label_map_type& objects,
    probes_label_map_type&  probes,
    BIMAPWalk** walk, ChessboardBiclusteringsPDFEval** pdfAdjust );

void BIMAPResultsSave( const char* filename, 
           const objects_label_map_type& objects, const probes_label_map_type&  probes,
           const BIMAPWalk& walk, const ChessboardBiclusteringsPDFEval* pdfAdjust );
