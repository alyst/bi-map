#include <boost/format.hpp>

#include <cemm/math/Distributions.h>

#include "cemm/bimap/CoOccurrenceGraph.h"

namespace cemm { namespace bimap {

/**
 * Calculates graph of elements co-occurrence based on
 * independent occurrence in observations
 * (hypergeometric distribution is used),
 * 
 * @return co-occurrence undirected graph as adjacency list
 */
co_occurrence_graph_t CreateCoOccurrenceGraph(
    const std::vector<observations_mask_t>& occurrences,        /** binary matrix of elements seen in observations */
    const observations_mask_t&              observationsMask,   /** mask of observations to use */
    prob_t                                  pValueThreshold     /** p-value of elements independence hypethesis,
                                                                    used as cutoff for co-occurence graph */
){
    size_t      nObservations = observationsMask.size();
    size_t      nElems = occurrences.size();
    log_prob_t  lpThreshold = log( pValueThreshold );
    log_prob_t  lpMean = 0;
    log_prob_t  lpVar = 0;
    size_t      nEdges = 0;
    size_t      nUsedObservations = observationsMask.count();

    co_occurrence_graph_t res( nElems );
    for ( elm_index_t i = 0; i < nElems; i++ ) {
        for ( elm_index_t j = i; j < nElems; j++ ) {
            size_t nCoOccurs = ( occurrences[i] & occurrences[j] & observationsMask ).count();
            size_t nIOccurs = ( occurrences[i] & observationsMask ).count();
            size_t nJOccurs = ( occurrences[j] & observationsMask ).count();
            // p-value that there could be a better overlap
            HypergeometricDistribution coOccurrence( nIOccurs, nUsedObservations - nIOccurs, nJOccurs );
            log_prob_t lpValue = coOccurrence.lnCdf_Q( nCoOccurs, true );
            if ( is_unset( lpValue ) ) lpValue = 0;
            lpMean += lpValue;
            lpVar += lpValue * lpValue;

            // add i-j, j-i link
            // if p-value is good, or too little data
            // or both elements are frequent
            if ( ( lpValue <= lpThreshold )
              || ( nObservations <= 4 )
              || ( nCoOccurs >= 0.75 * nObservations )
            ) {
                if ( i != j ) nEdges++;
                res[ i ].insert( j );
                res[ j ].insert( i );
            }
        }
    }
    lpMean /= nElems * ( nElems + 1 ) / 2;
    lpVar /= nElems * ( nElems + 1 ) / 2;
    lpVar -= lpMean * lpMean;
    lpVar = sqrt( lpVar );

    LOG_INFO( "n_elem=" << nElems << " n_obs=" << nObservations 
              << "(" << nUsedObservations << " used)"
              << " n_edges=" << nEdges
              << " ln_threshold=" << boost::format("%.3f") % lpThreshold
              << " E[ln(p)]=" << boost::format("%.3f") % lpMean 
              << " D[ln(p)]=" << boost::format("%.3f") % lpVar );
    return ( res );
}

} }