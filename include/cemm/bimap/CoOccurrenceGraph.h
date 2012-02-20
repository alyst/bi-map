#pragma once

#include "BasicTypedefs.h"

#include <vector>
#include <set>
#include <boost/dynamic_bitset.hpp>

namespace cemm { namespace bimap {

typedef boost::dynamic_bitset<> observations_mask_t;
typedef size_t elm_index_t;
typedef std::set<elm_index_t> elm_adjacency_list_t;
typedef std::vector<elm_adjacency_list_t> co_occurrence_graph_t;

co_occurrence_graph_t CreateCoOccurrenceGraph( const std::vector<observations_mask_t>& occurrences, 
                                               const observations_mask_t&              observationsMask, 
                                               prob_t pValueThreshold = 0.05 );

} }