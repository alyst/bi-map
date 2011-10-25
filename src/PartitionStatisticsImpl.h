#include "PartitionStatistics.h"

#include "IndexedPartitionsCollection.h"
#include "PartitionDataExtractor.h"

template<typename Part>
PartitionPartsPDF::PartitionPartsPDF(
    const IndexedPartitionsCollection<Part>& ptnColl
){
    typedef IndexedPartitionsCollection<Part> ptn_coll_type;
    typedef PartitionDataExtractor<Part> extractor_type;
    typedef typename extractor_type::partition_indexing ptn_indexing_type;
    typedef typename extractor_type::const_partition_iterator ptn_iterator;
    typedef typename extractor_type::const_part_iterator part_iterator;

    for ( ptn_iterator ptnit = extractor_type::PartitionsBegin( ptnColl.walk() );
          ptnit != extractor_type::PartitionsEnd( ptnColl.walk() ); ++ptnit
    ){
        log_prob_t  lnPartCountsProd = 0;
        size_t partsCnt = 0;
        for ( part_iterator pit = (*ptnit)->value().begin();
            pit != (*ptnit)->value().end(); ++pit
        ){
            lnPartCountsProd += log( ptnColl.partStats( (*pit)->serial() ).nsteps );
            partsCnt++;
        }
        _ptnPartsFreq[ (*ptnit)->serial() ] = lnPartCountsProd
                                      - partsCnt * log( ptnColl.walk().stepsCount() );
    }
}
