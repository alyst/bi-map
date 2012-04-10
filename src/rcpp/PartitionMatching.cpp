#include <cemm/bimap/BasicTypedefs.h>

#include <stdarg.h>
#include <cmath>

#include <boost/unordered_map.hpp>

#include <R_ext/Rdynload.h>
#include <cemm/RUtils.h>

#ifndef NDEBUG
//#define DYNLOAD_DEBUG
#endif

namespace cemm { namespace bimap {

struct PairsStats {
    size_t nCoCo;
    size_t nCoMm;
    size_t nMmCo;

    PairsStats()
    : nCoCo( 0 ), nCoMm( 0 ), nMmCo( 0 )
    {};

    PairsStats& operator+=( const PairsStats& that )
    {
        nCoCo += that.nCoCo;
        nCoMm += that.nCoMm;
        nMmCo += that.nMmCo;
        return ( *this );
    }
};

struct PartitionCollection {
    typedef std::string elm_label_t;
    typedef size_t elm_ix_t;
    typedef size_t clu_ix_t;
    typedef std::string ptn_ix_t;

    typedef std::set<elm_ix_t> elm_set_t;
    typedef boost::unordered_map<clu_ix_t, elm_set_t> clu2elmset_map_t;
    typedef boost::unordered_map<ptn_ix_t, clu2elmset_map_t> ptn_coll_t;
    typedef boost::unordered_map<elm_label_t, elm_ix_t> elm_label_map_t;
    typedef boost::unordered_map<ptn_ix_t, ptn_ix_t> ptn_map_t;

    elm_label_map_t     elmLabel2ix;    /// element label to element index map
    ptn_coll_t          ptnIx2clusters; /// partition index to its clusters map
    ptn_map_t           ptnIx2tmplPtnIx;/// partition index to template partition index
    std::size_t         _clustersCount; /// count of clusters in all partitions

    void read( const Rcpp::StringVector&    ptnId,
               const Rcpp::StringVector&    tmplPtnId,
               const Rcpp::IntegerVector&   cluId,
               const Rcpp::StringVector&    elmId,
               bool initElementIndexes = false );

    std::size_t elementsCount() const {
        return ( elmLabel2ix.size() );
    }

    std::size_t partitionsCount() const {
        return ( ptnIx2clusters.size() );
    }
    std::size_t clustersCount() const {
        return ( _clustersCount );
    }
    static PairsStats ClusterStats( const clu2elmset_map_t& partition, const elm_set_t& cluster );
};

void PartitionCollection::read(
   const Rcpp::StringVector&    ptnId,
   const Rcpp::StringVector&    tmplPtnId,
   const Rcpp::IntegerVector&   cluId,
   const Rcpp::StringVector&    elmId,
   bool                         initElementIndexes
){
    // reset maps
    if ( initElementIndexes ) elmLabel2ix.clear();
    ptnIx2clusters.clear();
    ptnIx2tmplPtnIx.clear();
    _clustersCount = 0;

    for ( int i = 0; i < ptnId.size(); i++ ) {
        // get the element
        elm_label_t elm = (elm_label_t)elmId[ i ];
        elm_label_map_t::const_iterator elmIt = elmLabel2ix.find( elm );
        elm_ix_t elmIx;
        if ( elmIt != elmLabel2ix.end() ) {
            elmIx = elmIt->second;
        } else {
            if ( initElementIndexes ) {
                // assign index to newly encountered element
                elmIx = elmLabel2ix.size();
                elmLabel2ix.insert( elmIt, elm_label_map_t::value_type( elm, elmIx ) );
            } else {
                THROW_RUNTIME_ERROR( "Element '" << elm << "' not found in template partition" );
            }
        }

        // get the partition
        ptn_ix_t ptnIx = (ptn_ix_t)ptnId[ i ];
        ptn_coll_t::iterator ptnIt = ptnIx2clusters.find( ptnIx );
        if ( ptnIt == ptnIx2clusters.end() ) {
            // create the new empty partition
            ptnIt = ptnIx2clusters.insert( ptnIt, ptn_coll_t::value_type( ptnIx, clu2elmset_map_t() ) );
        }
        // get the cluster index (local to the partition)
        clu_ix_t cluIx = cluId[ i ];
        // append cluster to the partition
        clu2elmset_map_t::iterator jt = ptnIt->second.find( cluIx );
        if ( jt != ptnIt->second.end() ) {
            // append element to the existing cluster
            jt->second.insert( elmIx );
        } else {
            // create new cluster (containing elmIx)
            elm_set_t eset;
            eset.insert( elmIx );
            ptnIt->second.insert( jt, clu2elmset_map_t::value_type( cluIx, eset ) );
            _clustersCount++;
        }
        // link partition and the template partition
        ptn_ix_t tmplPtnIx = (ptn_ix_t)tmplPtnId[ i ];
        ptn_map_t::const_iterator tmplPtnIt = ptnIx2tmplPtnIx.find( ptnIx );
        if ( tmplPtnIt == ptnIx2tmplPtnIx.end() ) {
            ptnIx2tmplPtnIx.insert( tmplPtnIt, ptn_map_t::value_type( ptnIx, tmplPtnIx ) );
        } else {
            if ( tmplPtnIt->second != tmplPtnIx ) {
                THROW_RUNTIME_ERROR( "Template partition Id '" << tmplPtnIx
                                     << "' is different from the one already stored: '" << tmplPtnIt->second << "'" );
            }
        }
    }
}

PairsStats PartitionCollection::ClusterStats(
    const clu2elmset_map_t& partition,
    const elm_set_t&        cluster
){
    const std::size_t clusterSize = cluster.size();
    PairsStats res;
    for ( clu2elmset_map_t::const_iterator partIt = partition.begin();
        partIt != partition.end(); ++partIt
    ){
        const elm_set_t& part = partIt->second;
        std::vector<elm_ix_t> isect;
        std::set_intersection( part.begin(), part.end(),
                               cluster.begin(), cluster.end(),
                               std::back_inserter( isect ) );
        const size_t isectSize = isect.size();
        res.nCoCo += isectSize * ( isectSize - 1 ); // matching pairs
        res.nCoMm += ( part.size() - isectSize ) * isectSize; // co-clustered in partition, but not in cluster
        res.nMmCo += ( clusterSize - isectSize ) * isectSize; // co-clustered in cluster, but mismatch in partition
    }
    return ( res );
}

/**
 *  Calculates number of matches and mismatches
 *  in co-clustered pairs between the "template" partition
 *  and a set of other partitions.
 */
RcppExport SEXP CountPartitionMismatches(
    SEXP    templatePartitionsCollectionExp,
    SEXP    partitionsCollectionExp,
    SEXP    templatePartitionIdColumnExp,
    SEXP    partitionIdColumnExp,
    SEXP    clusterIdColumnExp,
    SEXP    elementIdColumnExp
){
    BEGIN_RCPP

    std::string ptnIdColName = Rcpp::as<std::string>( partitionIdColumnExp );
    std::string cluIdColName = Rcpp::as<std::string>( clusterIdColumnExp );
    std::string elmIdColName = Rcpp::as<std::string>( elementIdColumnExp );
    std::string tmplPtnIdColName = Rcpp::as<std::string>( templatePartitionIdColumnExp );

    // read the template partition
    Rprintf( "Reading template partitions collection...\n" );
    PartitionCollection tmplPtnColl;
    {
        Rcpp::DataFrame ptnCollDf = templatePartitionsCollectionExp;
        tmplPtnColl.read( ptnCollDf[ tmplPtnIdColName ],
                          ptnCollDf[ tmplPtnIdColName ],
                          ptnCollDf[ cluIdColName ],
                          ptnCollDf[ elmIdColName ],
                          true );
    }
    size_t nElms = tmplPtnColl.elementsCount();

    Rprintf( "Reading partitions collection...\n" );
    PartitionCollection ptnColl;
    ptnColl.elmLabel2ix = tmplPtnColl.elmLabel2ix;
    {
        Rcpp::DataFrame ptnCollDf = partitionsCollectionExp;
        ptnColl.read( ptnCollDf[ ptnIdColName ],
                      ptnCollDf[ tmplPtnIdColName ],
                      ptnCollDf[ cluIdColName ],
                      ptnCollDf[ elmIdColName ],
                      false );
    }

    // calculate co-clustered and not co-clustered pairs
    Rprintf( "Calculating (%d partitions, %d templates)...\n",
             ptnColl.partitionsCount(), tmplPtnColl.partitionsCount() );
    Rcpp::DataFrame res;
    {
        Rcpp::StringVector ptnIdVec( ptnColl.partitionsCount() );
        Rcpp::StringVector tmplPtnIdVec( ptnColl.partitionsCount() );
        Rcpp::IntegerVector nCoCoVec( ptnColl.partitionsCount() );
        Rcpp::IntegerVector nCoMmVec( ptnColl.partitionsCount() );
        Rcpp::IntegerVector nMmCoVec( ptnColl.partitionsCount() );
        Rcpp::IntegerVector nMmMmVec( ptnColl.partitionsCount() );
        size_t i = 0;
        for ( PartitionCollection::ptn_coll_t::const_iterator ptnIt = ptnColl.ptnIx2clusters.begin();
              ptnIt != ptnColl.ptnIx2clusters.end(); ++ptnIt ) {
            PartitionCollection::ptn_map_t::const_iterator tmplIdIt = ptnColl.ptnIx2tmplPtnIx.find( ptnIt->first );
            if ( tmplIdIt == ptnColl.ptnIx2tmplPtnIx.end() ) {
                THROW_RUNTIME_ERROR( "Template ID not found for partition "
                                        << "'" << ptnIt->first << "'" );
            }
            PartitionCollection::ptn_coll_t::const_iterator tmplPtnIt = tmplPtnColl.ptnIx2clusters.find( tmplIdIt->second );
            if ( tmplPtnIt == tmplPtnColl.ptnIx2clusters.end() ) {
                THROW_RUNTIME_ERROR( "Template partition ID=" << tmplIdIt->second
                                        << " not found for partition "
                                        << "'" << tmplIdIt->first << "'" );
            }
            const PartitionCollection::clu2elmset_map_t& tmplPtn = tmplPtnIt->second;

            PairsStats ptnStats;
            size_t nPtnElems = 0;
            for ( PartitionCollection::clu2elmset_map_t::const_iterator tmplCluIt = tmplPtn.begin();
                  tmplCluIt != tmplPtn.end(); ++tmplCluIt
            ){
                nPtnElems += tmplCluIt->second.size();
                ptnStats += PartitionCollection::ClusterStats( ptnIt->second, tmplCluIt->second );
            }
            if ( nPtnElems > nElms ) {
                THROW_RUNTIME_ERROR( "Template partition '" << tmplPtnIt->first <<
                                     "' has " << nPtnElems << " elements, > " << nElms );
            }
            tmplPtnIdVec[ i ] = tmplIdIt->second;
            ptnIdVec[ i ] = ptnIt->first;
            nCoCoVec[ i ] = ptnStats.nCoCo / 2;
            nCoMmVec[ i ] = ptnStats.nCoMm / 2;
            nMmCoVec[ i ] = ptnStats.nMmCo / 2;
            nMmMmVec[ i ] = ( nPtnElems * ( nElms - 1 ) - ptnStats.nCoCo - ptnStats.nCoMm - ptnStats.nMmCo ) / 2;
            i++;
        }
        Rprintf( "Creating results data.frame\n" );
        res = Rcpp::DataFrame::create(
                Rcpp::Named( ptnIdColName, ptnIdVec ),
                Rcpp::Named( tmplPtnIdColName, tmplPtnIdVec ),
                Rcpp::Named( "co.tco", nCoCoVec ),
                Rcpp::Named( "co.tmismatch", nCoMmVec ),
                Rcpp::Named( "mismatch.tco", nMmCoVec ),
                Rcpp::Named( "mismatch.tmismatch", nMmMmVec ),
                Rcpp::Named( "stringsAsFactors", false )
             );
    }
    Rprintf( "Partitions Mismatches Calculation done...\n" );
    return ( res );
    END_RCPP
}

} }