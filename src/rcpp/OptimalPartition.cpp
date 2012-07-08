#include <cemm/bimap/BasicTypedefs.h>

#include <stdarg.h>
#include <cmath>

#include <boost/unordered_map.hpp>

#include <R_ext/Rdynload.h>
#include <cemm/RUtils.h>

#include <cemm/containers/PartitionOptimizer.h>

#ifndef NDEBUG
//#define DYNLOAD_DEBUG
#endif

namespace cemm { namespace bimap {

struct ClustersCollection {
    typedef std::string elm_label_t;
    typedef size_t elm_ix_t;
    typedef size_t clu_ix_t;
    typedef std::string ptn_ix_t;

    typedef double score_t;

    typedef std::set<elm_ix_t> elm_set_t;
    typedef cemm::containers::partition::PartitionOptimizer<elm_set_t> ptn_optimizer_t;

    typedef ptn_optimizer_t::part_map_t clu2elmset_map_t;
    typedef boost::unordered_map<elm_label_t, elm_ix_t> elm_label_map_t;
    typedef boost::unordered_map<clu_ix_t, score_t> clu2score_map_t;


    elm_label_map_t     elmLabel2ix;    /// element label to element index map
    clu2elmset_map_t    clusters;       /// cluster elements map
    clu2score_map_t     clustersScore;  /// cluster props map
    std::vector<elm_label_t>    elmIx2LabelMap;

    void read_clusters( const Rcpp::IntegerVector&   cluId,
                        const Rcpp::StringVector&    elmId );
    void read_score( const Rcpp::IntegerVector&    cluId,
                     const Rcpp::NumericVector&    score );

    std::size_t elementsCount() const {
        return ( elmLabel2ix.size() );
    }
};

void ClustersCollection::read_clusters(
   const Rcpp::IntegerVector&   cluId,
   const Rcpp::StringVector&    elmId
){
    LOG_INFO( "Reading clusters..." );
    for ( int i = 0; i < cluId.size(); i++ ) {
        // get the element
        elm_label_t elm = (elm_label_t)elmId[ i ];
        elm_label_map_t::const_iterator elmIt = elmLabel2ix.find( elm );
        elm_ix_t elmIx;
        if ( elmIt != elmLabel2ix.end() ) {
            elmIx = elmIt->second;
        } else {
            // assign index to newly encountered element
            elmIx = elmLabel2ix.size();
            elmLabel2ix.insert( elmIt, elm_label_map_t::value_type( elm, elmIx ) );
            elmIx2LabelMap.push_back( elm );
        }

        // get the cluster index (local to the partition)
        clu_ix_t cluIx = cluId[ i ];
        // append cluster to the partition
        clu2elmset_map_t::iterator jt = clusters.find( cluIx );
        if ( jt != clusters.end() ) {
            // append element to the existing cluster
            jt->second.insert( elmIx );
        } else {
            // create new cluster (containing elmIx)
            elm_set_t eset;
            eset.insert( elmIx );
            clusters.insert( jt, clu2elmset_map_t::value_type( cluIx, eset ) );
        }
    }
    LOG_INFO( clusters.size() << " cluster(s) read" );
}

void ClustersCollection::read_score(
   const Rcpp::IntegerVector&   cluId,
   const Rcpp::NumericVector&   score
){
    LOG_INFO( "Reading cluster properties..." );
    for ( int i = 0; i < cluId.size(); i++ ) {
        // get the cluster index (local to the partition)
        clu_ix_t cluIx = cluId[ i ];
        // append cluster to the partition
        clu2score_map_t::const_iterator jt = clustersScore.find( cluIx );
        if ( jt != clustersScore.end() ) {
            THROW_RUNTIME_ERROR( "Cluster #" << cluIx << " already has properties defined" );
        } else {
            clustersScore.insert( jt, std::make_pair( cluIx, score[i] ) );
        }
    }
    LOG_INFO( clustersScore.size() << " cluster(s) read" );
}

struct PartitionScore {
    typedef double score_type;
    typedef int aux_data_type;
    typedef cemm::containers::partition::part_id_set_t clu_id_set_t;

    const ClustersCollection&   cluColl;

    PartitionScore( const ClustersCollection& cluColl )
    : cluColl( cluColl )
    {}

    aux_data_type prepareAuxData( const ClustersCollection::elm_set_t& componentElements ) const
    {
        return ( 0 );
    }
    score_type operator()( const ClustersCollection::clu2elmset_map_t& allSubcomplexes,
                           const int& unused,
                           const clu_id_set_t& subPtn ) const
    {
        score_type  res = 0.0;
        for ( clu_id_set_t::const_iterator it = subPtn.begin(); it != subPtn.end(); ++it ) {
            ClustersCollection::clu2score_map_t::const_iterator jt = cluColl.clustersScore.find( *it );
            if ( jt == cluColl.clustersScore.end() ) {
                THROW_RUNTIME_ERROR( "Cluster #" << *it << " unknown" );
            }
            res += jt->second;
        }
        return ( res );
    }
};

/**
 *  Calculates number of matches and mismatches
 *  in co-clustered pairs between the "template" partition
 *  and a set of other partitions.
 */
RcppExport SEXP OptimalPartition(
    SEXP    clustersContentsDataFrameExp,
    SEXP    clustersPropsDataFrameExp,
    SEXP    clusterIdColumnExp,
    SEXP    elementIdColumnExp,
    SEXP    scoreColumnExp
){
    BEGIN_RCPP

    std::string cluIdColName = Rcpp::as<std::string>( clusterIdColumnExp );
    std::string elmIdColName = Rcpp::as<std::string>( elementIdColumnExp );
    std::string scoreColName = Rcpp::as<std::string>( scoreColumnExp );

    // read the template partition
    ClustersCollection clustersColl;
    {
        Rcpp::DataFrame clustersContentsDf = clustersContentsDataFrameExp;
        clustersColl.read_clusters( clustersContentsDf[ cluIdColName ],
                                    clustersContentsDf[ elmIdColName ] );
    }
    {
        Rcpp::DataFrame clustersPropsDf = clustersPropsDataFrameExp;
        clustersColl.read_score( clustersPropsDf[ cluIdColName ],
                                 clustersPropsDf[ scoreColName ] );
    }

    LOG_INFO( "Preprocessing clusters collection..." );
    ClustersCollection::ptn_optimizer_t ptnOptimizer( clustersColl.clusters.begin(), clustersColl.clusters.end() );

    LOG_INFO( "Calculating optimal partition..." );
    ClustersCollection::ptn_optimizer_t::new_partition_t bestPtn = ptnOptimizer.bestMatch( PartitionScore( clustersColl ) );
    LOG_INFO( "Optimal partition calculation done" );

    LOG_INFO( "Exporting partition" );
    {
        // new parts are ignored
        std::size_t nElems = 0;
        for ( std::size_t cluIx = 0; cluIx < bestPtn.first.size(); ++cluIx ) {
            nElems += clustersColl.clusters.find( bestPtn.first[cluIx] )->second.size();
        }
        Rcpp::StringVector  elmIdVec( nElems );
        Rcpp::IntegerVector cluIdVec( elmIdVec.size() );
        std::size_t i = 0;
        for ( std::size_t cluIx = 0; cluIx < bestPtn.first.size(); ++cluIx ) {
            const ClustersCollection::clu_ix_t cluId = bestPtn.first[cluIx];
            const ClustersCollection::elm_set_t& elmSet = clustersColl.clusters.find( cluId )->second;
            for ( ClustersCollection::elm_set_t::const_iterator elmIt = elmSet.begin();
                  elmIt != elmSet.end(); ++elmIt
            ){
                elmIdVec[i] = clustersColl.elmIx2LabelMap[ *elmIt ];
                cluIdVec[i] = cluId;
                i++;
            }
        }
        LOG_INFO( "Partition with " << nElems << " element(s) in " 
                  << bestPtn.first.size() << " cluster(s) exported" );
        return ( Rcpp::List::create(
            Rcpp::Named( elmIdColName ) = elmIdVec,
            Rcpp::Named( cluIdColName ) = cluIdVec,
            Rcpp::Named( "stringsAsFactors" ) = false
        ) );
    }

    END_RCPP
}

} }