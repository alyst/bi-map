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

    typedef std::string elm_label_t;
    typedef size_t elm_ix_t;
    typedef size_t clu_ix_t;
    typedef std::string ptn_ix_t;

    typedef std::set<elm_ix_t> elm_set_t;
    typedef boost::unordered_map<clu_ix_t, elm_set_t> clu2elmset_map_t;
    typedef boost::unordered_map<ptn_ix_t, clu2elmset_map_t> ptn_coll_t;
    typedef boost::unordered_map<elm_label_t, elm_ix_t> elm_label_map_t;
    typedef boost::unordered_map<ptn_ix_t, ptn_ix_t> ptn_map_t;

    size_t nElms;

    clu2elmset_map_t    tmplClu2ElmSet;
    elm_label_map_t     elmLabelMap;
    // read the template partition
    Rprintf( "Reading template partitions collection...\n" );
    ptn_coll_t tmplPtnColl;
    // read the partitions of the collection while simultaneously calculating the scores
    {
        Rcpp::DataFrame ptnCollDf = templatePartitionsCollectionExp;
        Rcpp::StringVector ptnCollPtnId = ptnCollDf[ tmplPtnIdColName ];
        Rcpp::IntegerVector ptnCollCluId = ptnCollDf[ cluIdColName ];
        Rcpp::StringVector ptnCollElmId = ptnCollDf[ elmIdColName ];
        for ( int i = 0; i < ptnCollPtnId.size(); i++ ) {
            ptn_ix_t ptnIx = (ptn_ix_t)ptnCollPtnId[ i ];
            clu_ix_t cluIx = ptnCollCluId[ i ];
            ptn_coll_t::iterator ptnIt = tmplPtnColl.find( ptnIx );
            if ( ptnIt == tmplPtnColl.end() ) {
                ptnIt = tmplPtnColl.insert( ptnIt, ptn_coll_t::value_type( ptnIx, clu2elmset_map_t() ) );
            }
            elm_label_t elm = (elm_label_t)ptnCollElmId[ i ];
            elm_label_map_t::const_iterator elmIt = elmLabelMap.find( elm );
            elm_ix_t elmIx;
            if ( elmIt != elmLabelMap.end() ) {
                elmIx = elmIt->second;
            } else {
                elmIx = elmLabelMap.size();
                elmLabelMap.insert( elmIt, elm_label_map_t::value_type( elm, elmIx ) );
            }
            clu2elmset_map_t::iterator jt = ptnIt->second.find( cluIx );
            if ( jt != ptnIt->second.end() ) {
                jt->second.insert( elmIx );
            } else {
                elm_set_t eset;
                eset.insert( elmIx );
                ptnIt->second.insert( jt, clu2elmset_map_t::value_type( cluIx, eset ) );
            }
        }
        nElms = elmLabelMap.size();
    }

    Rprintf( "Reading partitions collection...\n" );
    ptn_coll_t ptnColl;
    ptn_map_t  ptnToTmplPtnMap;
    // read the partitions of the collection while simultaneously calculating the scores
    {
        Rcpp::DataFrame ptnCollDf = partitionsCollectionExp;
        Rcpp::StringVector ptnCollPtnId = ptnCollDf[ ptnIdColName ];
        Rcpp::IntegerVector ptnCollCluId = ptnCollDf[ cluIdColName ];
        Rcpp::StringVector ptnCollElmId = ptnCollDf[ elmIdColName ];
        Rcpp::StringVector tmplPtnId = ptnCollDf[ tmplPtnIdColName ];

        for ( int i = 0; i < ptnCollPtnId.size(); i++ ) {
            ptn_ix_t ptnIx = (ptn_ix_t)ptnCollPtnId[ i ];
            clu_ix_t cluIx = ptnCollCluId[ i ];
            elm_label_t elm = (elm_label_t)ptnCollElmId[ i ];
            elm_label_map_t::const_iterator elmIt = elmLabelMap.find( elm );
            if ( elmIt == elmLabelMap.end() ) THROW_RUNTIME_ERROR( "Element '" << elm << "' not found in template partition" );
            elm_ix_t elmIx = elmIt->second;
            ptn_coll_t::iterator ptnIt = ptnColl.find( ptnIx );
            if ( ptnIt == ptnColl.end() ) {
                ptnIt = ptnColl.insert( ptnIt, ptn_coll_t::value_type( ptnIx, clu2elmset_map_t() ) );
            }
            clu2elmset_map_t::iterator jt = ptnIt->second.find( cluIx );
            if ( jt != ptnIt->second.end() ) {
                jt->second.insert( elmIx );
            } else {
                elm_set_t eset;
                eset.insert( elmIx );
                ptnIt->second.insert( jt, clu2elmset_map_t::value_type( cluIx, eset ) );
            }
            ptnToTmplPtnMap[ ptnIx ] = tmplPtnId[ i ];
        }
    }

    // calculate co-clustered and not co-clustered pairs
    Rprintf( "Calculating (%d partitions)...\n", ptnColl.size() );
    Rcpp::DataFrame res;
    {
        Rcpp::StringVector ptnIdVec( ptnColl.size() );
        Rcpp::StringVector tmplPtnIdVec( ptnColl.size() );
        Rcpp::IntegerVector nCoCoVec( ptnColl.size() );
        Rcpp::IntegerVector nCoMmVec( ptnColl.size() );
        Rcpp::IntegerVector nMmCoVec( ptnColl.size() );
        Rcpp::IntegerVector nMmMmVec( ptnColl.size() );
        size_t i = 0;
        for ( ptn_coll_t::const_iterator ptnIt = ptnColl.begin(); ptnIt != ptnColl.end(); ++ptnIt ) {
            ptn_map_t::const_iterator tmplIdIt = ptnToTmplPtnMap.find( ptnIt->first );
            if ( tmplIdIt == ptnToTmplPtnMap.end() ) {
                THROW_RUNTIME_ERROR( "Template ID not found for partition "
                                        << "'" << ptnIt->first << "'" );
            }
            ptn_coll_t::const_iterator tmplIt = tmplPtnColl.find( tmplIdIt->second );
            if ( tmplIt == tmplPtnColl.end() ) {
                THROW_RUNTIME_ERROR( "Template partition ID=" << tmplIdIt->second
                                        << " not found for partition "
                                        << "'" << tmplIdIt->first << "'" );
            }
            const clu2elmset_map_t& tmplPtn = tmplIt->second;

            ptn_ix_t ptnIx = ptnIt->first;
            size_t nCoCo = 0;
            size_t nCoMm = 0;
            size_t nMmCo = 0;
            size_t nPtnElems = 0;
            for ( clu2elmset_map_t::const_iterator tmplCluIt = tmplPtn.begin();
                    tmplCluIt != tmplPtn.end(); ++tmplCluIt
            ){
                nPtnElems += tmplCluIt->second.size();
                const elm_set_t& tmplClu = tmplCluIt->second;
                size_t tmplCluSize = tmplClu.size();
                for ( clu2elmset_map_t::const_iterator cluIt = ptnIt->second.begin();
                    cluIt != ptnIt->second.end(); ++cluIt
                ){
                    const elm_set_t& clu = cluIt->second;
                    size_t cluSize = clu.size();
                    std::vector<elm_ix_t> isect;
                    std::set_intersection( clu.begin(), clu.end(),
                                           tmplClu.begin(), tmplClu.end(),
                                           std::inserter( isect, isect.end() ) );
                    size_t isectSize = isect.size();
                    size_t coco = isectSize * ( isectSize - 1 );
                    nCoCo += coco;
                    size_t comm = ( cluSize - isectSize ) * isectSize;
                    nCoMm += comm; // co-clustered in collection, but not in template
                    size_t mmco = ( tmplCluSize - isectSize ) * isectSize;
                    nMmCo += mmco; // co-clustered in template, but mismatch in collection
                }
            }
            if ( nPtnElems > nElms ) {
                THROW_RUNTIME_ERROR( "Partition '" << ptnIt->first << 
                                     "' has " << nPtnElems << " elements, > " << nElms );
            }
            tmplPtnIdVec[ i ] = tmplIdIt->second;
            ptnIdVec[ i ] = ptnIx;
            nCoCoVec[ i ] = nCoCo / 2;
            nCoMmVec[ i ] = nCoMm / 2;
            nMmCoVec[ i ] = nMmCo / 2;
            nMmMmVec[ i ] = ( nPtnElems * ( nElms - 1 ) - nCoCo - nCoMm - nMmCo ) / 2;
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