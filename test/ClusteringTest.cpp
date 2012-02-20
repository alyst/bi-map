#include <gtest/gtest.h>

#include <boost/functional/hash.hpp>
#include <boost/lexical_cast.hpp>

#include <cemm/math/PitmanYorProcess.h>
#include <cemm/math/Distributions.h>

#include <cemm/mcmc/SingleElementMembershipStep.h>
#include <cemm/mcmc/SplitMergeStep.h>

#include <cemm/containers/LazyClone.h>
#include <cemm/containers/EntityIndexing.h>
#include <cemm/containers/EntityCollectionIndexing.h>
#include <cemm/containers/set_misc.h>

using namespace cemm::math;
using namespace cemm::containers;
using namespace cemm::mcmc;

struct Point {
    double x;
    double y;
    double z;
};

typedef std::set<int> element_set_t;

typedef EntityCollectionIndexing<element_set_t> ClusteringIndexing;
typedef ClusteringIndexing::collection_pointer_type IndexedClusteringPointer;

struct Cluster {
    element_set_t elements;
    Point center;
};

struct Clustering {
    std::vector<Cluster> clusters;
    std::vector<int> pointToCluster;
    IndexedClusteringPointer  pIndex;
    
    const IndexedClusteringPointer& index( ClusteringIndexing& indexing )
    {
        //ASSERT_TRUE( (bool)pIndex );
        ClusteringIndexing::collection_type    value;
        for ( size_t i = 0; i < clusters.size(); i++ ) {
            value.insert( indexing.elementIndexing().index( clusters[ i ].elements ) );
        }
        return ( pIndex = indexing.index( value ) );
    }

    void check() const {
        for ( size_t i = 0; i < pointToCluster.size(); i++ ) {
            int cluIx = pointToCluster[ i ];
            if ( cluIx < 0 || cluIx >= (int)clusters.size() ) {
                throw std::runtime_error( "Invalid cluster index" );
            }
            if ( clusters[ cluIx ].elements.find( i ) == clusters[ cluIx ].elements.end() ) {
                throw std::runtime_error( "Point's cluster not found" );
            }
        }
    }

    void updateClusterIndexes() {
        std::vector<int> cluIxMap( clusters.size() );
        int oldCluIx = 0;
        for ( size_t newCluIx = 0; newCluIx < clusters.size(); oldCluIx++ ) {
            if ( clusters[ newCluIx ].elements.empty() ) {
                cluIxMap[ oldCluIx ] = -1;
                clusters.erase( clusters.begin() + newCluIx );
            }
            else {
                cluIxMap[ oldCluIx ] = newCluIx++;
            }
        }
        if ( cluIxMap.size() > clusters.size() ) for ( size_t i = 0; i < pointToCluster.size(); i++ ) {
            int& cluIx = pointToCluster[ i ];
            BOOST_ASSERT( cluIxMap[ cluIx ] >= 0 );
            cluIx = cluIxMap[ cluIx ];
        }
        check();
    }

    void putToCluster( int elmIx, int cluIx )
    {
        int oldCluIx = pointToCluster[ elmIx ];
        if ( oldCluIx != cluIx ) {
            clusters[ oldCluIx ].elements.erase( elmIx );
            if ( cluIx == (int)clusters.size() ) {
                clusters.push_back( Cluster() );
            }
            clusters[ cluIx ].elements.insert( elmIx );
            pointToCluster[ elmIx ] = cluIx;
         }
    }
};

struct Data {
    std::vector<Point> points;
};

class Partition: public LazyClone<Clustering> {
public:
    typedef int element_index_type;
    typedef std::set<element_index_type> element_index_set_type;

    typedef int cluster_index_type;
    typedef std::set<cluster_index_type> cluster_index_set_type;

    typedef Point& cluster_params_proxy_type;
    typedef Point cluster_params_type;

protected:
    typedef LazyClone<Clustering> super_type;

public:
    static const element_index_type ELEMENT_NA = -1;

    class ElementsSetProxy: public LazyClone<element_index_set_type> {
    protected:
        typedef LazyClone<element_index_set_type> super_type;

    public:
        ElementsSetProxy( const element_set_t& elements )
        : super_type( elements )
        {
        }
        ElementsSetProxy( const Partition& ptn, element_index_type elmIx = ELEMENT_NA )
        : super_type( boost::make_shared<element_set_t>( element_set_t() ), false )
        {
            if ( elmIx != ELEMENT_NA ) wrapped().insert( elmIx );
        }

        size_t size() const {
            return ( wrapped().size() );
        }

        element_index_type operator[]( int index ) const {
            element_index_set_type::const_iterator it = wrapped().begin();
            std::advance( it, index );
            return ( it != wrapped().end() ? *it : -1 );
        }
        bool contains( element_index_type elmIx ) const {
            return ( wrapped().find( elmIx ) != wrapped().end() );
        }
        ElementsSetProxy& operator+=( const ElementsSetProxy& that )
        {
            unshare();
            const locked_ptr_type thatpElms = that.lock();
            wrapped().insert( thatpElms->begin(), thatpElms->end() );
            return ( *this );
        }
        ElementsSetProxy& operator+=( element_index_type elmIx )
        {
            unshare();
            wrapped().insert( elmIx );
            return ( *this );
        }
        ElementsSetProxy& operator-=( const ElementsSetProxy& that )
        {
            unshare();
            const locked_ptr_type thatpElms = that.lock();
            wrapped().erase( thatpElms->begin(), thatpElms->end() );
            return ( *this );
        }
        ElementsSetProxy& operator-=( element_index_type elmIx )
        {
            unshare();
            wrapped().erase( elmIx );
            return ( *this );
        }
    };

    typedef ElementsSetProxy elements_set_proxy_type;
    class ClusterProxy {
    private:
        friend class Partition;

        Clustering& clus;
        const int cluIx;

        const Cluster& cluster() const {
            return ( clus.clusters[ cluIx ] );
        }

        ClusterProxy( Clustering& clus, int cluIx )
        : clus( clus ), cluIx( cluIx )
        {
            BOOST_ASSERT( cluIx >= 0 );
            BOOST_ASSERT( cluIx < (int)clus.clusters.size() );
        }

    public:
        const size_t size() const {
            return ( cluster().elements.size() );
        }
        const elements_set_proxy_type items() const {
            return ( ElementsSetProxy( cluster().elements ) );
        }
        const Point& params(
            const element_index_set_type& mask = element_index_set_type()
        ) const {
            return ( cluster().center );
        }
        cluster_params_proxy_type params(
            const element_index_set_type& mask = element_index_set_type()
        ){
            return ( clus.clusters[ cluIx ].center );
        }
        std::string label() const {
            return ( boost::lexical_cast<std::string>( cluIx ) );
        }
        cluster_index_type index() const {
            return ( cluIx );
        }
    };

    typedef ClusterProxy cluster_proxy_type;

    const static cluster_index_type ClusterNA = -1;

    Partition( Clustering& clustering )
    : super_type( clustering )
    {
    }

    virtual ~Partition()
    {}

    size_t clustersCount() const {
        return ( wrapped().clusters.size() );
    }

    size_t elementsCount() const {
        return ( wrapped().pointToCluster.size() );
    }

    cluster_index_type clusterIndex( element_index_type elmIx ) const {
        return ( wrapped().pointToCluster[ elmIx ] );
    }

    const cluster_proxy_type cluster( cluster_index_type cluIx ) const {
        return ( ClusterProxy( const_cast<Clustering&>( wrapped() ), cluIx ) );
    }

    cluster_proxy_type cluster( cluster_index_type cluIx ) {
        return ( ClusterProxy( wrapped(), cluIx ) );
    }

    void mergeClusters( int clu1Ix, int clu2Ix )
    {
        SCOPED_TRACE( "mergeClusters()" );
        if ( clu1Ix > clu2Ix ) std::swap( clu1Ix, clu2Ix );
        unshare();
        const element_set_t& clu2Elms = wrapped().clusters[ clu2Ix ].elements;
        wrapped().clusters[ clu1Ix ].elements.insert( clu2Elms.begin(), clu2Elms.end() );

        for ( size_t k = 0; k < wrapped().pointToCluster.size(); k++ ) {
            int& cluIx = wrapped().pointToCluster[ k ];
            if ( cluIx == clu2Ix ) {
                cluIx = clu1Ix;
            }
        }
        wrapped().clusters[ clu2Ix ].elements = element_set_t();
        EXPECT_NO_THROW( wrapped().check() );
    }

    int exchangeElements(
        cluster_index_type clu1Ix, 
        cluster_index_type clu2Ix, 
        const element_index_set_type& clu2Elements
    ){
        SCOPED_TRACE( "exchangeElements()" );
        unshare();
        cluster_index_type res = clu2Ix;
        for ( element_index_set_type::const_iterator it = clu2Elements.begin(); it != clu2Elements.end(); ++it ) {
            int& cluIx = wrapped().pointToCluster[ *it ]; 
            if ( cluIx == clu1Ix ) {
                if ( res == ClusterNA ) {
                    res = wrapped().clusters.size();
                    wrapped().clusters.push_back( Cluster() );
                }
                wrapped().clusters[ res ].elements.insert( *it );
                wrapped().clusters[ clu1Ix ].elements.erase( *it );
                cluIx = res;
            }
            else {
                BOOST_ASSERT( clu2Ix == cluIx );
            }
        }
        if ( clu2Ix != ClusterNA ) {
            element_set_t newClu1Items = wrapped().clusters[ clu2Ix ].elements;
            for ( element_set_t::const_iterator it = clu2Elements.begin(); it != clu2Elements.end(); ++it ) {
                newClu1Items.erase( *it );
            }
            for ( element_set_t::const_iterator it = newClu1Items.begin(); it != newClu1Items.end(); ++it ) {
                int& cluIx = wrapped().pointToCluster[ *it ]; 
                //ASSERT_EQ( clu2Ix, cluIx );
                wrapped().clusters[ clu2Ix ].elements.erase( *it );
                wrapped().clusters[ clu1Ix ].elements.insert( *it );
                cluIx = clu1Ix;
            }
        }
        EXPECT_NO_THROW( wrapped().check() );
        return ( res );
    }

    void putToCluster(
        element_index_type              elmIx,
        cluster_index_type              cluIx
    ){
        SCOPED_TRACE( "putToCuster()" );
        int oldCluIx = clusterIndex( elmIx );
        if ( oldCluIx == cluIx ) return;
        unshare();
        wrapped().putToCluster( elmIx, cluIx );
        EXPECT_NO_THROW( wrapped().check() );
    }
};

const double CluPointsSigma = 0.5;
const double CluCentersSigma = 4.0;

struct PartitionStats {
    const Data& data;
    const PitmanYorProcess& clusPrior;

    struct GroupSizes {
        const Partition& ptn;

        GroupSizes( const Partition& ptn )
        : ptn( ptn )
        {}

        size_t size() const {
            return ( ptn.clustersCount() );
        }
        size_t operator[]( size_t cluIx ) const {
            return ( ptn.cluster( cluIx ).size() );
        }
    };

    double llh( const Partition& ptn, const Partition::element_index_type elmIx, const Point& params ) const
    {
        return ( llh( elmIx, params ) );
    }

    double llh( const Partition::element_index_type elmIx, const Point& params ) const
    {
        double res = 0;
        const Point& objPt = data.points[ elmIx ];
        res += gaussian_pdf_ln( objPt.x - params.x, CluPointsSigma );
        res += gaussian_pdf_ln( objPt.y - params.y, CluPointsSigma );
        res += gaussian_pdf_ln( objPt.z - params.z, CluPointsSigma );
        return ( res );
    }

    double llh( const Partition& clus, const Partition::element_index_set_type& objs, const Point& params ) const
    {
        double res = 0;
        for ( Partition::element_index_set_type::const_iterator objIt = objs.begin(); objIt != objs.end(); ++objIt ) {
            res += llh( *objIt, params );
        }
        return ( res );
    }

    double llhDelta( const Partition& clus,
                const std::vector<Partition::elements_set_proxy_type>& newClusters,
                const std::vector<Point>& newParams,
                const Partition::cluster_index_set_type& oldClusters
    ) const {
        double res = 0;
        for ( size_t newCluIx = 0; newCluIx < newClusters.size(); newCluIx++ ) {
            const Partition::element_index_set_type objs = newClusters[ newCluIx ];
            for ( Partition::element_index_set_type::const_iterator objIt = objs.begin(); objIt != objs.end(); ++objIt ) {
                res += llh( *objIt, newParams[ newCluIx ] );
            }
        }
        for ( Partition::cluster_index_set_type::const_iterator oldCluIxIt = oldClusters.begin();
              oldCluIxIt != oldClusters.end(); ++oldCluIxIt
        ){
            const Partition::element_index_set_type objs = clus.cluster( *oldCluIxIt ).items();
            for ( Partition::element_index_set_type::const_iterator objIt = objs.begin(); objIt != objs.end(); ++objIt ) {
                res += llh( *objIt, clus.cluster( *oldCluIxIt ).params() );
            }
        }
        return ( res );
    }

    double llh( const Partition& clus ) const
    {
        double res = 0;
        for ( int i = 0; i < (int)clus.clustersCount(); i++ ) {
            const Partition::ClusterProxy clu = clus.cluster( i );
            res += llh( clus, clu.items(), clu.params() );
        }
        return ( res );
    }

    double lpp( const Partition& clus ) const
    {
        double res = clusPrior.lnP( GroupSizes( clus ) );
        for ( Partition::cluster_index_type cluIx = 0; cluIx < (int)clus.clustersCount(); cluIx++ ) {
            const Point& center = clus.cluster( cluIx ).params();
            res += gaussian_pdf_ln( center.x - 0, CluCentersSigma );
            res += gaussian_pdf_ln( center.y - 0, CluCentersSigma );
            res += gaussian_pdf_ln( center.z - 0, CluCentersSigma );
        }
        return ( res );
    }

    double paramsLPP( const Partition& clust, const Point& params,
                      const Partition::element_index_set_type& clusterObjs,
                      const Partition::element_index_set_type& sampledObjs ) const
    {
        double res = 0;
        res += gaussian_pdf_ln( params.x - 0, CluCentersSigma );
        res += gaussian_pdf_ln( params.y - 0, CluCentersSigma );
        res += gaussian_pdf_ln( params.z - 0, CluCentersSigma );
        return ( res );
    }

    double clusterAssignmentLPP( size_t clusterSize, size_t clustersCount, size_t totalElements ) const
    {
        return ( log( clusPrior.clusterAssignmentPrior( clusterSize, clustersCount, totalElements - 1 ) ) );
    }

    PartitionStats( const Data& data, const PitmanYorProcess& clusPrior )
    : data( data )
    , clusPrior( clusPrior )
    {}
};

struct ClusterSampler {
    typedef Point params_type;

    const gsl_rng* rng;
    const Data& data;

    Point defaults() const {
        Point a;
        a.x = a.y = a.z = 0.0;
        return ( a );
    }

    Point defaults( const Partition& ptn ) const {
        return ( defaults() );
    }

    bool operator()( Point& params, 
                     const Partition::element_index_set_type& clusterObjects, 
                     const Partition::element_index_set_type& sampledObjects, 
                     bool overwrite, bool posterior ) const
    {
        if ( overwrite ) {
            if ( posterior ) {
                Partition::element_index_set_type allObjs( clusterObjects );
                allObjs.insert( sampledObjects.begin(), sampledObjects.end() );
                Point sum;
                sum.x = sum.y = sum.z = 0.0;
                //Point disp = sum;
                for ( Partition::element_index_set_type::const_iterator it = allObjs.begin(); it != allObjs.end(); ++it ) {
                    const Point& objPt = data.points[ *it ];
                    sum.x += objPt.x;
                    sum.y += objPt.y;
                    sum.z += objPt.z;
                    //disp.x += objPt.x * objPt.x;
                    //disp.y += objPt.y * objPt.y;
                    //disp.z += objPt.z * objPt.z;
                }
//                 Point mean;
//                 mean.x = sum.x / allObjs.size();
//                 mean.y = sum.y / allObjs.size();
//                 mean.z = sum.z / allObjs.size();
//                 disp.x /= allObjs.size();
//                 disp.y /= allObjs.size();
//                 disp.z /= allObjs.size();
//                 disp.x -= mean.x * mean.x;
//                 disp.y -= mean.y * mean.y;
//                 disp.z -= mean.z * mean.z;

//                 Point postMean;
                Point postInvDisp;
                postInvDisp.x = ( allObjs.size() * CluCentersSigma * CluCentersSigma + CluPointsSigma * CluPointsSigma );
                postInvDisp.y = ( allObjs.size() * CluCentersSigma * CluCentersSigma + CluPointsSigma * CluPointsSigma );
                postInvDisp.z = ( allObjs.size() * CluCentersSigma * CluCentersSigma + CluPointsSigma * CluPointsSigma );
                Point postMean;
                postMean.x = ( 0 + sum.x * CluCentersSigma * CluCentersSigma ) / postInvDisp.x;
                postMean.y = ( 0 + sum.y * CluCentersSigma * CluCentersSigma ) / postInvDisp.y;
                postMean.z = ( 0 + sum.z * CluCentersSigma * CluCentersSigma ) / postInvDisp.z;
                Point postDisp;
                postDisp.x = ( CluCentersSigma * CluCentersSigma * CluPointsSigma * CluPointsSigma ) / postInvDisp.x;
                postDisp.y = ( CluCentersSigma * CluCentersSigma * CluPointsSigma * CluPointsSigma ) / postInvDisp.y;
                postDisp.z = ( CluCentersSigma * CluCentersSigma * CluPointsSigma * CluPointsSigma ) / postInvDisp.z;
                params.x = postMean.x + gsl_ran_gaussian( rng, sqrt( postDisp.x ) );
                params.y = postMean.y + gsl_ran_gaussian( rng, sqrt( postDisp.y ) );
                params.z = postMean.z + gsl_ran_gaussian( rng, sqrt( postDisp.z ) );
            }
            else {
                params.x = gsl_ran_gaussian( rng, CluCentersSigma );
                params.y = gsl_ran_gaussian( rng, CluCentersSigma );
                params.z = gsl_ran_gaussian( rng, CluCentersSigma );
            }
            return ( true );
        }
        return ( false );
    }

    bool operator()( const Partition& ptn, Point& params, 
                     const Partition::element_index_set_type& clusterObjects, 
                     const Partition::element_index_set_type& sampledObjects, 
                     bool overwrite, bool posterior ) const
    {
        return ( operator()( ptn, params, clusterObjects, sampledObjects, 
                             overwrite, posterior ) );
    }

    double transitionLP( const Partition& before, const Partition& after, size_t cluIx ) const
    {
        PartitionStats stats( data, PitmanYorProcess( 0.5, 0.5 ) );
        const Partition::ClusterProxy cluBefore = before.cluster( cluIx );
        const Partition::ClusterProxy cluAfter = after.cluster( cluIx );
        double res = gaussian_pdf_ln( cluAfter.params().x - 0, CluCentersSigma )
                   - gaussian_pdf_ln( cluBefore.params().x - 0, CluCentersSigma )
                   + gaussian_pdf_ln( cluAfter.params().y - 0, CluCentersSigma )
                   - gaussian_pdf_ln( cluBefore.params().y - 0, CluCentersSigma )
                   + gaussian_pdf_ln( cluAfter.params().z - 0, CluCentersSigma )
                   - gaussian_pdf_ln( cluBefore.params().z - 0, CluCentersSigma )
                   + stats.llh( after, cluBefore.items(), cluAfter.params() )
                   - stats.llh( before, cluBefore.items(), cluBefore.params() );
        return ( res );
    }

    ClusterSampler( const gsl_rng* rng, const Data& data )
    : rng( rng )
    , data( data )
    {}
};

TEST( MCMC_Clustering, DISABLED_generic )
{
    gsl_rng* rng = gsl_rng_alloc( gsl_rng_default );
    Clustering origModel;
    PitmanYorProcess    cluPrior( 0.3, 0.2 ); /// prior for split-merge steps and random generation
    PitmanYorProcess    cluDPrior( 0.01, 0.0 ); /// modified prior is used for single element sampling

    const int PointsCount = 100;

    // generate model
    origModel.pointToCluster.resize( PointsCount );
    PitmanYorSample priorClu = cluPrior.random( rng, PointsCount );
    for ( int cluIx = 0; cluIx < (int)priorClu.clustersCount(); cluIx++ ) {
        Cluster clu;
        clu.elements.insert( priorClu.cluster( cluIx ).begin(), priorClu.cluster( cluIx ).end() );
        clu.center.x = gsl_ran_gaussian( rng, CluCentersSigma );
        clu.center.y = gsl_ran_gaussian( rng, CluCentersSigma );
        clu.center.z = gsl_ran_gaussian( rng, CluCentersSigma );

        origModel.clusters.push_back( clu );
        for ( PitmanYorSample::sample_set_type::const_iterator it = priorClu.cluster( cluIx ).begin(); it != priorClu.cluster( cluIx ).end(); ++it ) {
            origModel.pointToCluster[ *it ] = cluIx;
        }
    }
    EXPECT_NO_THROW( origModel.check() );
    {
        std::ostringstream msg;
        msg << "Generated clusters: " << origModel.clusters.size() << ": ";
        for ( int cluIx = 0; cluIx < (int)origModel.clusters.size(); cluIx++ ) {
            const Cluster& clu = origModel.clusters[ cluIx ];
            msg << clu.elements.size() << "(x=" << clu.center.x << " y=" << clu.center.y << " z=" << clu.center.z << ") ";
        }
        LOG_INFO( msg.str() );
    }

    // generate data
    Data data;
    data.points.resize( PointsCount );
    for ( int ptIx = 0; ptIx < PointsCount; ptIx++ ) {
        int cluIx = origModel.pointToCluster[ ptIx ];
        const Cluster& clu = origModel.clusters[ cluIx ];
        data.points[ ptIx ].x = clu.center.x + gsl_ran_gaussian( rng, CluPointsSigma );
        data.points[ ptIx ].y = clu.center.y + gsl_ran_gaussian( rng, CluPointsSigma );
        data.points[ ptIx ].z = clu.center.z + gsl_ran_gaussian( rng, CluPointsSigma );
    }

    ClusteringIndexing  indexing( 10000, 30000 );
    std::vector<Clustering> walk;

    {
    // generate initial clustering
    Clustering clustering;
    clustering.pointToCluster.resize( PointsCount, 0 );
    clustering.clusters.resize( 1 );
    clustering.clusters[0].center.x = clustering.clusters[0].center.y = 
    clustering.clusters[0].center.z = 0.0;
    for ( size_t i = 0; i < clustering.pointToCluster.size(); i++ ) {
        clustering.clusters[0].elements.insert( i );
    }
    EXPECT_NO_THROW( clustering.check() );

    const int IterationsCount = 5000;
    LOG_INFO( "Making " << IterationsCount << " iterations..." );
    for ( int i = 0; i < IterationsCount; i++ ) {
        if ( i % 1 == 0 ) {
            int elmIx1 = gsl_rng_uniform_int( rng, PointsCount );
            int elmIx2;
            while ( ( elmIx2 = gsl_rng_uniform_int( rng, PointsCount ) ) == elmIx1 ) {};
            LOG_DEBUG2_IF( clustering.pointToCluster[ elmIx1 ] == clustering.pointToCluster[ elmIx2 ],
                           "Currently merged elements " << elmIx1 << " and " << elmIx2 );
            LOG_DEBUG2_IF( origModel.pointToCluster[ elmIx1 ] != origModel.pointToCluster[ elmIx2 ],
                           "Originally split elements " << elmIx1 << " and " << elmIx2 );
            typedef SplitMergeSamplingStep<Partition, PartitionStats, ClusterSampler> ClusteringSplitMerge;
            ClusteringSplitMerge csm( rng, PartitionStats( data, cluPrior ), ClusterSampler( rng, data ) );
            ClusteringSplitMerge::result_type res = csm( Partition( clustering ), elmIx1, elmIx2 );
            if ( !res.ptn.isRef() ) {
                clustering = res.ptn;
            }
            clustering.updateClusterIndexes();
        }
        if ( ( i % 5 == 3 ) ) {
            for ( size_t j = 0; j < data.points.size(); j++ ) 
            {
                //int j = gsl_rng_uniform_int( rng, clustering.pointToCluster.size() );
                Point   newCenter;
                int curCluIx = clustering.pointToCluster[ j ];
                ASSERT_GE( curCluIx, 0 );
                PartitionStats ps( data, cluDPrior );
                int newCluIx = ::SampleClusterOfElement( rng, Partition( clustering ), ps,
                                                         j, ClusterSampler( rng, data ),
                                                         ClusterOfElementStepParams( 2, 2 ),
                                                         SamplingTransform(), ps.llh( clustering ) + ps.lpp( clustering ),
                                                         &newCenter );
                ASSERT_GE( newCluIx, 0 );
                if ( newCluIx != curCluIx ) {
                    clustering.putToCluster( j, newCluIx );
                    clustering.clusters[ newCluIx ].center = newCenter;
                    clustering.updateClusterIndexes();
                }
            }
        }
        #if 0
        for ( int cluIx = 0; cluIx < clustering.clusters.size(); cluIx++ ) {
            ClusterSampler( rng, data )( clustering.clusters[ cluIx ].center, 
                                                    clustering.clusters[ cluIx ].elements, 
                                                    clustering.clusters[ cluIx ].elements, true, true );
        }
        #endif
        if ( ( i > 2000 ) && ( ( ( i + 1 ) % 20 ) == 0 ) ) {
            Clustering copy = clustering;
            copy.index( indexing );
            walk.push_back( copy );
        }
    }
    }
    
    LOG_INFO( "Recorded " << walk.size() << " clusterings" );
    int  bestCluIx = -1;
    size_t bestCluRefCnt = 0;
    for ( size_t i = 0; i < walk.size(); i++ ) {
        if ( walk[ i ].pIndex->refCount() >= bestCluRefCnt ) {
            bestCluIx = i;
            bestCluRefCnt = walk[ i ].pIndex->refCount();
        }
    }
    LOG_INFO( "Best clustering " << bestCluIx << " referenced " << bestCluRefCnt << " times" );

    const Clustering& bestClu = walk[ bestCluIx ];

    {
        std::ostringstream msg;
        msg << "Inferred clusters: " << bestClu.clusters.size() << ": ";
        for ( size_t cluIx = 0; cluIx < bestClu.clusters.size(); cluIx++ ) {
            const Cluster& clu = bestClu.clusters[ cluIx ];
            msg << clu.elements.size() << "(x=" << clu.center.x << " y=" << clu.center.y << " z=" << clu.center.z << ") ";
        }
        LOG_INFO( msg.str() );
    }

    size_t matches = 0;
    size_t wrongMerge = 0;
    size_t wrongSplit = 0;
    for ( int i = 0; i < PointsCount; i++ ) {
        int origCluI = origModel.pointToCluster[ i ];
        int cluI = bestClu.pointToCluster[ i ];
        for ( int j = i + 1; j < PointsCount; j++ ) {
            int origCluJ = origModel.pointToCluster[ j ];
            int cluJ = bestClu.pointToCluster[ j ];
            if ( ( cluI == cluJ ) && ( origCluJ == origCluI ) ) {
                //LOG_INFO( "match " << cluI << " - " << origCluJ );
                matches++;
            }
            else if ( ( cluI == cluJ ) && ( origCluJ != origCluI ) 
            ){
                //LOG_INFO( "w.merge " << cluI << " - " << origCluJ );
                wrongMerge++;
            }
            else if ( ( cluI != cluJ ) && ( origCluJ == origCluI ) 
            ){
                //LOG_INFO( "w.split " << cluI << " - " << origCluJ );
                wrongSplit++;
            }
        }
    }
    LOG_INFO( "Matching pairs: " << matches );
    LOG_INFO( "Wrong merges: " << wrongMerge );
    LOG_INFO( "Wrong splits: " << wrongSplit );
    int totalOrigPairs = 0;
    for ( size_t i = 0; i < origModel.clusters.size(); i++ ) {
        int cluSize = origModel.clusters[i].elements.size();
        totalOrigPairs += cluSize * ( cluSize - 1 ) / 2;
    }
    EXPECT_GE( (double) matches / totalOrigPairs, 0.75 );
    EXPECT_LE( (double) wrongMerge / totalOrigPairs, 0.25 );
    EXPECT_LE( (double) wrongSplit / totalOrigPairs, 0.25 );
    EXPECT_GE( (double) bestCluRefCnt / walk.size(), 0.33 );
    gsl_rng_free( rng );
}
