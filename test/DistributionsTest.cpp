#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

#include <math/Distributions.h>
#include <math/LagrangianPoisson.h>
#include <math/GenericDiscreteDistribution.h>

#include "TestCommon.h"

#define LN_PDF_EPS ( 1E+4 * std::numeric_limits<log_prob_t>::epsilon() )

TEST( LogMath, ln_binomial_coef )
{
    EXPECT_NEAR( ln_binomial_coef( 3000, 1 ), log( 3000 ), LN_PDF_EPS );
}

TEST( Distributions, Geometric )
{
    GeometricDistribution a = GeometricDistribution::ByFailureRate( 0.1 );
    EXPECT_NEAR( GeometricDistribution::ByFailureRate( 0.1 ).lnPdf( 3 ),
                 GeometricDistribution::BySuccessRate( 0.9 ).lnPdf( 3 ), LN_PDF_EPS );
    EXPECT_NEAR( GeometricDistribution::ByFailureRate( 0.9 ).lnCdf_P( 4 ),
                 GeometricDistribution::BySuccessRate( 0.1 ).lnCdf_P( 4 ), LN_PDF_EPS );
    EXPECT_NEAR( GeometricDistribution::ByFailureRate( 0.9 ).lnCdf_Q( 4 ),
                 GeometricDistribution::BySuccessRate( 0.1 ).lnCdf_Q( 4 ),
                 LN_PDF_EPS );
    EXPECT_EQ( 1.0,
               exp( GeometricDistribution::ByFailureRate( 0.1 ).lnCdf_P( 4 ) ) +
               exp( GeometricDistribution::ByFailureRate( 0.1 ).lnCdf_Q( 4 ) ) );
    EXPECT_NEAR( gsl_ran_geometric_pdf( 3, 0.1 ),
                 exp( GeometricDistribution::BySuccessRate( 0.1, 1 ).lnPdf( 3 ) ),
                 LN_PDF_EPS );
    EXPECT_EQ( gsl_ran_geometric_pdf( 0, 0.1 ),
               exp( GeometricDistribution::BySuccessRate( 0.1, 1 ).lnPdf( 0 ) ) );
    EXPECT_NEAR( gsl_cdf_geometric_P( 3, 0.1 ),
                 exp( GeometricDistribution::BySuccessRate( 0.1, 0 ).lnCdf_P( 2 ) ),
                 LN_PDF_EPS );
    EXPECT_NEAR( gsl_cdf_geometric_Q( 3, 0.1 ),
                 exp( GeometricDistribution::BySuccessRate( 0.1, 0 ).lnCdf_Q( 2 ) ),
                 LN_PDF_EPS );
    EXPECT_NEAR( GeometricDistribution::ByFailureRate( 0.1, 1 ).lnPdf( 3 ),
                 GeometricDistribution::ByFailureRate( 0.1, 0 ).lnPdf( 2 ), LN_PDF_EPS );
    EXPECT_NEAR( GeometricDistribution::BySuccessRate( 0.1, 2 ).lnPdf( 4 ),
                 GeometricDistribution::ByFailureRate( 0.9, 1 ).lnPdf( 3 ), LN_PDF_EPS );
}

TEST( Distributions, Hypergeometric )
{
    EXPECT_NEAR( log( gsl_ran_hypergeometric_pdf( 1, 3, 7, 3 ) ),
                 HypergeometricDistribution( 3, 7, 3 ).lnPdf( 1 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_hypergeometric_Q( 1, 3, 7, 3 ) ),
                 HypergeometricDistribution( 3, 7, 3 ).lnCdf_Q( 1 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_hypergeometric_Q( 1, 3, 7, 3 ) + gsl_ran_hypergeometric_pdf( 1, 3, 7, 3 ) ),
                 HypergeometricDistribution( 3, 7, 3 ).lnCdf_Q( 1, true ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_hypergeometric_Q( 10, 30, 70, 30 ) ),
                 HypergeometricDistribution( 30, 70, 30 ).lnCdf_Q( 10 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_hypergeometric_Q( 20, 30, 70, 30 ) ),
                 HypergeometricDistribution( 30, 70, 30 ).lnCdf_Q( 20 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_hypergeometric_Q( 25, 30, 70, 30 ) ),
                 HypergeometricDistribution( 30, 70, 30 ).lnCdf_Q( 25 ), LN_PDF_EPS );
    EXPECT_NEAR( log( 1.0 / 100 ),
                 HypergeometricDistribution( 1, 99, 1 ).lnPdf( 1 ), LN_PDF_EPS );
    EXPECT_NEAR( log( 99.0 / 100 ),
                 HypergeometricDistribution( 1, 99, 1 ).lnPdf( 0 ), LN_PDF_EPS );
    EXPECT_NEAR( log( 98.0 / 100 ),
                 HypergeometricDistribution( 1, 99, 2 ).lnPdf( 0 ), LN_PDF_EPS );
    // test extreme values
    EXPECT_NEAR( -log( 184756 ),
                 HypergeometricDistribution( 10, 10, 10 ).lnPdf( 0 ), LN_PDF_EPS );
    EXPECT_NEAR( -log( 184756 ),
                 HypergeometricDistribution( 10, 10, 10 ).lnPdf( 10 ), LN_PDF_EPS );
    EXPECT_EQ( -std::numeric_limits<double>::infinity(),
               HypergeometricDistribution( 10, 10, 10 ).lnCdf_P( 0, false ) );
    EXPECT_EQ( -std::numeric_limits<double>::infinity(),
               HypergeometricDistribution( 10, 10, 10 ).lnCdf_Q( 10, false ) );
    EXPECT_EQ( -log( 184756 ),
               HypergeometricDistribution( 10, 10, 10 ).lnCdf_Q( 10, true ) );
    EXPECT_EQ( -log( 184756 ),
               HypergeometricDistribution( 10, 10, 10 ).lnCdf_P( 0, true ) );
    EXPECT_NEAR( 0,
                 HypergeometricDistribution( 10, 10, 10 ).lnCdf_Q( 0, true ),
                 LN_PDF_EPS );
    EXPECT_NEAR( 0,
                 HypergeometricDistribution( 10, 10, 10 ).lnCdf_P( 10, true ), 
                 LN_PDF_EPS  );
}

TEST( Distributions, FisherHypergeometric )
{
    EXPECT_NEAR( log( gsl_ran_hypergeometric_pdf( 1, 3, 7, 3 ) ),
                 FisherHypergeometricDistribution( 3, 7, 3, 0 ).lnPdf( 1 ), LN_PDF_EPS );

    EXPECT_NEAR( -0.697107581776,
                 FisherHypergeometricDistribution( 3, 7, 3, log( 2.0 ) ).lnPdf( 1 ), LN_PDF_EPS );
    EXPECT_NEAR( -0.148894425938,
                 FisherHypergeometricDistribution( 3, 7, 3, log( 2.0 ) ).lnCdf_Q( 1, true ), LN_PDF_EPS );
    EXPECT_NEAR( -0.148894425938,
                 FisherHypergeometricDistribution( 3, 7, 3, log( 2.0 ) ).lnCdf_Q( 0, false ), LN_PDF_EPS );

    EXPECT_NEAR( -2.20678444172,
                 FisherHypergeometricDistribution( 30, 70, 30, log( 2.0 ) ).lnPdf( 10 ), 1E-10 );
    EXPECT_NEAR( -0.241945613981,
                 FisherHypergeometricDistribution( 30, 70, 30, log( 2.0 ) ).lnCdf_Q( 10, false ), LN_PDF_EPS );

    EXPECT_NEAR( -8.46215439545,
                 FisherHypergeometricDistribution( 30, 70, 30, log( 0.2 ) ).lnPdf( 10 ), 1E-10 );
    EXPECT_NEAR( -10.2898694358,
                 FisherHypergeometricDistribution( 30, 70, 30, log( 0.2 ) ).lnCdf_Q( 10, false ), 1E-8 );

    EXPECT_NEAR( -957.92570027089762461,
                 FisherHypergeometricDistribution( 300, 2700, 300, log( 1.0 ) ).lnPdf( 299 ), 1E-8 );
    EXPECT_NEAR( -957.92569903633048545,
                 FisherHypergeometricDistribution( 300, 2700, 300, log( 1.0 ) ).lnCdf_Q( 299, true ), 1E-8 );
    EXPECT_NEAR( -593.02779139450798255,
                 FisherHypergeometricDistribution( 300, 2700, 300, log( 1.0 ) ).lnCdf_Q( 250, true ), 1E-8 );

    // test extreme values
    EXPECT_NEAR( log( 9765625 ) - log( 3038019876 ),
                 FisherHypergeometricDistribution( 10, 10, 10, log( 0.2 ) ).lnPdf( 0 ), LN_PDF_EPS );
    EXPECT_NEAR( -log( 3038019876 ),
                 FisherHypergeometricDistribution( 10, 10, 10, log( 0.2 ) ).lnPdf( 10 ), LN_PDF_EPS );
    EXPECT_EQ( -std::numeric_limits<double>::infinity(),
               FisherHypergeometricDistribution( 10, 10, 10, log( 0.2 ) ).lnCdf_P( 0, false ) );
    EXPECT_EQ( -std::numeric_limits<double>::infinity(),
               FisherHypergeometricDistribution( 10, 10, 10, log( 0.2 ) ).lnCdf_Q( 10, false ) );
    EXPECT_NEAR( 0, FisherHypergeometricDistribution( 10, 10, 10, log( 0.2 ) ).lnCdf_Q( 0, true ), LN_PDF_EPS );
    EXPECT_NEAR( 0, FisherHypergeometricDistribution( 10, 10, 10, log( 0.2 ) ).lnCdf_P( 10, true ), LN_PDF_EPS );
    EXPECT_EQ( log( 9765625 ) - log( 3038019876 ),
               FisherHypergeometricDistribution( 10, 10, 10, log( 0.2 ) ).lnCdf_P( 0, true ) );
    EXPECT_EQ( -log( 3038019876 ),
               FisherHypergeometricDistribution( 10, 10, 10, log( 0.2 ) ).lnCdf_Q( 10, true ) );
}

TEST( Distributions, Beta )
{
    EXPECT_NEAR( log( gsl_ran_beta_pdf( 0.5, 2 + 1, 5 + 1 ) ),
        BetaDistribution( 2, 5 )( 0.5 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_beta_pdf( 0.1, 2 + 1, 5 + 1 ) ),
        BetaDistribution( 2, 5 )( 0.1 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_beta_pdf( 0.8, 2 + 1, 5 + 1 ) ),
        BetaDistribution( 2, 5 )( 0.8 ), LN_PDF_EPS );

    EXPECT_NEAR( log( gsl_ran_beta_pdf( 0.5, 2 + 1, 3 + 1 ) ),
        BetaDistribution( 2, 3 )( 0.5 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_beta_pdf( 0.1, 2 + 1, 3 + 1 ) ),
        BetaDistribution( 2, 3 )( 0.1 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_beta_pdf( 0.8, 2 + 1, 3 + 1 ) ),
        BetaDistribution( 2, 3 )( 0.8 ), LN_PDF_EPS );

    EXPECT_NEAR( log( gsl_ran_beta_pdf( 0.5, 1 + 1, 5 + 1 ) ),
        BetaDistribution( 1, 5 )( 0.5 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_beta_pdf( 0.1, 1 + 1, 5 + 1 ) ),
        BetaDistribution( 1, 5 )( 0.1 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_beta_pdf( 0.8, 1 + 1, 5 + 1 ) ),
        BetaDistribution( 1, 5 )( 0.8 ), LN_PDF_EPS );
}

TEST( Distributions, LagrangianPoisson )
{
    ln_factorial_cache_fill( 100 );

    EXPECT_NEAR( log( gsl_cdf_poisson_Q( 2, 5.0 ) ),
        LagrangianPoissonDistribution( log( 5.0 ) ).lnCdf_Q( 2 ), LN_PDF_EPS );
    EXPECT_NEAR( 1.0,
        exp( LagrangianPoissonDistribution( log( 5.0 ) ).lnCdf_Q( 2 ) ) +
        exp( LagrangianPoissonDistribution( log( 5.0 ) ).lnCdf_P( 2 ) ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_poisson_Q( 2, 3.0 ) ),
        LagrangianPoissonDistribution( log( 3.0 ) ).lnCdf_Q( 2 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_poisson_Q( 5, 5.0 ) ),
        LagrangianPoissonDistribution( log( 5.0 ) ).lnCdf_Q( 5 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_poisson_Q( 0, 3.0 ) ),
        LagrangianPoissonDistribution( log( 3.0 ) ).lnCdf_Q( 0 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_poisson_Q( 2, 3.0 ) ),
        LagrangianPoissonDistribution( log( 3.0 ) ).lnCdf_Q( 2 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_cdf_poisson_Q( 5, 3.0 ) ),
        LagrangianPoissonDistribution( log( 3.0 ) ).lnCdf_Q( 5 ), LN_PDF_EPS );

    EXPECT_NEAR( log( gsl_ran_poisson_pdf( 0, 5.0 ) ), 
        LagrangianPoissonDistribution( log( 5.0 ) ).lnPdf( 0 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_poisson_pdf( 2, 5.0 ) ), 
        LagrangianPoissonDistribution( log( 5.0 ) ).lnPdf( 2 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_poisson_pdf( 5, 5.0 ) ), 
        LagrangianPoissonDistribution( log( 5.0 ) ).lnPdf( 5 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_poisson_pdf( 0, 3.0 ) ), 
        LagrangianPoissonDistribution( log( 3.0 ) ).lnPdf( 0 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_poisson_pdf( 2, 3.0 ) ), 
        LagrangianPoissonDistribution( log( 3.0 ) ).lnPdf( 2 ), LN_PDF_EPS );
    EXPECT_NEAR( log( gsl_ran_poisson_pdf( 5, 3.0 ) ), 
        LagrangianPoissonDistribution( log( 3.0 ) ).lnPdf( 5 ), LN_PDF_EPS );

    EXPECT_NEAR( -3.6188758248672007492, LagrangianPoissonDistribution( log( 0.2 ), 0.1 ).lnPdf( 2 ), LN_PDF_EPS );
    EXPECT_NEAR( -0.2, LagrangianPoissonDistribution( log( 0.2 ), 0.1 ).lnPdf( 0 ), LN_PDF_EPS );
    EXPECT_NEAR( -1.0, LagrangianPoissonDistribution( log( 1.0 ), 0.1 ).lnPdf( 0 ), LN_PDF_EPS );

    EXPECT_NEAR( -10, LagrangianPoissonDistribution( log( 10 ), 0.1 ).lnPdf( 0 ), LN_PDF_EPS );
    EXPECT_NEAR( -5.1248865857628291433, LagrangianPoissonDistribution( log( 10 ), 0.1 ).lnPdf( 3 ), LN_PDF_EPS );
    EXPECT_NEAR( -28.480722004770693150, LagrangianPoissonDistribution( log( 10 ), 0.1 ).lnPdf( 50 ), LN_PDF_EPS );

    EXPECT_NEAR( -10, LagrangianPoissonDistribution( log( 10 ), 0.7 ).lnPdf( 0 ), LN_PDF_EPS );
    EXPECT_NEAR( -6.6027634710286185086, LagrangianPoissonDistribution( log( 10 ), 0.7 ).lnPdf( 3 ), LN_PDF_EPS );
    EXPECT_NEAR( -4.6487198600333182713, LagrangianPoissonDistribution( log( 10 ), 0.7 ).lnPdf( 50 ), LN_PDF_EPS );
}

class LagrangianPoissonTest: public ::testing::TestWithParam< std::pair<double, double> > {
protected:
    gsl_rng* rndGen;

    LagrangianPoissonTest()
    : rndGen( gsl_rng_alloc( gsl_rng_default ) )
    {
    }

    virtual ~LagrangianPoissonTest()
    {
        gsl_rng_free( rndGen );
    }
};

TEST_P( LagrangianPoissonTest, RandomSample )
{
    LagrangianPoissonDistribution lp( GetParam().first, GetParam().second );
    GenericDiscreteDistribution gdd( lp, 1E-5 );
    size_t sum = 0;
    size_t n = 100000;
    for ( size_t i = 0; i < n; i++ ) {
        sum += gdd.random( rndGen );
    }
    EXPECT_NEAR( lp.mean(), (double)sum / n, 5E-3 * lp.mean() );
}

INSTANTIATE_TEST_CASE_P( PoissonRandomSample,
                         LagrangianPoissonTest,
                        ::testing::Values( std::make_pair( 1.0, 0.0 ),
                                           std::make_pair( 2.0, 0.0 ),
                                           std::make_pair( 4.0, 0.0 ),
                                           std::make_pair( 8.0, 0.0 ),
                                           std::make_pair( 1.0, 0.1 ),
                                           std::make_pair( 2.0, 0.1 ),
                                           std::make_pair( 4.0, 0.1 ),
                                           std::make_pair( 8.0, 0.1 ),
                                           std::make_pair( 1.0, 0.3 ),
                                           std::make_pair( 2.0, 0.3 ),
                                           std::make_pair( 4.0, 0.3 ),
                                           std::make_pair( 8.0, 0.3 )
                        ) );
