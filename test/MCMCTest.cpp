#include <gsl/gsl_sf_gamma.h>

#include <math/Distributions.h>

#include <mcmc/MetropolisHastingsStep.h>
#include <mcmc/TransitionDistributions.h>

#include "TestCommon.h"

/*************************************************/
TEST_F2( GslRngTestF, MCMC_MH_steps, DISABLED_weibull_mean_sampling )
{
    double  sum2 = 0.0;
    double  sum = 0.0;
    double  val = 1.0;
    
    const   int     totalCycles = 2000000;
    const   int     samplesPerSmallSum = 1000;
    const   int     samplePeriod = 50;
    const   int     burnInCycles = 100;
    int     nsum = 0;
    const   double  wa = 5;
    const   double  wb = 2;

    for ( int i = 0; i < totalCycles; i++ )
    {
        val = MetropolisHastringsPosteriorSample<double>( 
            rndGen, 
            PositiveGaussianTransitionDistribution( 15.0 ),
            //GammaTransitionDistribution( 0.001 ),
            //ExponentialTransitionDistribution( /*10*/ ),
            //UniformTransitionDistribution( 0, 5000 ),
            WeibullDistribution( wa, wb ),
            DegeneratedDistribution<double>(),
            val
        );
        if ( i > burnInCycles ) {
            if ( i % samplePeriod ) {
                nsum++;
                sum += val;
                if ( nsum % samplesPerSmallSum == 0 ) {
                    sum2 += sum;
                    sum = 0.0;
                }
            }
        }
    }
    sum2 += sum;
    sum2 /= nsum;
    EXPECT_NEAR( sum2, wa * gsl_sf_gamma( 1 + 1.0 / wb ), 0.01 );
}
