# BIMAP method testing
# 
# Author: Alexey Stukalov
###############################################################################

bimap_scripts_path <- ".."
source( paste( bimap_scripts_path, "BIMAP.R", sep = .Platform$file.sep ) )
source( paste( bimap_scripts_path, "BIMAP_plot.R", sep = .Platform$file.sep ) )
source( paste( bimap_scripts_path, "BIMAP_gen.R", sep = .Platform$file.sep ) )

require( compoisson )

compare_comp <- function( vals, shape )
{
    plot( vals, lapply( vals, function(x) ComPoissonLnZeta(log(x),shape) ), 
        type="l", col="blue" )
   # lines( vals, lapply( vals, function(x) com.compute.log.z(x^shape,shape)) )
}

compare_comp2 <- function( vals, shape )
{
    plot( vals, lapply( vals, function(x) ( ComPoissonLnZeta(log(x),shape) - com.compute.log.z(x^shape,shape) ) ), 
        type="l", col="blue" )
}
lapply( 200:350, function(x) ComPoissonLnZeta(log(x),2.9) )

compare_comp( 1:750, 1.0 )
compare_comp2( 1:790, 1.0 )

#########################################
# test0
#########################################

test0 <- list(
    proteins = data.frame(
        label = c('A'),
        sequence.length = as.integer( c(10) )
    ),
    samples = data.frame( 
        label = c('1'),
        bait_label = c('A')
    ),
    msruns = data.frame(
        label = c( '1-1', '1-2', '1-3', '1-4', '1-5', '1-6' ),
        sample_label = c( '1', '1', '1', '1', '1', '1' )
    ),
    objects.clusters = data.frame(
        objects.cluster.serial = c( 1 ),
        object = c( 'A' ),
        object.multiple = as.integer( c( 1 ) )
    ),
    states.clusters = data.frame(
        states.cluster.serial = c( 1 ),
        state = c( '1' )
    ),
    crossClusters = data.frame(
        objects.cluster.serial = c( 1 ),
        states.cluster.serial = c( 1 )
    ),
    stateSignals = data.frame(
        objects.cluster.serial = c( 1 ),
        state = c('1' ),
        state.signal = c( 0 )
    )
)

test0$signals <- MSSignalsMatrix(
    test0$objects.clusters,
    test0$states.clusters,
    test0$crossClusters,
    test0$stateSignals
)

test0$measurements.matrix <- MSDataMatrixRandom( 
    test0$proteins, test0$msruns, test0$signals, 
    seq.length.factor = 1.0,
    signal.shape = 3.0, noise.shape = 1.0 )
test0$measurements.frame <- MSMatrixToFrame( test0$measurements.matrix )

test0$walk <- BIMAPMassSpecData(
    test0$proteins, test0$samples, test0$msruns, test0$measurements.frame,
    eesampler.burnInIterations = 5000,
    eesampler.turbinesCount = 3,
    samplingIterations = 30000,
    ini.objects.partition = test0$objects.clusters,
    ini.states.partition = test0$states.clusters,
    ini.crossClusters = test0$crossClusters,
    object.membership.updatePeriod = 0,
    state.membership.updatePeriod = 0,
    cell.membership.updatePeriod = 0,
    prior.sequence.length.factor = 0.5,
    prior.objects.clustering.concentration = 0.1,
    prior.objects.clustering.discount = 0.0,
    prior.states.clustering.concentration = 0.1,
    prior.states.clustering.discount = 0.0,
    prior.baseline.sigma = 10.0,
    prior.noise.gamma.shape = 3.0,
    prior.noise.gamma.scale = 3.0,
    prior.signal.shape = 1.0,
    prior.noise.shape = 1.0,
    object.signal.temperature = 100.0,
    hyperprior.clusterStateSignal.variance.shape = 1, 
    hyperprior.clusterStateSignal.variance.scale = 0.1
)

test0$walk.clusterings.frame <- BIMAP.walk.data.frame( test0$walk ) 
test0$walk.priors.frame <- BIMAP.walk.priors.data.frame( test0$walk ) 

test0$walk.signals.frame <- BIMAP.walk.signals.data.frame( test0$walk ) 

plot.cell.signal( test0$walk.signals.frame, 'A', '1' )

plot.signal_noise_sigma( test0$walk.clusterings.frame, test0$walk.priors.frame )

plot( test0$walk.clusterings.frame$serial )
plot( test0$walk.cluster.frame$iteration, test0$walk.cluster.frame$cluster )

plot.bimap( find_clustering( test0$walk, 0 ), 
    test0$proteins, test0$samples, test0$msruns 
)

plot_cluster_samples( test0$walk.cluster.frame, 2 )

0.5 * log( subset( test1$proteins, label == 'A' )$sequence_length )
mean( subset( test1$walk.cluster.frame, cluster == 2 & iteration > 10000 )$signal )

#########################################
# test1
#########################################

test1 <- list(
    proteins = data.frame(
        label = c('A','B','C'),
        sequence.length = as.integer( c(50,50,50) )
    ),
    samples = data.frame( 
        label = c('1'),
        bait_label = c('A')
    ),
    msruns = data.frame(
        label = c( '1-1', '1-2', '1-3', '1-4', '1-5', '1-6' ),
        sample_label = c( '1', '1', '1', '1', '1', '1' )
    ),
    objects.clusters = data.frame(
        objects.cluster.serial = c( 1, 1, 2 ),
        object = c( 'A', 'B', 'C' ),
        object.multiple = as.integer( c( 1, 1, 1 ) )
    ),
    states.clusters = data.frame(
        states.cluster.serial = c( 1 ),
        state = c( '1' )
    ),
    crossClusters = data.frame(
        objects.cluster.serial = c( 1, 2 ),
        states.cluster.serial = c( 1, 1 )
    ),
    stateSignals = data.frame(
        objects.cluster.serial = c( 1, 2 ),
        state = c('1', '1' ),
        state.signal = c( 2, 8 )
    )
)

test1$signals <- MSSignalsMatrix(
    test1$objects.clusters,
    test1$states.clusters,
    test1$crossClusters,
    test1$stateSignals
)

test1$measurements.matrix <- MSDataMatrixRandom( 
    test1$proteins, test1$msruns, test1$signals, 
    seq.length.factor = 0.5,
    signal.shape = 2, noise.shape = 1.0 )
test1$measurements.frame <- MSMatrixToFrame( test1$measurements.matrix )

test1$walk <- BIMAPMassSpecData(
    test1$proteins, test1$samples, test1$msruns, test1$measurements.frame,
    burnInIterations = 0, samplingIterations = 8000,
    ini.objects.partition = test1$objects.clusters,
    ini.states.partition = test1$states.clusters,
    ini.crossClusters = test1$crossClusters,
    object.membership.updatePeriod = 25,
    state.membership.updatePeriod = 0,
    cell.membership.updatePeriod = 0,
    prior.sequence.length.factor = 0.5,
    prior.objects.clustering.concentration = 0.1,
    prior.objects.clustering.discount = 0.01,
    prior.states.clustering.concentration = 0.1,
    prior.states.clustering.discount = 0.01,
    prior.baseline.sigma = 10,
    prior.signal.shape = 2.0,
    prior.noise.shape = 1.0,
    hyperprior.clusterStateSignal.variance.shape = 1, 
    hyperprior.clusterStateSignal.variance.scale = 0.1
)

test1$walk.clusterings.frame <- BIMAP.walk.data.frame( test1$walk ) 
test1$walk.priors.frame <- BIMAP.walk.priors.data.frame( test1$walk ) 
#test1$walk.proteins.clusters.frame <- BIMAP.walk.proteins.clusters.data.frame( test1$walk ) 
#test1$walk.samples.clusters.frame <- BIMAP.walk.samples.clusters.data.frame( test1$walk ) 

test1$walk.signals.frame <- BIMAP.walk.signals.data.frame( test1$walk ) 

plot.cell.signal( test1$walk.signals.frame, 'A', '1' )
plot.cell.signal( test1$walk.signals.frame, 'C', '1' )

plot( subset( test1$walk.signals.frame, object == 'B' & state == '1' )$objects.cluster.serial )
plot( subset( test1$walk.signals.frame, object == 'C' & state == '1' )$objects.cluster.serial )

plot.signal_noise_sigma( test1$walk.clusterings.frame, test1$walk.priors.frame )

plot( test1$walk.clusterings.frame$serial )
plot( test1$walk.cluster.frame$iteration, test1$walk.cluster.frame$cluster )

plot.bimap( find_clustering( test1$walk, 1 ), 
    test1$proteins, test1$samples, test1$msruns, test1$walk.signals.frame 
)

plot_cluster_samples( test1$walk.cluster.frame, 2 )

0.5 * log( subset( test1$proteins, label == 'A' )$sequence_length )
mean( subset( test1$walk.cluster.frame, cluster == 2 & iteration > 10000 )$signal )

###################################################
# test2
###################################################
test2 <- list(
    proteins = data.frame(
        label = c('A','B','C','D'),
        sequence.length = as.integer( c(50,50,50,50) )
    ),
    samples = data.frame( 
        label = c('1','2','3','4'),
        bait_label = c('A','B','C','D')
    ),
    msruns = data.frame(
        label = c( '1-1', '1-2', '1-3', '1-4',
                   '2-1', '2-2', '2-3', '2-4',
                   '3-1', '3-2', '3-3', '3-4',
                   '4-1', '4-2', '4-3', '4-4' ),
        sample_label = c( '1', '1', '1', '1',
                          '2', '2', '2', '2',
                          '3', '3', '3', '3',
                          '4', '4', '4', '4'
                          )
                      ),
    objects.clusters = data.frame(
        objects.cluster.serial = c( 1, 1, 2, 2 ),
        object = c( 'A', 'B', 'C', 'D' ),
        object.multiple = as.integer( c( 1, 1, 1, 1 ) )
    ),
    states.clusters = data.frame(
        states.cluster.serial = c( 1, 1, 2, 2 ),
        state = c( '1', '2', '3', '4' )
    ),
    crossClusters = data.frame(
        objects.cluster.serial = c( 1, 1, 2 ),
        states.cluster.serial = c( 1, 2, 2 )
    ),
    stateSignals = data.frame(
        objects.cluster.serial = c( 1, 1, 2, 2, 2, 2 ),
        state = c( '1', '2', '1', '2', '3', '4' ),
        state.signal = c( 2, 8, 1, 1, 2, 8 )
    )
)

test2$signals <- MSSignalsMatrix(
    test2$objects.clusters,
    test2$states.clusters,
    test2$crossClusters,
    test2$stateSignals
)

test2$measurements.matrix <- MSDataMatrixRandom( 
    test2$proteins, test2$msruns, test2$signals, 
    seq.length.factor = 0.5,
    signal.shape = 2.0,
    noise.shape = 1.0, noise.peak = -3 )
test2$measurements.frame <- MSMatrixToFrame( test2$measurements.matrix )

test2$walk <- BIMAP.walk.eval(
    test2$proteins, test2$samples, test2$msruns, test2$measurements.frame,
    file = '/home/astukalov/test2walk.xml',
    eesampler.burnInIterations = 5000,
    eesampler.turbinesCount = 5,
    eesampler.ringsCount = 30,
    eesampler.generateRate = 0.1, 
    samplingIterations = 25000,
    sampleRate.signal = 0.3,
    sampleRate.object.membership = 0.01,
    sampleRate.objects.splitMerge = 0.06,
    sampleRate.state.membership = 0.01,
    sampleRate.states.splitMerge = 0.06,
    signal.sequence.length.factor = 0.3,
    signal.temperature = 10.0,
    prior.objects.clustering.concentration = 0.02,
    prior.objects.clustering.discount = 0.01,
    prior.states.clustering.concentration = 0.02,
    prior.states.clustering.discount = 0.01,
    prior.noise.gamma.shape = 3.0,
    prior.noise.gamma.scale = 0.3,
    prior.signal.shape = 1.0,
    prior.noise.shape = 1.0,
    prior.crossCluster.enabled = 0.1,
    hyperprior.signal.variance.scale = 0.1,
    hyperprior.signal.variance.shape = 10,
    hyperprior.baseline.scale = 2
)

test2$walk.clusterings.frame <- BIMAP.walk.data.frame( test2$walk ) 
test2$walk.priors.frame <- BIMAP.walk.priors.data.frame( test2$walk ) 
#test2$walk.proteins.clusters.frame <- BIMAP.walk.proteins.clusters.data.frame( test2$walk ) 
#test1$walk.samples.clusters.frame <- BIMAP.walk.samples.clusters.data.frame( test1$walk ) 

test2$walk.signals.frame <- BIMAP.walk.signals.data.frame( test2$walk ) 

plot.cell.signal( test2$walk.signals.frame, 'A', '1' )
plot.cell.signal( test2$walk.signals.frame, 'B', '2' )
plot.cell.signal( test2$walk.signals.frame, 'C', '3' )
plot.cell.signal( test2$walk.signals.frame, 'D', '4' )

plot( subset( test1$walk.signals.frame, object == 'B' & state == '1' )$objects.cluster.serial )
plot( subset( test1$walk.signals.frame, object == 'C' & state == '1' )$objects.cluster.serial )

plot.signal_noise_sigma( test2$walk.clusterings.frame, test2$walk.priors.frame )

plot( as.factor( test2$walk.signals.frame$object.multiple ) )

par( oma = c(0,0,0,0), mar = c(2,2,0,0) )
plot( test2$walk.clusterings.frame$serial )
plot( test2$walk.clusterings.frame$objects.partition.serial )
plot( test2$walk.clusterings.frame$states.partition.serial )
plot( test2$walk.cluster.frame$iteration, test2$walk.cluster.frame$cluster )

par( oma = c(0,0,0,0), mar = c(0,2,2,0) )
plot.bimap( BIMAP.walk.lastClustering( test2$walk, 0 ) 
    , test2$proteins, test2$samples, test2$msruns
    , test2$walk.signals.frame
)

plot_cluster_samples( test1$walk.cluster.frame, 2 )

0.5 * log( subset( test1$proteins, label == 'A' )$sequence_length )
mean( subset( test1$walk.cluster.frame, cluster == 2 & iteration > 10000 )$signal )

plot_ms_data( test.measurements.matrix )

plot( test.walk.clusterings.frame$serial )

plot.bimap( find_clustering( test.walk, 0 ), 
    test.prot_info, test.sample_info, test.msrun_info 
)

plot_cluster_samples( test.walk.cluster.frame, 2 )

proteins.data.frame( test.clusters[[0]] )


#########################################
# test3
#########################################

test3 <- list(
    proteins = data.frame(
        label = c( 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J' ),
        sequence.length = as.integer( c( 10, 20, 15, 10, 12, 12, 17, 9, 10, 11 ) )
    ),
    samples = data.frame( 
        label = c( '1', '2', '3', '4', '5' ),
        bait_label = c( 'A', 'B', 'C', 'D', 'E' )
    ),
    msruns = data.frame(
        label = c( '1-1', '2-1', '2-2', '2-3', '3-1', '3-2', '4-1', '4-2', '5-1', '5-2', '5-3', '5-4' ),
        sample_label = c( '1', '2', '2', '2', '3', '3', '4', '4', '5', '5', '5', '5' )
    ),
    objects.clusters = data.frame(
        objects.cluster.serial = c( 1, 1, 2, 2, 2, 3, 4, 4, 5, 5 ),
        object = c( 'A', 'F', 'B', 'C', 'G', 'D', 'E', 'H', 'I', 'J' ),
        object.multiple = as.integer( c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ) )
    ),
    states.clusters = data.frame(
        states.cluster.serial = c( 1, 1, 2, 2, 3 ),
        state = c( '1', '2', '3', '4', '5' )
    ),
    crossClusters = data.frame(
        objects.cluster.serial = c( 1, 2, 2 ),
        states.cluster.serial = c( 1, 1, 2 )
    ),
    stateSignals = data.frame(
        objects.cluster.serial = c( 1, 2, 2 ),
        state = c( '1', '2', '1', '2', '3', '4' ),
        state.signal = c( 2, 8, 4, 2, 1, 1 )
    )
)

test3$signals <- MSSignalsMatrix(
    test3$objects.clusters,
    test3$states.clusters,
    test3$crossClusters,
    test3$stateSignals
)

test3$measurements.matrix <- MSDataMatrixRandom( 
    test3$proteins, test3$msruns, test3$signals, 
    seq.length.factor = 0.5,
    signal.shape = 2, noise.shape = 1.0 )
test3$measurements.frame <- MSMatrixToFrame( test3$measurements.matrix )

test3$walk <- BIMAPMassSpecData(
    test3$proteins, test3$samples, test3$msruns, test3$measurements.frame,
    clustersCount = 4, burnInIterations = 0, samplingIterations = 80000,
    hyperprior.state.membership = 0.1,
    hyperprior.object.membership = 0.1,
    hyperprior.clusterStateSignal.variance.scale = 0.1,
    hyperprior.noise.sigma = 20.0,
    hyperprior.baseline.sigma = 20.0
)

test3$walk.clusterings.frame <- BIMAP.walk.data.frame( test3$walk ) 
test3$walk.priors.frame <- BIMAP.walk.priors.data.frame( test3$walk ) 
test3$walk.cluster.frame <- BIMAP.walk.cluster.data.frame( test3$walk ) 

plot.signal_noise_sigma( test3$walk.clusterings.frame, test3$walk.priors.frame )


plot( test3$walk.clusterings.frame$serial )
plot( test3$walk.cluster.frame$iteration, test3$walk.cluster.frame$cluster )

cluster.frame <- unique( subset( test3$walk.cluster.frame, select = c( iteration, cluster ) ) )

plot( as.factor( cluster.frame$cluster ) )

find_clustering( test3$walk, 16 )
plot_biclustering( find_clustering( test3$walk, 0 ), 
    test3$proteins, test3$samples, test3$msruns 
)

plot_cluster_samples( test.walk.cluster.frame, 17 )
