# BIMAP R interface basic definitions.
# 
# Author: Alexey Stukalov (astukalov@cemm.oeaw.ac.at)
###############################################################################

require( 'plyr' )
require( 'Rcpp' )

# load RBIMAP library
# Note: RBIMAP.libpath should be set before executing this script, e.g.:
#     RBIMAP.libpath <- file.path( bimap_scripts_path, "build/release/src/R" )
#dyn.unload( file.path( RBIMAP.libpath, paste("libRBIMAP", .Platform$dynlib.ext, sep="")) ) 
dyn.load( file.path( RBIMAP.libpath, paste("libRBIMAP", .Platform$dynlib.ext, sep="")), type = "Call" ) 

#' Imports MS data into format compatible to submission to BIMAP 
#' @param ms_data 
#' @param protein_info 
#' @param msrun.multipliers 
#' @param sample_column 
#' @param msrun_column 
#' @param bait_column 
#' @param prey_column 
#' @param sc_column 
#' @param pc_column 
#' @param protein_ac_column 
#' @param protein_seqlength_column 
#' @returnType 
#' @return 
#' @author astukalov
#' @export
BIMAP.msdata.import <- function( ms_data, protein_info, msrun.multipliers = NULL,
    sample_column = 'sample', msrun_column = 'msrun',
    bait_column = 'bait_ac', prey_column = 'prey_ac',
    sc_column = 'sc', pc_column = 'pc',
    sample_extra_columns = c(),
    protein_ac_column = 'primaryac', protein_seqlength_column = 'seqlength',
    protein_extra_columns = c()
){
    ms_data_noglob <- ms_data[ ms_data[, msrun_column ] != '#glob#', ]
    measurements.df <- data.frame(
            msrun = as.character( ms_data_noglob[, msrun_column ] ),
            prey_ac = as.character( ms_data_noglob[, prey_column ] ),
            sc = as.integer( ms_data_noglob[, sc_column ] ),
            stringsAsFactors = FALSE
        )
    if ( !is.null( pc_column ) ) {
        measurements.df$pc <- as.integer( ms_data_noglob[, pc_column ] )
    }
    # remove measurements duplicates -- only one per msrun-prey pair is allowed
    mes_pkeys <- list( measurements.df$msrun, measurements.df$prey_ac )
    measurements.df.x <- data.frame(
        msrun = as.vector( tapply( measurements.df$msrun, mes_pkeys, function( x ) x[1] ) ),
        prey_ac= as.vector( tapply( measurements.df$prey_ac, mes_pkeys, function( x ) x[1] ) ),
        sc = as.vector( tapply( measurements.df$sc, mes_pkeys, max ) ),
        stringsAsFactors = FALSE
    )
    if ( !is.null( pc_column ) ) {
        measurements.df.x$pc <- as.vector( tapply( measurements.df$pc, mes_pkeys, max ) )
    }
    # remove non-measurements
    measurements.df <- subset( measurements.df.x, sc > 0 )
    samples.colnames <- unique( c( sample_column, bait_column,
                                   intersect( colnames(ms_data_noglob),
                                              sample_extra_columns ) ) )
    if ( sample_column == bait_column ) samples.colnames <- c( sample_column, samples.colnames ) # handle case sample = bait
    print( samples.colnames )
    print( head( ms_data_noglob ) )
    samples.df <- unique( subset( ms_data_noglob, select = samples.colnames ) )
    samples.colnames[1:2] <- c( 'sample', 'bait_ac' )
    colnames( samples.df ) <- samples.colnames 
    rownames( samples.df ) <- samples.df$sample
    samples.df <- samples.df[ order( samples.df$bait_ac, samples.df$sample ), ]
    msruns.df <- as.data.frame( unique( subset( ms_data_noglob, select = unique( c( msrun_column, sample_column ) ) ) ),
                                stringsAsFactors = FALSE )
    msruns.df <- data.frame( 
            msrun = as.character( msruns.df[, msrun_column ] ),
            sample = as.character( msruns.df[, sample_column ] ),
            stringsAsFactors = FALSE )
    msruns.df$bait_ac <- samples.df[ msruns.df$sample, 'bait_ac' ]
    msruns.df <- msruns.df[ order( msruns.df$bait_ac, msruns.df$sample, msruns.df$msrun ), ]
    msruns.df$bait_ac <- NULL
    rownames( msruns.df ) <- msruns.df$msrun
    if ( !is.null(msrun.multipliers) ) {
        mean.mult <- mean( msrun.multipliers[ msruns.df$msrun ], na.rm = TRUE )
        msruns.df$multiplier <- msrun.multipliers[ msruns.df$msrun ] / mean.mult
        msruns.df[ is.na(msruns.df$multiplier), 'multiplier' ] <- 1
    } else {
        msruns.df$multiplier <- 1
    }
    # prepare proteins frame
    protein_extra_columns <- setdiff( unique( protein_extra_columns ), protein_seqlength_column )
    proteins.df <- protein_info[ protein_info[ , protein_ac_column ] %in% 
                union( as.character( ms_data_noglob[ , prey_column ] ), 
                    as.character( ms_data_noglob[ , bait_column ] ) ), 
                                 c( protein_ac_column, protein_seqlength_column,
                                    protein_extra_columns ) ]
    colnames( proteins.df ) <- c( 'protein_ac', 'seqlength', protein_extra_columns )
    # group all protein records referring to the same protein AC
    # (this is important if proteins.df contains multiple isosofrms)
    protein_ids <- plyr:::id_var( proteins.df[,c('protein_ac')], drop = TRUE )
    protein_indices <- plyr:::split_indices( seq_len(nrow(proteins.df)), protein_ids, n = attr( protein_ids, "n" ) )
    protein_first_index <- vapply( protein_indices, function(indices) indices[[1]], integer(1) )
    # sequence length is the maximum in AC group
    seqlength <- vapply( protein_indices,
                         function(indices) max( proteins.df[indices,'seqlength'] ),
                         integer(1) )
    # leave one (first) record per protein AC
    proteins.df <- proteins.df[ protein_first_index, ]
    proteins.df$seqlength <- seqlength
    #print( as.character( ms_data_noglob$prey_ac ) )
    rownames( proteins.df ) <- proteins.df$protein_ac 

    exp_design.df <- merge( samples.df, msruns.df )[,c('bait_ac','sample','msrun','multiplier')]
    rownames( exp_design.df ) <- exp_design.df$msrun
    return ( list(  measurements = measurements.df,
                    samples = samples.df, 
                    msruns = msruns.df, 
                    exp_design = exp_design.df,
                    proteins = proteins.df ) )
}

#' Calculated MS runs multipliers based on total spectral counts.
#' MS run multiplier is a sum of protein spectra divided
#' by a mean of such sum across all technical replicates of the same sample.
#' @param ms_data MS data dataframe with SC
#' @param sample_column column for ID of biological sample
#' @param msrun_column column for ID of technical replicate
#' @param prey_column column for protein AC code.
#'        If MS data contains multiple entries for the protein, the maximal SC is used.
#' @param sc_column column for spectral counts
#' @returnType vector
#' @return multipliers vector with the names being MS run IDs
#' @author Alexey Stukalov
#' @export
BIMAP.msrun.multipliers <- function( ms_data, sample_column = 'sample', msrun_column = 'msrun',
                                              prey_column = 'prey_ac', sc_column = 'sc' )
{
    res <- ddply( ms_data, sample_column, function( sample_data ) {
        res <- ddply( sample_data, msrun_column, function( msrun_data ) {
                sc_data <- tapply( msrun_data[,sc_column], msrun_data[,prey_column], max )
                data.frame( sc.sum = sum( sc_data ),
                            stringsAsFactors = FALSE )
        } )
        res$multiplier <- res$sc.sum / mean( res$sc.sum )
        return ( res )
    } )
    t <- res$multiplier
    names( t ) <- res$msrun
    return ( t )
}

#' Saves AP-MS data into single .xml file that could be read by BIMAP-sampler --input_file
#' @param filename name of output file
#' @author Alexey Stukalov
BIMAP.msdata.load <- function(
    data_path = NULL, filename = NULL,
    format = c( 'OPAData', 'CSV' )
){
    if ( format == 'OPAData'){
        return ( .Call( "OPADataLoad", filename ) )
    } else {
        message( "Reading BI-MAP data from ", data_path, "..." )
        bimap.data <- list()
        bimap.data$measurements <- read.table(
            file = file.path( data_path, 'measurements.txt' ),
            header = TRUE, sep = '\t',
            stringsAsFactors = FALSE )

        bimap.data$proteins <- read.table(
            file = file.path( data_path, 'proteins.txt' ),
            header = TRUE, sep = '\t',
            stringsAsFactors = FALSE )
        rownames( bimap.data$proteins ) <- bimap.data$proteins$protein_ac

        bimap.data$exp_design <- read.table(
            file = file.path( data_path, 'exp_design.txt' ),
            header = TRUE, sep = '\t',
            stringsAsFactors = FALSE )

        bimap.data$samples <- unique( bimap.data$exp_design[,c('bait_ac','sample')] )
        rownames( bimap.data$samples ) <- bimap.data$samples$sample

        bimap.data$msruns <- unique( bimap.data$exp_design[,c('sample','msrun', 'multiplier')] )
        rownames( bimap.data$msruns ) <- bimap.data$msruns$msrun
    }
    return ( bimap.data )
}

#' Saves AP-MS data into single .xml file that could be read by BIMAP-sampler --input_file
#' @param filename name of output file
#' @author Alexey Stukalov
BIMAP.msdata.save <- function( bimap.data, 
    data_path = NULL, filename = NULL,
    format = 'OPA',
    map.baits.to.objects = TRUE
){
    if ( format == 'OPA' ){
        return ( .Call( "OPADataSave", filename,
                bimap.data$proteins, bimap.data$samples, bimap.data$msruns,
                bimap.data$measurements,
                list( map.baits.to.objects = map.baits.to.objects ) ) )
    } else {
        message( "Writing BI-MAP data to ", data_path, "..." )
        write.table( bimap.data$measurements[,intersect(colnames(bimap.data$measurements),
                    c('msrun','prey_ac','sc','pc'))],
            file = file.path( data_path, 'measurements.txt' ),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t' )
        write.table( bimap.data$exp_design[,c('bait_ac','sample','msrun','multiplier')],
            file = file.path( data_path, 'exp_design.txt' ),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t' )
        write.table( rbind( bimap.data$proteins[,c('protein_ac','seqlength')],
                data.frame( protein_ac = 'nobait', seqlength = 1 ) ),
            file = file.path( data_path, 'proteins.txt' ),
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t' )
    }
}

#' Appends extra information to BI-MAP input data structure
#' @param bimap.data 
#' @param protein.info (optional) extra protein information dataframe 
#' @param sample.info (optional) extra samples information dataframe
#' @param msrun.info (optional) extra MS runs information dataframe
#' @author Alexey Stukalov
BIMAP.msdata.extra_info <- function( bimap.data,
                                     proteins.info = NULL,
                                     samples.info = NULL,
                                     msruns.info = NULL )
{
    if ( !is.null( proteins.info ) ) {
        bimap.data$proteins <- cbind( bimap.data$proteins, proteins.info[ rownames( bimap.data$proteins ),
                                      setdiff( colnames( proteins.info ), colnames( bimap.data$proteins ) ) ])
    }
    if ( !is.null( samples.info ) ) {
        bimap.data$samples <- cbind( bimap.data$samples, samples.info[ rownames( bimap.data$samples ),
                                      setdiff( colnames( samples.info ), colnames( bimap.data$samples ) ) ])
    }
    if ( !is.null( msruns.info ) ) {
        bimap.data$msruns <- cbind( bimap.data$msruns, msruns.info[ rownames( bimap.data$msruns ),
                                      setdiff( colnames( msruns.info ), colnames( bimap.data$msruns ) ) ])
    }
    return ( bimap.data )
}

#' Runs BIMAP method within R session
#' @param protein_info 
#' @param sample_info 
#' @param msrun_info 
#' @param measurements 
#' @param walk.samples 
#' @param walk.create.RObject 
#' @param walk.file 
#' @return 
#' @author Alexey Stukalov
BIMAP.mcmcwalk.eval <- function(
    protein_info,
    sample_info,
    msrun_info,
    measurements,
    map.baits.to.objects = TRUE,
    walk.samples = 1000,
    walk.create.RObject = TRUE,
    walk.file = NULL,
    eesampler.burnInIterations = 1000,
    eesampler.ladderAdjustPeriod = 200,
    eesampler.turbinesCount = 2,
    eesampler.levelsCount = 2,
    eesampler.energyTolerance = 0.01,
    eesampler.jumpRate = 0.01,
    eesampler.generateRate = 0.005,
    eesampler.temperatureMultiplier = 1.8,
    eesampler.detentionIterations = 29,
    sampleRate.object.membership = 0.4,
    sampleRate.objects.splitMerge = 0.4,
    sampleRate.probe.membership = 0.3,
    sampleRate.probes.splitMerge = 0.3,
    sampleRate.block.flip = 0.05,
    sampleRate.object.multiple = 0.125,
    sampleRate.signal = 0.1,
    samplePeriod.priors = 197,
    samplePeriod.chessboardBiclustering = 23,
    block.resamples = 12,
    ini.objects.partition = NULL,
    ini.probes.partition = NULL,
    ini.blocks = NULL,
    signal.sequence.length.factor = 0.5,
    signal.shape = 0.0,
    prior.objects.clustering.concentration = 0.1,
    prior.objects.clustering.discount = 0.0,
    prior.probes.clustering.concentration = 0.1,
    prior.probes.clustering.discount = 0.0,
    hyperprior.baseline = 0.0,
    hyperprior.baseline.scale = 1.0,
    hyperprior.signal.variance.shape = 1.0,
    hyperprior.signal.variance.scale = 0.1,
    prior.true_misses = 1000,
    prior.false_hits = 1,
    prior.block.enabled = 0.1,
    precomp.object.freq.threshold = 0.6,
    precomp.probe.freq.threshold = 0.6,
    objects.components.threshold = 0.1,
    probes.components.threshold = 0.1,
    objects.clot.threshold = 0.9,
    stable.objects.clusters.threshold = NA,
    stable.probes.clusters.threshold = NA,
    log.particles.file = NULL,
    log.eeJumps.file = NULL
){
    params <- list(
        map.baits.to.objects = map.baits.to.objects,
        file = file,
        walk.create.RObject = walk.create.RObject,
        walk.samples = walk.samples,
        walk.file = walk.file,
        eesampler.burnInIterations = eesampler.burnInIterations, 
        eesampler.ladderAdjustPeriod = eesampler.ladderAdjustPeriod, 
        eesampler.jumpRate = eesampler.jumpRate, 
        eesampler.generateRate = eesampler.generateRate, 
        eesampler.levelsCount = eesampler.levelsCount,
        eesampler.turbinesCount = eesampler.turbinesCount,
        eesampler.energyTolerance = eesampler.energyTolerance,
        eesampler.temperatureMultiplier = eesampler.temperatureMultiplier,
        eesampler.detentionIterations = eesampler.detentionIterations,
        sampleRate.object.membership = sampleRate.object.membership,
        sampleRate.objects.splitMerge = sampleRate.objects.splitMerge,
        sampleRate.probe.membership = sampleRate.probe.membership,
        sampleRate.probes.splitMerge = sampleRate.probes.splitMerge,
        sampleRate.block.flip = sampleRate.block.flip,
        sampleRate.object.multiple = sampleRate.object.multiple,
        sampleRate.signal = sampleRate.signal,
        samplePeriod.priors = samplePeriod.priors,
        samplePeriod.chessboardBiclustering = samplePeriod.chessboardBiclustering,
        block.resamples = block.resamples,
        ini.objects.partition = ini.objects.partition,
        ini.probes.partition = ini.probes.partition,
        ini.blocks = ini.blocks,
        signal.sequence.length.factor = signal.sequence.length.factor,
        signal.shape = signal.shape,
        prior.objects.clustering.concentration = prior.objects.clustering.concentration,
        prior.objects.clustering.discount = prior.objects.clustering.discount,
        prior.probes.clustering.concentration = prior.probes.clustering.concentration,
        prior.probes.clustering.discount = prior.probes.clustering.discount,
        hyperprior.baseline = hyperprior.baseline,
        hyperprior.baseline.scale = hyperprior.baseline.scale,
        hyperprior.signal.variance.shape = hyperprior.signal.variance.shape,
        hyperprior.signal.variance.scale = hyperprior.signal.variance.scale,
        prior.true_misses = prior.true_misses,
        prior.false_hits = prior.false_hits,
        prior.block.enabled = prior.block.enabled,
        precomp.object.freq.threshold = precomp.object.freq.threshold,
        precomp.probe.freq.threshold = precomp.probe.freq.threshold,
        objects.components.threshold = objects.components.threshold,
        probes.components.threshold = probes.components.threshold,
        objects.clot.threshold = objects.clot.threshold,
        stable.objects.clusters.threshold = stable.objects.clusters.threshold,
        stable.probes.clusters.threshold = stable.probes.clusters.threshold,
        log.particles.file = log.particles.file,
        log.eeJumps.file = log.eeJumps.file
    )
    res <- .Call( "BIMAPWalkEval", 
                  protein_info, sample_info, msrun_info, 
                  subset( measurements, select = intersect( c( 'prey_ac', 'msrun', 'sc', 'pc' ), colnames(measurements) ), msrun != '#glob#' ), 
                  params )
    return ( res )
}

# MPI (parallelized) version of BIMAP sampler 
BIMAP.mcmcwalk.eval_MPI <- function(
    protein_info,
    sample_info,
    msrun_info,
    measurements,
    mpirun.options = c( "-n", 2 ),
    config.file = '/home/astukalov/projects/data/egfr/BIMAP_config.ini',
    program.path = RBIMAP.libpath,
    walk.create.RObject = TRUE,
    walk.file = NULL
){
    tempOPAfilename <- paste( tempfile( "osadata_" ), ".xml.gz", sep='' )
    OPAData.save( tempOPAfilename, protein_info, sample_info, msrun_info, measurements )
    sampler_options <- c( '--input_file', tempOPAfilename )
    if ( !is.null(config.file) ) {
        sampler_options <- c( sampler_options, "--config_file", config.file )
    }
    tempWalk.file = NULL
    if ( is.null(walk.file) ) {
        # would generate temporarily walk file
        tempWalk.file <- paste( tempfile( "BIMAPwalk_" ), ".xml.gz", sep='' )
        walk.file <- tempWalk.file
    }
    sampler_options <- c( sampler_options, "--output_file", walk.file )
    cmdline <- paste( "mpiexec", paste( mpirun.options, collapse=' ' ),
                      file.path( program.path, "MPIBIMAPSampler" ),
                      paste( sampler_options, collapse = ' ' ), sep = ' ' )
    message( 'Running BI-MAP MPI sampler: ', cmdline, '...' )
    system( cmdline )
    message( 'BI-MAP MPI sampler finished' )

    unlink( tempOPAfilename ) # delete temporarily data file

    # load result
    res <- NULL
    if ( walk.create.RObject ) {
        res <- BIMAP.walk.load( walk.file )
    }
    if ( !is.null( tempWalk.file ) ) unlink( tempWalk.file ) # delete temporary result file (saved to dataframe)
    return ( res )
}

# Loads BIMAPWalk from file,
# optionally identifying independent components with given threshold
BIMAP.mcmcwalk.load <- function ( filename, 
    objects.components.threshold = NA,
    probes.components.threshold = NA,
    objects.clot.threshold = NA,
    stable.objects.clusters.threshold = NA,
    stable.probes.clusters.threshold = NA
){
    return ( .Call( "BIMAPWalkLoad", filename, 
                    objects.components.threshold,
                    probes.components.threshold,
                    objects.clot.threshold,
                    stable.objects.clusters.threshold,
                    stable.probes.clusters.threshold
 ) )
}

setClass( "BIMAPCluster",
    representation = representation( 
        serial = "integer",
        proteins = "list",
        samples = "list"
    )
)

BIMAPCluster <- function( proteins, samples, serial = 1 )
{
    return ( new( "BIMAPCluster", serial = as.integer( serial ), proteins = proteins, samples = samples ) )
}

setGeneric( "proteins.data.frame", function( .Object, ... ) standardGeneric( "proteins.data.frame" ) )

setMethod( "proteins.data.frame",
    signature( .Object = "BIMAPCluster" ), 
    function( .Object )
    {
        nprots <- length( .Object@proteins )
        return( data.frame( cluster = rep( .Object@serial, nprots ), 
                            protein = names( .Object@proteins ), 
                            count = as.numeric( .Object@proteins ) ) )
    }
)

setClass( "BIMAPClustering",
    representation = representation( 
        serial = "integer",
        objects.partition.serial = "integer",
        probes.partition.serial = "integer",
        objects.clusters = "data.frame",
        probes.clusters = "data.frame",
        blocks = "data.frame",
        probes.signals = "data.frame",
        baseline.peak = "numeric", baseline.shape = "numeric",
        noise.peak = "numeric", noise.shape = "numeric"
    )
)

setClass( "BIMAPWalk",
    representation = representation( 
        prob.weights = "data.frame",
        clusterings.walk = "data.frame",
        priors.walk = "data.frame",
        clusterings = "data.frame",
        objects.clusters = "data.frame", # elements
        objects.clusters.info = "data.frame",
        objects.partitions = "data.frame",
        objects.subpartitions = "data.frame",
        objects.components = "data.frame",
        objects.data = "data.frame",
        probes.clusters = "data.frame", # elements
        probes.clusters.info = "data.frame",
        probes.partitions = "data.frame",
        probes.subpartitions = "data.frame",
        probes.components = "data.frame",
        blocks = "data.frame",
        blocks.freq = "data.frame",
        signals = "data.frame",
        stable.blocks.scores = "data.frame"
    )
)

MSMatrixToFrame <- function( msmatrix )
{
    res <- NULL
    for ( msrun in colnames( msmatrix ) ) {
        for ( protein in rownames( msmatrix) ) {
            df <- data.frame( msrun = c( msrun ), prey_ac = c( protein ), sc = c( msmatrix[[protein,msrun]] ) )
            if ( is.null( res ) ) {
                res <- df
            } else {
                res <- rbind( res, df )
            }
        }
    }
    return ( res )
}


BIMAP.mcmcwalk.biclusterings_dataframe <- function( walk )
{
    data.frame(
        iteration = as.numeric( names( walk@clusterings ) ),
        serial = factor( as.character( lapply( walk@clusterings, function( clus ) clus@serial ) ) ),
        objects.partition.serial = factor( as.character( lapply( walk@clusterings, function( clus ) clus@objects.partition.serial ) ) ),
        probes.partition.serial = factor( as.character( lapply( walk@clusterings, function( clus ) clus@probes.partition.serial ) ) ),
        baseline.peak = as.numeric( lapply( walk@clusterings, function( clus ) clus@baseline.peak ) ),
        noise.peak = as.numeric( lapply( walk@clusterings, function( clus ) clus@noise.peak ) )
    )
}

BIMAP.mcmcwalk.priors_dataframe <- function( walk )
{
    data.frame(
        iteration = as.numeric( names( walk@priors ) ),
        #noise.sigma = as.numeric( lapply( walk@priors, function( prior ) prior@noise.sigma ) ),
        baseline.signal = as.numeric( lapply( walk@priors, function( prior ) prior@baseline.signal ) ),
        clusters.signal.sigma = as.numeric( lapply( walk@priors, function( prior ) prior@clusters.signal.sigma ) )
    )
}

BIMAP.mcmcwalk.blocks_dataframe <- function( walk )
{
    res <- do.call( 'rbind', lapply( names( walk@clusterings ), function ( it ) {
        do.call( 'rbind', lapply( walk@clusterings[[ it ]]@blocks, function ( cluster ) {
            cluframe <- data.frame( clustering@blocks )
            cluframe <- cbind( cluframe, iteration = as.integer(rep(as.numeric(it), nrow(cluframe))) )
            return ( cluframe )
        } ) )
    } ) )
    res$iteration <- as.factor( res$iteration )
    res$objects.cluster.serial <- as.factor( res$objects.cluster.serial )
    res$probes.cluster.serial <- as.factor( res$probes.cluster.serial )
    return ( res )
}

BIMAP.mcmcwalk.signals_dataframe <- function( walk, per.object = TRUE )
{
    res <- do.call( 'rbind', lapply( names( walk@clusterings ), function( it ) {
        clustering <- walk@clusterings[[ it ]]
        cluframe <- data.frame( biclustering.serial = rep( clustering@serial, 
                                length( clustering@probes.signals$objects.cluster.serial ) ),
                                objects.cluster.serial = clustering@probes.signals$objects.cluster.serial,
                                probe = clustering@probes.signals$probe,
                                probe.signal = clustering@probes.signals$probe.signal
                              )
        if ( per.object ) {
            objptnframe <- data.frame(
                objects.cluster.serial = clustering@objects.clusters$objects.cluster.serial,
                object = clustering@objects.clusters$object,
                object.multiple = clustering@objects.clusters$object.multiple
            )
            cluframe <- merge( cluframe, objptnframe, by = c('objects.cluster.serial' ) )
        }
        cluframe <- cbind( cluframe, iteration = rep(as.integer(it), nrow(cluframe)) )
        return ( cluframe )
    } ) )
    res$iteration <- as.factor( res$iteration )
    res$probe <- as.factor( res$probe )
    res$objects.cluster.serial <- as.factor( res$objects.cluster.serial )
    if ( per.object ) {
        res$object <- as.factor( res$object )
    }
    return ( res )
}

