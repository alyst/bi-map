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

BIMAP.distances <- function(
    protein_info, sample_info,
    msrun_info, measurements = NULL,
    signal.sequence.length.factor = 0.5,
    signal.shape = 0.0,
    precomp.object.freq.threshold = 0.6,
    precomp.probe.freq.threshold = 0.6
){
    params <- list(
        signal.sequence.length.factor = signal.sequence.length.factor,
        signal.shape = signal.shape,
        precomp.object.freq.threshold = precomp.object.freq.threshold,
        precomp.probe.freq.threshold = precomp.probe.freq.threshold
    )
    return ( .Call( "CalcDistances",
            protein_info, sample_info, msrun_info, measurements,
            params ) )
}

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
    protein_ac_column = 'primaryac', protein_seqlength_column = 'seqlength'
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
    proteins.df <- protein_info[ protein_info[ , protein_ac_column ] %in% 
                union( as.character( ms_data_noglob[ , prey_column ] ), 
                    as.character( ms_data_noglob[ , bait_column ] ) ), 
                                         c( protein_ac_column, protein_seqlength_column ) ]
    proteins.df <- ddply( proteins.df, c( protein_ac_column ), function( prot_data ) {
                data.frame( seqlength = max( prot_data[,protein_seqlength_column] ),
                            stringsAsFactors = FALSE )
    } )
    #print( as.character( ms_data_noglob$prey_ac ) )
    colnames( proteins.df ) <- c( 'protein_ac', 'seqlength' )
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
    format = c( 'OPA', 'CSV' )
){
    if ( format == 'OPA' ){
        return ( .Call( "OPADataSave", filename,
                bimap.data$proteins, bimap.data$samples, bimap.data$msruns,
                bimap.data$measurements ) )
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
    log.particles.file = NULL,
    log.eeJumps.file = NULL
){
    params <- list(
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
	objects.clot.threshold = NA )
{
    return ( .Call( "BIMAPWalkLoad", filename, 
                    objects.components.threshold,
                    probes.components.threshold,
                    objects.clot.threshold ) )
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
        signals = "data.frame"
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

BIMAP.mcmcwalk.lastClustering <- function( walk, serial )
{
    for ( it in rev(names(walk@clusterings)) ) {
        clu <- walk@clusterings[[it]] 
        if ( clu@serial == serial ) {
            return ( clu )
        }
    }
    warning(paste('Clustering with serial #',serial,' not found',sep=''))
}

#' Order clusterings according to their posterior probability.
#' @param walk 
#' @param n number of top probable clusterings to return
#' @returnType 
#' @return data.frame with clusters information
#' @author astukalov
#' @export
BIMAP.mcmcwalk.biclusterings_order <- function( walk, n = 30, by = c( 'parts.pdf', 'components.pdf' ) ) {
    clus.df <- walk@clusterings
    rownames( clus.df ) <- clus.df$clustering_id
    if ( by == 'parts.pdf' ) {
       clu_pdf <- clus.df$objects.parts.lnpdf + clus.df$probes.parts.lnpdf
    } else if ( by == 'components.pdf' ){
       clu_pdf <- clus.df$objects.lnpdf + clus.df$probes.lnpdf
    }
    clus.df <- clus.df[ order( clu_pdf,
            #clus.df$probes.lnpdf,
            clus.df$total.lnpdf,
            decreasing = TRUE ), ]
    if ( is.numeric( n ) ) {
        clus.df <- clus.df[ 1:min(n, nrow(clus.df)), ]
    }
    return ( clus.df )
}

BIMAP.mcmcwalk.clusters_counts <- function( walk, entities = c( 'objects', 'probes' ) ) {
    partition.serial.colname <- paste( entities, 'partition.serial', sep = '.' )
    cluster.serial.colname <- paste( entities, 'cluster.serial', sep = '.' )
    partitions.walk <- subset( merge( walk@clusterings.walk, walk@clusterings, by = 'clustering.serial' ), select = partition.serial.colname )
    clusters.walk <- subset( merge( partitions.walk, slot( walk, paste( entities, 'partitions', sep='.' ) ), by = partition.serial.colname ), select = c( cluster.serial.colname ) )
    return ( table( clusters.walk[, cluster.serial.colname ] ) )
}

BIMAP.mcmcwalk.clusters_size <- function( walk, entities = c( 'objects','probes' ) ) {
    cluster.serial.colname <- paste( entities, 'cluster.serial', sep = '.' )
    clusters.contents <- slot( walk, paste( entities, 'clusters', sep='.' ) )
    cluster.sizes.df <- ddply( clusters.contents, c( cluster.serial.colname), nrow )
    clusters.sizes <- cluster.sizes.df[,2]
    names(clusters.sizes) <- cluster.sizes.df[,1]
    return ( clusters.sizes )
}

# prepare pivot table of signal means and SDs
BIMAP.signal_stats <- function( signals.frame )
{
    fix.array <- function( arr ) {
        res <- as.matrix( arr )
        dimnames( res ) <- dimnames( arr )[1:2]
        return ( res )
    }
    # prepare pivot table of signal means and SDs
    signals.pivot.mean <- fix.array( daply( signals.frame,# .drop_o = FALSE,
        .(proteins.cluster, samples.cluster ),
        function( rows ) mean( rows$signal ) ) )
    
    signals.pivot.sd <- fix.array( daply( signals.frame,# .drop_o = FALSE,
        .(proteins.cluster, samples.cluster ),
        function( rows ) sd( rows$signal ) ) )
    return ( c( list( signals.mean = signals.pivot.mean,
                      signals.sd = signals.pivot.sd ) ) )
}

#' Extracts chessboard biclustering with specified ID from random walk.
#' @param bimap.walk sampling walk
#' @param bimapId ID of chessboard biclustering to extract
#' @param onblock.threshold "on" blocks, are those whose "on" probe frequency is greater than this threshold
#' @param extract.signals extract signals of individual blocks
#' @returnType 
#' @return 
#' @author astukalov
BIMAP.mcmcwalk.extract_biclustering <- function( bimap.walk, bimapId,
    onblock.threshold = 0.6,
    extract.signals = TRUE
){
    steps <- subset( bimap.walk@clusterings.walk, clustering.serial == bimapId )$step
    if ( length(steps) == 0 ) {
        stop( paste("No chessboard biclustering with ID=", bimapId, "found") )
    }

    # get object and probe clusters of chessboard biclustering
    opId <- subset( bimap.walk@clusterings, clustering.serial == bimapId )$objects.partition.serial
    ppId <- subset( bimap.walk@clusterings, clustering.serial == bimapId )$probes.partition.serial
    ocIds <- as.character( subset( bimap.walk@objects.partitions, objects.partition.serial == opId )$objects.cluster.serial )
    pcIds <- as.character( subset( bimap.walk@probes.partitions, probes.partition.serial == ppId )$probes.cluster.serial )

    if ( is.numeric(onblock.threshold) ) {
        # build consensus blocks
        blockStats <- subset( bimap.walk@blocks.freq, objects.cluster.serial %in% ocIds & probes.cluster.serial %in% pcIds )
        blocks <- subset( blockStats, enabled >= onblock.threshold * total,
                             select=c( 'objects.cluster.serial', 'probes.cluster.serial' ) )
    } else {
        # get non-empty blocks of the clustering
        blocks <- subset( bimap.walk@blocks, clustering.serial == bimapId,
                          select = c("objects.cluster.serial", "probes.cluster.serial") )
    }

    nBlocks <- nrow( blocks )
    if ( nBlocks == 0 ) {
        stop( paste("No on-blocks found in clustering ID=", bimapId ) )
    } else {
        message( nBlocks, ' on-blocks found of ', length( ocIds ), 'x', length( pcIds ), ' possible' )
    }
    colnames( blocks ) <- c( 'proteins.cluster', 'samples.cluster' )
    blocks$proteins.cluster <- as.character( blocks$proteins.cluster )
    blocks$samples.cluster <- as.character( blocks$samples.cluster )

    proteins.clusters <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% ocIds,
        select = c( 'objects.cluster.serial', 'object' ) )
    colnames( proteins.clusters ) <- c( 'proteins.cluster', 'protein_ac' )
    proteins.clusters$proteins.cluster <- as.character( proteins.clusters$proteins.cluster )
    proteins.clusters$protein_ac <- as.character( proteins.clusters$protein_ac )
    rownames( proteins.clusters ) <- proteins.clusters$protein_ac
    samples.clusters <- subset( bimap.walk@probes.clusters, probes.cluster.serial %in% pcIds,
        select = c( 'probes.cluster.serial', 'probe' ) )
    colnames( samples.clusters ) <- c( 'samples.cluster', 'sample' )
    samples.clusters$sample <- as.character( samples.clusters$sample )
    samples.clusters$samples.cluster <- as.character( samples.clusters$samples.cluster )
    rownames( samples.clusters ) <- samples.clusters$sample

    res <- list(
            steps = steps,
            blocks = blocks,
            proteins.clusters = proteins.clusters,
            proteins.clusters.info = bimap.walk@objects.clusters.info[ unique( proteins.clusters$proteins.cluster ),
                                                                    c( 'objects.cluster.serial', 'size',
                                                                       'nsteps', 'nsteps.included', 'avg.pairs.cooccur') ],
            samples.clusters = samples.clusters,
            samples.clusters.info = bimap.walk@probes.clusters.info[ unique( samples.clusters$samples.cluster ),
                                                                 c( 'probes.cluster.serial', 'size',
                                                                    'nsteps', 'nsteps.included', 'avg.pairs.cooccur') ]
    )
    colnames( res$proteins.clusters.info ) <- c( 'proteins.cluster', 'size', 'nsteps', 'nsteps_included', 'avg_pairs_cooccur' )
    res$proteins.clusters.info$avg_pairs_freq <- res$proteins.clusters.info$avg_pairs_cooccur / nrow( bimap.walk@clusterings.walk )
    colnames( res$samples.clusters.info ) <- c( 'samples.cluster', 'size', 'nsteps', 'nsteps_included', 'avg_pairs_cooccur' )
    res$samples.clusters.info$avg_pairs_freq <- res$samples.clusters.info$avg_pairs_cooccur / nrow( bimap.walk@clusterings.walk )
    if ( extract.signals ) {
        message( 'Extracting signals for on-blocks of bi-clustering #', bimapId, '...' )
        # get subframe with all samples for signals of given bimap
        # (but samples might be from other bimaps, which contain the same clusters)
        block_ids <- paste( blocks$proteins.cluster, blocks$samples.cluster )
        signals.subframe <- subset( bimap.walk@signals,
                paste( objects.cluster.serial, probes.cluster.serial ) %in% block_ids )[c('step', 'objects.cluster.serial', 'probes.cluster.serial', 'signal')]
        message( nrow(signals.subframe), " signals extracted" )
        colnames( signals.subframe ) <- c( 'step', 'proteins.cluster', 'samples.cluster', 'signal' )
        signals.subframe$proteins.cluster <- as.character( signals.subframe$proteins.cluster )
        signals.subframe$samples.cluster <- as.character( signals.subframe$samples.cluster )
        message( "Calculating on-blocks signals statistics..." )
        signal_stats <- BIMAP.signal_stats( signals.subframe )
        res$blocks$signal.mean <- apply( res$blocks, 1,
                function( block ) signal_stats$signals.mean[ block[['proteins.cluster']], block[['samples.cluster']] ] )
        res$blocks$signal.sd = apply( res$blocks, 1,
                function( block ) signal_stats$signals.sd[ block[['proteins.cluster']], block[['samples.cluster']] ] )
        res <- c( res, signal_stats )
        res$signals.subframe = signals.subframe
        message( 'On-blocks signals statistics calculation done' )
    }
    return ( res )
}

BIMAP.biclustering.filter <- function( bimap.props, protein_acs, sample_acs )
{
    # filter elements and clusters
    res <- list( #proteins = subset( bimap.props$proteins, protein_ac %in% protein_acs ),
                 #samples = subset( bimap.props$samples, sample %in% sample_acs ),
                 proteins.clusters = subset( bimap.props$proteins.clusters, protein_ac %in% protein_acs ),
                 samples.clusters = subset( bimap.props$samples.clusters, sample %in% sample_acs ) )
    samples.cluster_acs <- unique( res$samples.clusters$samples.cluster )
    proteins.cluster_acs <- unique( res$proteins.clusters$proteins.cluster )
    # remove blocks corresponding to removed clusters
    res$blocks <- subset( bimap.props$blocks,
                           proteins.cluster %in% proteins.cluster_acs &
                           samples.cluster %in% samples.cluster_acs )
    res$signals.subframe <- subset( bimap.props$signals.subframe,
                                  proteins.cluster %in% proteins.cluster_acs &
                                  samples.cluster %in% samples.cluster_acs )
    signal_stats <- BIMAP.signal_stats( res$signals.subframe )
    for ( i in 1:nrow(res$blocks) ) {
        blockRow <- res$blocks[ i, ] 
        res$blocks$signal.mean = signal_stats$signals.mean[ blockRow$proteins.cluster,
                                                             blockRow$samples.cluster ] 
        res$blocks$signal.sd = signal_stats$signals.sd[ blockRow$proteins.cluster,
                                                        blockRow$samples.cluster ]
    }
    res$signals.mean <- signal_stats$signals.mean
    res$signals.sd <- signal_stats$signals.sd
    return ( res )
}

BIMAP.mcmcwalk.objects_cluster_labels <- function( walk )
{
    df <- data.frame( serial = c(), label = c() )
    for ( clustering in walk@clusterings ) {
        if ( length( clustering@proteins.clusters$objects.cluster.serial ) == 0 ) next;
        ocs <- as.data.frame( clustering@proteins.clusters )
        df <- unique( rbind( df, ddply( ocs, .(objects.cluster.serial), function( objclu ) {
            objs <- objclu[,'object']
            objs <- objs[ order(objs) ]
            data.frame( serial = objclu[[1,'objects.cluster.serial']],
                        label = paste( objs, collapse = ',' ) )
        } ) ) )
    }
    res <- as.character( df$label )
    names( res ) <- as.integer( df$serial )
    return ( res )
}

BIMAP.mcmcwalk.objects_partition_labels <- function( walk )
{
    res <- list()
    for ( clustering in walk@clusterings ) {
        if ( length( clustering@proteins.clusters$objects.cluster.serial ) == 0 ) next;
        ocs <- as.data.frame( clustering@proteins.clusters )
        clunames <- as.character( ddply( ocs, .(objects.cluster.serial), function( objclu ) {
                        objs <- objclu[,'object']
                        objs <- objs[ order(objs) ]
                        data.frame( serial = objclu[[1,'objects.cluster.serial']],
                                    label = paste( objs, collapse = '' ) )
                    } )$label )
        clunames <- clunames[ order( clunames ) ]
        res[[ as.character(clustering@objects.partition.serial) ]] = paste( clunames, collapse = ' ' )
    }
    return ( res )
}

BIMAP.mcmcwalk.probes_cluster_labels <- function( walk )
{
    df <- data.frame( serial = c(), label = c() )
    for ( clustering in walk@clusterings ) {
        if ( length( clustering@samples.clusters$probes.cluster.serial ) == 0 ) next;
        ocs <- as.data.frame( clustering@samples.clusters )
        df <- unique( rbind( df, ddply( ocs, .(probes.cluster.serial), function( probesclu ) {
                probes <- probesclu[,'probe']
                probes <- probes[ order(probes ) ]
                data.frame( serial = probesclu[[1,'probes.cluster.serial']], label = paste( probes, collapse = ',' ) )
            } ) ) )
    }
    res <- as.character( df$label )
    names( res ) <- as.integer( df$serial ) 
    return ( res )
}

BIMAP.mcmcwalk.probes_partition_labels <- function( walk )
{
    res <- list()
    for ( clustering in walk@clusterings ) {
        if ( length( clustering@samples.clusters$probes.cluster.serial ) == 0 ) next;
        ocs <- as.data.frame( clustering@samples.clusters )
        clunames <- as.character( ddply( ocs, .(probes.cluster.serial), function( probeclu ) {
                    probes <- probeclu[,'probe']
                    probes <- probes[ order(probes) ]
                    data.frame( serial = probeclu[[1,'probes.cluster.serial']], label = paste( probes, collapse = '' ) )
                } )$label )
        clunames <- clunames[ order( clunames ) ]
        res[[ as.character(clustering@probes.partition.serial) ]] = paste( clunames, collapse = ' ' )
    }
    return ( res )
}


#' Extracts table of frequencies of given protein interactions.
#' @param bimap.walk 
#' @param protein.ac 
#' @param protein.info 
#' @param min.freq 
#' @returnType 
#' @return 
#' @author astukalov
#' @export
BIMAP.mcmcwalk.interactors_frame <- function( bimap.walk, protein.ac, protein.info =-NULL, min.freq = 0.01 )
{
    protein.clusters <- subset( bimap.walk@proteins.clusters, object == protein.ac )$objects.cluster.serial
    proteins.clusterings.walk <- subset( merge( bimap.walk@clusterings.walk, bimap.walk@clusterings ), select = c('objects.partition.serial') ) 
    protein.clusters.walk <- subset( merge( proteins.clusterings.walk, bimap.walk@objects.partitions, by = 'objects.partition.serial' ), 
        objects.cluster.serial %in% protein.clusters, select = c( 'objects.cluster.serial' ) )
    interactors.walk <- subset( merge( bimap.walk@proteins.clusters, protein.clusters.walk, by = c('objects.cluster.serial') ), select = c( 'object' ) )
    res <- table( interactors.walk$object ) / nrow( bimap.walk@clusterings.walk )
    res <- as.table( sort( res[ res >= min.freq ], decreasing = TRUE ) )
    if ( !is.null( protein.info ) ) {
        names( res ) <- protein.info[ names( res ), 'short_id' ]
    }
    return ( res )
}

#' Extract most stable clusters using greedy approach.
#' @param bimap.walk 
#' @param ms_data original MS dataset
#' @param sample_col sample column in MS dataset
#' @param prey_col prey AC column in MS dataset
#' @param min.avg.nsample average number of samples the cluster is seen
#' @param min.size minimal size of cluster
#' @param min.included.freq minimal frequency the cluster is seen (itself or as subcomponent) in partitions
#' @param size.weight weight of cluster in scoring
#' @returnType dataframe of clusters
#' @return 
#' @author astukalov
#' @export
BIMAP.mcmcwalk.extract_stable_clusters <- function( bimap.walk, bimap.data,
        min.avg.nsample = 3, min.size = 2,
        min.included.freq = 0.95, size.weight = 0.2,
        allow.intersections = FALSE )
{
    min.nsteps.included = nrow( bimap.walk@clusterings.walk ) * min.included.freq
    good.clusters.info <- subset( bimap.walk@objects.clusters.info, size >= min.size & nsteps.included >= min.nsteps.included )
    good.clusters.info$score <- ( good.clusters.info$nsteps.included / nrow( bimap.walk@clusterings.walk ) - min.included.freq ) / ( 1.0 - min.included.freq ) +
                  size.weight * ( good.clusters.info$size - min.size )
    #print( good.clusters.info )
    good.clusters <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% good.clusters.info$objects.cluster.serial )
    avg.nsamples <- ddply( unique( merge( good.clusters, merge( bimap.data$measurements, bimap.data$exp_design ),
                                          by.x = 'object', by.y = 'prey_ac' )
                            [ c( 'objects.cluster.serial', 'object', 'sample' ) ] ), .( objects.cluster.serial ), function( clu_ms_data ) {
                data.frame( avg.nsample = min( table( clu_ms_data$object ) ),
                        stringsAsFactors = TRUE )
            } )
    good.clusters.info <- merge( good.clusters.info, avg.nsamples, by = c( 'objects.cluster.serial' ) )
    good.clusters.info <- subset( good.clusters.info, avg.nsample >= min.avg.nsample )
    good.clusters <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% good.clusters.info$objects.cluster.serial )
    intersecting.pairs <- unique( merge( good.clusters, good.clusters, by = c( 'object' ),
                    all = FALSE, sort = FALSE )
                    [ c( 'objects.cluster.serial.x', 'objects.cluster.serial.y' ) ] )
    if ( !allow.intersections ) {
        best.clusters.serials <- c()
        while ( nrow( good.clusters.info ) > 0 ) {
            best.cluster <- good.clusters.info[ order( good.clusters.info$score, decreasing = TRUE )[1], 'objects.cluster.serial' ]
            best.clusters.serials <- c( best.clusters.serials, best.cluster )
            intersects.with <- subset( intersecting.pairs, objects.cluster.serial.x == best.cluster )$objects.cluster.serial.y
            good.clusters.info <- subset( good.clusters.info, !( objects.cluster.serial %in% intersects.with ) )
        }
        res <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% best.clusters.serials )
    } else {
        res <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% good.clusters.info$objects.cluster.serial )
    }
    res <- res[ order( res$objects.cluster.serial, res$object ), ]
    return ( res )
}
