# BIMAP R interface basic definitions.
# 
# Author: Alexey Stukalov (astukalov@cemm.oeaw.ac.at)
###############################################################################

require( 'plyr' )
require( 'Rcpp' )

# load RBIMAP library
#dyn.unload( file.path( RBIMAP.libpath, paste("libRBIMAP", .Platform$dynlib.ext, sep="")) ) 
dyn.load( file.path( RBIMAP.libpath, paste("libRBIMAP", .Platform$dynlib.ext, sep="")), type = "Call" ) 

#' Saves AP-MS data into single .xml file that could be read by BIMAP-sampler --input_file
#' @param filename name of output file
#' @param protein_info protein information dataframe (AC, Sequence Length)
#' @param sample_info  sample information dataframe
#' @param msrun_info MS runs information dataframe
#' @param measurements MS measurements dataframe
#' @return result of call to external OPADataSave()
#' @author Alexey Stukalov
BIMAP.save_apms_data <- function(
    filename,
    protein_info, sample_info,
    msrun_info, measurements = NULL )
{
    return ( .Call( "OPADataSave", filename,
                    protein_info, sample_info, msrun_info,
                    measurements ) )
}

BIMAP.distances <- function(
    protein_info, sample_info,
    msrun_info, measurements = NULL,
    signal.sequence.length.factor = 0.5,
    signal.shape = 0.0,
    precomp.object.freq.threshold = 0.6,
    precomp.state.freq.threshold = 0.6
){
    params <- list(
        signal.sequence.length.factor = signal.sequence.length.factor,
        signal.shape = signal.shape,
        precomp.object.freq.threshold = precomp.object.freq.threshold,
        precomp.state.freq.threshold = precomp.state.freq.threshold
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
BIMAP.import_msdata <- function( ms_data, protein_info, msrun.multipliers = NULL,
    sample_column = 'sample', msrun_column = 'msrun',
    bait_column = 'bait_ac', prey_column = 'prey_ac',
    sc_column = 'sc', pc_column = 'pc',
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
    samples.df <- unique( subset( ms_data_noglob, select = unique( c( sample_column, bait_column ) ) ) )
    samples.df <- data.frame( 
            sample = as.character( samples.df[, sample_column ] ), 
            bait_ac = as.character( samples.df[, bait_column ] ),
            stringsAsFactors = FALSE )
    rownames( samples.df ) <- samples.df$sample 
    msruns.df <- as.data.frame( unique( subset( ms_data_noglob, select = unique( c( msrun_column, sample_column ) ) ) ),
                                stringsAsFactors = FALSE )
    msruns.df <- data.frame( 
            msrun = as.character( msruns.df[, msrun_column ] ), 
            sample = as.character( msruns.df[, sample_column ] ),
            stringsAsFactors = FALSE )
    rownames( msruns.df ) <- msruns.df$msrun 
    if ( !is.null(msrun.multipliers) ) {
        mean.mult <- mean( msrun.multipliers[ msruns.df$msrun ], na.rm = TRUE )
        msruns.df$multiplier <- msrun.multipliers[ msruns.df$msrun ] / mean.mult
        msruns.df[ is.na(msruns.df$multiplier), 'multiplier' ] <- 1
    }
    proteins.df <- unique( protein_info[ protein_info[ , protein_ac_column ] %in% 
                union( as.character( ms_data_noglob[ , prey_column ] ), 
                    as.character( ms_data_noglob[ , bait_column ] ) ), 
                                         c( protein_ac_column, protein_seqlength_column ) ] ) 
         
    #print( as.character( ms_data_noglob$prey_ac ) )
    colnames( proteins.df ) <- c( 'protein_ac', 'seqlength' )
    rownames( proteins.df ) <- proteins.df$protein_ac 

    return ( list(  measurements = measurements.df,
                    samples = samples.df, 
                    msruns = msruns.df, 
                    proteins = proteins.df ) )
}

BIMAP.walk.eval <- function(
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
    sampleRate.state.membership = 0.3,
    sampleRate.states.splitMerge = 0.3,
    sampleRate.block.flip = 0.05,
    sampleRate.object.multiple = 0.125,
    sampleRate.signal = 0.1,
    samplePeriod.priors = 197,
    samplePeriod.chessboardBiclustering = 23,
    block.resamples = 12,
    ini.objects.partition = NULL,
    ini.states.partition = NULL,
    ini.blocks = NULL,
    signal.sequence.length.factor = 0.5,
    signal.shape = 0.0,
    prior.objects.clustering.concentration = 0.1,
    prior.objects.clustering.discount = 0.0,
    prior.states.clustering.concentration = 0.1,
    prior.states.clustering.discount = 0.0,
    hyperprior.baseline = 0.0,
    hyperprior.baseline.scale = 1.0,
    hyperprior.signal.variance.shape = 1.0,
    hyperprior.signal.variance.scale = 0.1,
    prior.true_misses = 1000,
    prior.false_hits = 1,
    prior.block.enabled = 0.1,
    precomp.object.freq.threshold = 0.6,
    precomp.state.freq.threshold = 0.6,
    objects.components.threshold = 0.1,
    states.components.threshold = 0.1,
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
        sampleRate.state.membership = sampleRate.state.membership,
        sampleRate.states.splitMerge = sampleRate.states.splitMerge,
        sampleRate.block.flip = sampleRate.block.flip,
        sampleRate.object.multiple = sampleRate.object.multiple,
        sampleRate.signal = sampleRate.signal,
        samplePeriod.priors = samplePeriod.priors,
        samplePeriod.chessboardBiclustering = samplePeriod.chessboardBiclustering,
        block.resamples = block.resamples,
        ini.objects.partition = ini.objects.partition,
        ini.states.partition = ini.states.partition,
        ini.blocks = ini.blocks,
        signal.sequence.length.factor = signal.sequence.length.factor,
        signal.shape = signal.shape,
        prior.objects.clustering.concentration = prior.objects.clustering.concentration,
        prior.objects.clustering.discount = prior.objects.clustering.discount,
        prior.states.clustering.concentration = prior.states.clustering.concentration,
        prior.states.clustering.discount = prior.states.clustering.discount,
        hyperprior.baseline = hyperprior.baseline,
        hyperprior.baseline.scale = hyperprior.baseline.scale,
        hyperprior.signal.variance.shape = hyperprior.signal.variance.shape,
        hyperprior.signal.variance.scale = hyperprior.signal.variance.scale,
        prior.true_misses = prior.true_misses,
        prior.false_hits = prior.false_hits,
        prior.block.enabled = prior.block.enabled,
        precomp.object.freq.threshold = precomp.object.freq.threshold,
        precomp.state.freq.threshold = precomp.state.freq.threshold,
        objects.components.threshold = objects.components.threshold,
        states.components.threshold = states.components.threshold,
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
BIMAP.walk.MPI_eval <- function(
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
    tempOSAfilename <- paste( tempfile( "osadata_" ), ".xml.gz", sep='' )
    OSAData.save( tempOSAfilename, protein_info, sample_info, msrun_info, measurements )
    sampler_options <- c( '--input_file', tempOSAfilename )
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
    print( paste('Running MPI sampler: ', cmdline ) )
    system( cmdline )
    print( paste( 'Finished MPI sampling' ) )

    unlink( tempOSAfilename ) # delete temporarily data file

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
BIMAP.walk.load <- function ( filename, 
    objects.components.threshold = NA,
    states.components.threshold = NA )
{
    return ( .Call( "BIMAPWalkLoad", filename, 
                    objects.components.threshold, states.components.threshold ) )
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
        states.partition.serial = "integer",
        objects.clusters = "data.frame",
        states.clusters = "data.frame",
        blocks = "data.frame",
        states.signals = "data.frame",
        baseline.peak = "numeric", baseline.shape = "numeric",
        noise.peak = "numeric", noise.shape = "numeric"
    )
)

setClass( "BIMAPWalk",
    representation = representation( 
        clusterings.walk = "data.frame",
        priors.walk = "data.frame",
        clusterings = "data.frame",
        objects.clusters = "data.frame", # elements
        objects.clusters.info = "data.frame",
        objects.partitions = "data.frame",
        objects.subpartitions = "data.frame",
        objects.components = "data.frame",
        objects.data = "data.frame",
        states.clusters = "data.frame", # elements
        states.clusters.info = "data.frame",
        states.partitions = "data.frame",
        states.subpartitions = "data.frame",
        states.components = "data.frame",
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


BIMAP.walk.data.frame <- function( walk )
{
    data.frame(
        iteration = as.numeric( names( walk@clusterings ) ),
        serial = factor( as.character( lapply( walk@clusterings, function( clus ) clus@serial ) ) ),
        objects.partition.serial = factor( as.character( lapply( walk@clusterings, function( clus ) clus@objects.partition.serial ) ) ),
        states.partition.serial = factor( as.character( lapply( walk@clusterings, function( clus ) clus@states.partition.serial ) ) ),
        baseline.peak = as.numeric( lapply( walk@clusterings, function( clus ) clus@baseline.peak ) ),
        noise.peak = as.numeric( lapply( walk@clusterings, function( clus ) clus@noise.peak ) )
    )
}

BIMAP.walk.priors.data.frame <- function( walk )
{
    data.frame(
        iteration = as.numeric( names( walk@priors ) ),
        #noise.sigma = as.numeric( lapply( walk@priors, function( prior ) prior@noise.sigma ) ),
        baseline.signal = as.numeric( lapply( walk@priors, function( prior ) prior@baseline.signal ) ),
        clusters.signal.sigma = as.numeric( lapply( walk@priors, function( prior ) prior@clusters.signal.sigma ) )
    )
}

BIMAP.walk.crosscluster.data.frame <- function( walk )
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
    res$states.cluster.serial <- as.factor( res$states.cluster.serial )
    return ( res )
}

BIMAP.walk.signals.data.frame <- function( walk, per.object = TRUE )
{
    res <- do.call( 'rbind', lapply( names( walk@clusterings ), function( it ) {
        clustering <- walk@clusterings[[ it ]]
        cluframe <- data.frame( cross.clustering.serial = rep( clustering@serial, 
                                length( clustering@states.signals$objects.cluster.serial ) ),
                                objects.cluster.serial = clustering@states.signals$objects.cluster.serial,
                                state = clustering@states.signals$state,
                                state.signal = clustering@states.signals$state.signal
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
    res$state <- as.factor( res$state )
    res$objects.cluster.serial <- as.factor( res$objects.cluster.serial )
    if ( per.object ) {
        res$object <- as.factor( res$object )
    }
    return ( res )
}

BIMAP.walk.lastClustering <- function( walk, serial )
{
    for ( it in rev(names(walk@clusterings)) ) {
        clu <- walk@clusterings[[it]] 
        if ( clu@serial == serial ) {
            return ( clu )
        }
    }
    warning(paste('Clustering with serial #',serial,' not found',sep=''))
}

BIMAP.walk.clusters.counts <- function( walk, entities = c( 'objects', 'states' ) ) {
    partition.serial.colname <- paste( entities, 'partition.serial', sep = '.' )
    cluster.serial.colname <- paste( entities, 'cluster.serial', sep = '.' )
    partitions.walk <- subset( merge( walk@clusterings.walk, walk@clusterings, by = 'clustering.serial' ), select = partition.serial.colname )
    clusters.walk <- subset( merge( partitions.walk, slot( walk, paste( entities, 'partitions', sep='.' ) ), by = partition.serial.colname ), select = c( cluster.serial.colname ) )
    return ( table( clusters.walk[, cluster.serial.colname ] ) )
}

BIMAP.walk.clusters.size <- function( walk, entities = c( 'objects','states' ) ) {
    cluster.serial.colname <- paste( entities, 'cluster.serial', sep = '.' )
    clusters.contents <- slot( walk, paste( entities, 'clusters', sep='.' ) )
    cluster.sizes.df <- ddply( clusters.contents, c( cluster.serial.colname), nrow )
    clusters.sizes <- cluster.sizes.df[,2]
    names(clusters.sizes) <- cluster.sizes.df[,1]
    return ( clusters.sizes )
}

BIMAP.walk.greedy.partition <- function( walk, entities = c( 'objects','states' ), score = function( size, counts ) return ( counts ) ) {
    clusters.sizes <- BIMAP.walk.clusters.size( walk, entities )
    clusters.counts <- BIMAP.walk.clusters.counts( walk, entities )
    clusters.info = data.frame(
            serial = as.integer( names( clusters.sizes ) ),
            size = clusters.sizes,
            counts = clusters.counts[ names( clusters.sizes ) ],
            available = TRUE,
            stringsAsFactors = FALSE
        )
    rownames( clusters.info ) <- clusters.info$serial 
    cluster.serial.colname <- paste( entities, 'cluster.serial', sep = '.' )
    clusters.contents <- slot( walk, paste( entities, 'clusters', sep='.' ) )
    entity.colname <- ifelse(entities=='objects','object','state')
    clusters.info$score <- as.numeric( apply( clusters.info, 1, function (row ) score( as.integer(row[['size']]), as.integer(row[['counts']]) ) ) )

    entities.info <- ddply( clusters.contents, c( entity.colname ), function( shared.clusters ) {
        data.frame( nclusters = nrow( shared.clusters ), available = TRUE, stringsAsFactors = FALSE )
    } )
    rownames( entities.info ) <- entities.info[,entity.colname]
    #print(entities.info)

    res <- list()

    while ( any(entities.info$available ) ) {
        avail.clusters.info <- subset( clusters.info, available )
        if ( nrow( avail.clusters.info ) == 0 ) {
            warning( 'No valid clusters left, partition incomplete' )
            break
        }
        priority.entities <- as.character( subset( entities.info, available & nclusters == 1 )[,entity.colname] )
        if ( length(priority.entities) > 0 ) {
            priority.clusters <- unique( clusters.contents[ clusters.contents[,entity.colname] %in% priority.entities, cluster.serial.colname ] )
            #print( 'Priority entries:' )
            #print( priority.entities )
            #print( priority.clusters )
            avail.clusters.info <- subset( avail.clusters.info, serial %in% priority.clusters )
            if ( nrow( avail.clusters.info ) == 0 ) {
                warning( 'internal error, no available clusters' )
                break
            }
        }
        best.cluster.serial <- avail.clusters.info$serial[[ order( avail.clusters.info$score, decreasing = TRUE )[[1]] ]]
        #print( best.cluster.serial )
        res <- unlist( append( res, c( best.cluster.serial ) ) )
        cluster.elems <- as.character( clusters.contents[ clusters.contents[, cluster.serial.colname] == best.cluster.serial, entity.colname ] )
        #print( cluster.elems )
        sharing.clusters <- unique( clusters.contents[ clusters.contents[, entity.colname] %in% cluster.elems, cluster.serial.colname ] )
        sharing.clusters.elems <- unique( as.character( clusters.contents[ clusters.contents[, cluster.serial.colname] %in% sharing.clusters, entity.colname ] ) )
        cluster.elems.mask <- entities.info[,entity.colname] %in% cluster.elems
        sharing.clusters.elems.mask <- entities.info[,entity.colname] %in% sharing.clusters.elems
        sharing.clusters.mask <- clusters.info$serial %in% sharing.clusters
        if ( !all( entities.info[ cluster.elems.mask, 'available' ]) ) {
            warning('internal error, entities available flags')
            break
        }
        entities.info[ sharing.clusters.elems.mask, 'nclusters' ] <- entities.info[ sharing.clusters.elems.mask, 'nclusters' ] - 1
        entities.info[ cluster.elems.mask, 'available' ] <- FALSE
        clusters.info[ sharing.clusters.mask, 'available' ] <- FALSE
    }
    if ( length( res ) > 0 ) res <- sort( res )
    return ( res )
}

# prepare pivot table of signal means and SDs
BIMAP.signal_stats <- function( signals.frame )
{
    # prepare pivot table of signal means and SDs
    signals.pivot.mean <- daply( signals.frame,
        .(proteins.cluster, samples.cluster ),
        function( rows ) mean( rows$signal ) )
    signals.pivot.sd <- daply( signals.frame,
        .(proteins.cluster, samples.cluster ),
        function( rows ) sd( rows$signal ) )
    return ( c( list( signals.mean = signals.pivot.mean,
                      signals.sd = signals.pivot.sd ) ) )
}

#' Extracts chessboard biclustering with specified ID from random walk.
#' @param bimap.walk sampling walk
#' @param bimapId ID of chessboard biclustering to extract
#' @param onblock.threshold "on" blocks, are those whose "on" state frequency is greater than this threshold
#' @param extract.signals extract signals of individual cross-clusters
#' @returnType 
#' @return 
#' @author astukalov
BIMAP.extract_clustering <- function( bimap.walk, bimapId,
    onblock.threshold = 0.6,
    extract.signals = TRUE
){
    steps <- subset( bimap.walk@clusterings.walk, clustering.serial == bimapId )$step
    if ( length(steps) == 0 ) {
        stop( paste("No chessboard biclustering with ID=", bimapId, "found") )
    }
    
    # get object and state clusters of chessboard biclustering
    opId <- subset( bimap.walk@clusterings, clustering.serial == bimapId )$objects.partition.serial
    spId <- subset( bimap.walk@clusterings, clustering.serial == bimapId )$states.partition.serial
    ocIds <- as.character( subset( bimap.walk@objects.partitions, objects.partition.serial == opId )$objects.cluster.serial )
    scIds <- as.character( subset( bimap.walk@states.partitions, states.partition.serial == spId )$states.cluster.serial )

    if ( is.numeric(block.threshold) ) {
        # build consensus cross-clusters
        ccStats <- subset( bimap.walk@blocks.freq, objects.cluster.serial %in% ocIds & states.cluster.serial %in% scIds )
        crossClus <- subset( ccStats, enabled >= block.threshold * total,
                             select=c( 'objects.cluster.serial', 'states.cluster.serial' ) )
    } else {
        # get non-empty cross-clusters of the clustering
        crossClus <- subset( bimap.walk@blocks, clustering.serial == bimapId,
                             select = c("objects.cluster.serial", "states.cluster.serial") )
    }

    cluCounts <- nrow( crossClus )
    if ( cluCounts == 0 ) {
        stop( paste("No cross-clusters found in clustering ID=", bimapId ) )
    }
    colnames( crossClus ) <- c( 'proteins.cluster', 'samples.cluster' )
    crossClus$proteins.cluster <- as.character( crossClus$proteins.cluster )
    crossClus$samples.cluster <- as.character( crossClus$samples.cluster )

    proteins.clusters <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% ocIds,
        select = c( 'objects.cluster.serial', 'object' ) )
    colnames( proteins.clusters ) <- c( 'proteins.cluster', 'protein_ac' )
    proteins.clusters$proteins.cluster <- as.character( proteins.clusters$proteins.cluster )
    proteins.clusters$protein_ac <- as.character( proteins.clusters$protein_ac )
    rownames( proteins.clusters ) <- proteins.clusters$protein_ac
    samples.clusters <- subset( bimap.walk@states.clusters, states.cluster.serial %in% scIds,
        select = c( 'states.cluster.serial', 'state' ) )
    colnames( samples.clusters ) <- c( 'samples.cluster', 'sample' )
    samples.clusters$sample <- as.character( samples.clusters$sample )
    samples.clusters$samples.cluster <- as.character( samples.clusters$samples.cluster )
    rownames( samples.clusters ) <- samples.clusters$sample

    res <- list(
            steps = steps,
            cross.clusters = crossClus,
            proteins.clusters = proteins.clusters,
            proteins.clusters.info = bimap.walk@objects.clusters.info[ unique( proteins.clusters$proteins.cluster ),
                                                                    c( 'objects.cluster.serial', 'size',
                                                                       'nsteps', 'nsteps.included', 'avg.pairs.cooccur') ],
            samples.clusters = samples.clusters,
            samples.clusters.info = bimap.walk@states.clusters.info[ unique( samples.clusters$samples.cluster ),
                                                                 c( 'states.cluster.serial', 'size',
                                                                    'nsteps', 'nsteps.included', 'avg.pairs.cooccur') ]
    )
    colnames( res$proteins.clusters.info ) <- c( 'proteins.cluster', 'size', 'nsteps', 'nsteps_included', 'avg_pairs_cooccur' )
    res$proteins.clusters.info$avg_pairs_freq <- res$proteins.clusters.info$avg_pairs_cooccur / nrow( bimap.walk@clusterings.walk )
    colnames( res$samples.clusters.info ) <- c( 'samples.cluster', 'size', 'nsteps', 'nsteps_included', 'avg_pairs_cooccur' )
    res$samples.clusters.info$avg_pairs_freq <- res$samples.clusters.info$avg_pairs_cooccur / nrow( bimap.walk@clusterings.walk )
    if ( extract.signals ) {
        print( paste( 'Extracting signals for clustering #', bimapId, '...', sep='') )
        # get subframe with all samples for signals of given bimap
        # (but samples might be from other bimaps, which contain the same clusters)
        cc_str <- paste( crossClus$proteins.cluster, crossClus$samples.cluster )
        signals.subframe <- subset( bimap.walk@signals,
                paste( objects.cluster.serial, states.cluster.serial ) %in% cc_str )[c('step', 'objects.cluster.serial', 'states.cluster.serial', 'signal')]
        colnames( signals.subframe ) <- c( 'step', 'proteins.cluster', 'samples.cluster', 'signal' )
        signals.subframe$proteins.cluster <- as.character( signals.subframe$proteins.cluster )
        signals.subframe$samples.cluster <- as.character( signals.subframe$samples.cluster )
        signal_stats <- BIMAP.signal_stats( signals.subframe )
        res$cross.clusters$signal.mean <- apply( res$cross.clusters, 1,
                function( cc ) signal_stats$signals.mean[ cc[['proteins.cluster']], cc[['samples.cluster']] ] )
        res$cross.clusters$signal.sd = apply( res$cross.clusters, 1,
                function( cc ) signal_stats$signals.sd[ cc[['proteins.cluster']], cc[['samples.cluster']] ] )
        res <- c( res, signal_stats )
        res$signals.subframe = signals.subframe
        print( 'Signals extracted' )
    }
    return ( res )
}

BIMAP.filter_clustering <- function( bimap.props, protein_acs, sample_acs )
{
    # filter elements and clusters
    res <- list( #proteins = subset( bimap.props$proteins, protein_ac %in% protein_acs ),
                 #samples = subset( bimap.props$samples, sample %in% sample_acs ),
                 proteins.clusters = subset( bimap.props$proteins.clusters, protein_ac %in% protein_acs ),
                 samples.clusters = subset( bimap.props$samples.clusters, sample %in% sample_acs ) )
    samples.cluster_acs <- unique( res$samples.clusters$samples.cluster )
    proteins.cluster_acs <- unique( res$proteins.clusters$proteins.cluster )
    # remove cross clusters corresponding to removed clusters
    res$cross.clusters <- subset( bimap.props$cross.clusters,
                                  proteins.cluster %in% proteins.cluster_acs &
                                  samples.cluster %in% samples.cluster_acs )
    res$signals.subframe <- subset( bimap.props$signals.subframe,
                                  proteins.cluster %in% proteins.cluster_acs &
                                  samples.cluster %in% samples.cluster_acs )
    signal_stats <- BIMAP.signal_stats( res$signals.subframe )
    for ( i in 1:nrow(res$cross.clusters) ) {
        ccrow <- res$cross.clusters[ i, ] 
        res$cross.clusters$signal.mean = signal_stats$signals.mean[ ccrow$proteins.cluster,
                                                                   ccrow$samples.cluster ] 
        res$cross.clusters$signal.sd = signal_stats$signals.sd[ ccrow$proteins.cluster,
                                                               ccrow$samples.cluster ]
    }
    res$signals.mean <- signal_stats$signals.mean
    res$signals.sd <- signal_stats$signals.sd
    return ( res )
}

BIMAP.objects.cluster.labels <- function( walk )
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

BIMAP.objects.partition.labels <- function( walk )
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

BIMAP.states.cluster.labels <- function( walk )
{
    df <- data.frame( serial = c(), label = c() )
    for ( clustering in walk@clusterings ) {
        if ( length( clustering@samples.clusters$states.cluster.serial ) == 0 ) next;
        ocs <- as.data.frame( clustering@samples.clusters )
        df <- unique( rbind( df, ddply( ocs, .(states.cluster.serial), function( statesclu ) {
                states <- statesclu[,'state']
                states <- states[ order(states ) ]
                data.frame( serial = statesclu[[1,'states.cluster.serial']], label = paste( states, collapse = ',' ) )
            } ) ) )
    }
    res <- as.character( df$label )
    names( res ) <- as.integer( df$serial ) 
    return ( res )
}

BIMAP.states.partition.labels <- function( walk )
{
    res <- list()
    for ( clustering in walk@clusterings ) {
        if ( length( clustering@samples.clusters$states.cluster.serial ) == 0 ) next;
        ocs <- as.data.frame( clustering@samples.clusters )
        clunames <- as.character( ddply( ocs, .(states.cluster.serial), function( stateclu ) {
                    states <- stateclu[,'state']
                    states <- states[ order(states) ]
                    data.frame( serial = stateclu[[1,'states.cluster.serial']], label = paste( states, collapse = '' ) )
                } )$label )
        clunames <- clunames[ order( clunames ) ]
        res[[ as.character(clustering@states.partition.serial) ]] = paste( clunames, collapse = ' ' )
    }
    return ( res )
}

BIMAP.interactors.table <- function( bimap.walk, protein.ac, protein.info =-NULL, min.freq = 0.01 )
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
BIMAP.extract.stable.clusters <- function( bimap.walk,
        ms_data, sample_col = 'sample', prey_col = 'prey_ac',
        min.avg.nsample = 3, min.size = 2,
        min.included.freq = 0.95, size.weight = 0.2 )
{
    min.nsteps.included = nrow( bimap.walk@clusterings.walk ) * min.included.freq
    good.clusters.info <- subset( bimap.walk@objects.clusters.info, size >= min.size & nsteps.included >= min.nsteps.included )
    good.clusters.info$score <- ( good.clusters.info$nsteps.included / nrow( bimap.walk@clusterings.walk ) - min.included.freq ) / ( 1.0 - min.included.freq ) +
                  size.weight * ( good.clusters.info$size - min.size )
    #print( good.clusters.info )
    good.clusters <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% good.clusters.info$objects.cluster.serial )
    avg.nsamples <- ddply( unique( merge( good.clusters, ms_data, by.x = 'object', by.y = prey_col )
                            [ c( 'objects.cluster.serial', 'object', sample_col ) ] ), .( objects.cluster.serial ), function( clu_ms_data ) {
                data.frame( avg.nsample = min( table( clu_ms_data$object ) ),
                        stringsAsFactors = TRUE )
            } )
    good.clusters.info <- merge( good.clusters.info, avg.nsamples, by = c( 'objects.cluster.serial' ) )
    good.clusters.info <- subset( good.clusters.info, avg.nsample >= min.avg.nsample )
    good.clusters <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% good.clusters.info$objects.cluster.serial )
    intersecting.pairs <- unique( merge( good.clusters, good.clusters, by = c( 'object' ),
                    all = FALSE, sort = FALSE )
                    [ c( 'objects.cluster.serial.x', 'objects.cluster.serial.y' ) ] )
    best.clusters.serials <- c()
    while ( nrow( good.clusters.info ) > 0 ) {
        best.cluster <- good.clusters.info[ order( good.clusters.info$score, decreasing = TRUE )[1], 'objects.cluster.serial' ]
        best.clusters.serials <- c( best.clusters.serials, best.cluster )
        intersects.with <- subset( intersecting.pairs, objects.cluster.serial.x == best.cluster )$objects.cluster.serial.y
        good.clusters.info <- subset( good.clusters.info, !( objects.cluster.serial %in% intersects.with ) )
    }
    res <- subset( bimap.walk@objects.clusters, objects.cluster.serial %in% best.clusters.serials )
    res <- res[ order( res$objects.cluster.serial, res$object ), ]
    return ( res )
}
