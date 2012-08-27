source( file.path( bimap_scripts_path, "BIMAP.R" ) )
source( file.path( bimap_scripts_path, "BIMAP_math.R" ) )

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
BIMAP.mcmcwalk.biclusterings_order <- function( walk, n = 30, by = 'parts.pdf' ) {
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
#' @param min.avg.noccur average number of samples the cluster is seen
#' @param min.size minimal size of cluster
#' @param min.included.freq minimal frequency the cluster is seen (itself or as subcomponent) in partitions
#' @param size.weight weight of cluster in scoring
#' @returnType dataframe of clusters
#' @return 
#' @author astukalov
#' @export
BIMAP.mcmcwalk.extract_stable_clusters <- function( bimap.walk, bimap.data,
		min.avg.noccur = 3, min.size = 2,
		min.included.freq = 0.95, size.weight = 0.2,
		allow.intersections = FALSE,
		use.object.clusters = TRUE,
		greedy.algorithm = FALSE )
{
	elem_col <- ifelse( use.object.clusters, 'object', 'probe' )
	cluster_col <- ifelse( use.object.clusters, 'objects.cluster.serial', 'probes.cluster.serial' )
	
	clusters.info <- if ( use.object.clusters ) bimap.walk@objects.clusters.info
			else bimap.walk@probes.clusters.info
	clusters <- if ( use.object.clusters ) bimap.walk@objects.clusters
			else bimap.walk@probes.clusters
	clusters$element <- clusters[,elem_col]
	clusters$cluster.serial <- clusters[,cluster_col]
	clusters.info$cluster.serial <- clusters.info[,cluster_col]
	min.nsteps.included = nrow( bimap.walk@clusterings.walk ) * min.included.freq
	good.clusters.info <- subset( clusters.info, size >= min.size & nsteps.included >= min.nsteps.included )
	good.clusters.info$score <- ( good.clusters.info$nsteps.included / nrow( bimap.walk@clusterings.walk ) ) *
			( good.clusters.info$size / max( good.clusters.info$size ) )^size.weight
	#print( good.clusters.info )
	good.clusters <- subset( clusters, cluster.serial %in% good.clusters.info$cluster.serial )
	avg.noccur <- ddply( unique( merge( good.clusters, merge( bimap.data$measurements, bimap.data$exp_design ),
							by.x = 'element',
							by.y = ifelse( use.object.clusters, 'prey_ac', 'sample' ) )
							[ c( 'cluster.serial', 'element', ifelse( use.object.clusters, 'sample', 'prey_ac' ) ) ] ),
			.( cluster.serial ), function( clu_ms_data ) {
				# average number of samples (preys for probes) where object is seen
				data.frame( avg.noccur = min( table( clu_ms_data$element ) ),
						stringsAsFactors = FALSE )
			} )
	good.clusters.info <- merge( good.clusters.info, avg.noccur, by = c( 'cluster.serial' ) )
	good.clusters.info <- subset( good.clusters.info, avg.noccur >= min.avg.noccur )
	good.clusters <- subset( clusters, cluster.serial %in% good.clusters.info$cluster.serial )
	message( 'Total ', length( unique( good.clusters$element ) ), ' element(s) in ',
			length( unique( good.clusters$cluster ) ), ' stable cluster(s)' )
	if ( !allow.intersections ) {
		if ( greedy.algorithm ) {
			intersecting.pairs <- unique( merge( good.clusters, good.clusters, by = c( 'element' ),
							all = FALSE, sort = FALSE )
							[ c( 'cluster.serial.x', 'cluster.serial.y' ) ] )
			best.clusters.serials <- c()
			while ( nrow( good.clusters.info ) > 0 ) {
				best.cluster <- good.clusters.info[ order( good.clusters.info$score, decreasing = TRUE ), 'cluster.serial' ][1]
				best.clusters.serials <- c( best.clusters.serials, best.cluster )
				intersects.with <- subset( intersecting.pairs, cluster.serial.x == best.cluster )$cluster.serial.y
				good.clusters.info <- subset( good.clusters.info, !( cluster.serial %in% intersects.with ) )
			}
			res <- subset( clusters, cluster.serial %in% best.clusters.serials )
		} else {
			res <- BIMAP.optimal_partition( good.clusters, good.clusters.info, 'cluster.serial', 'element', 'score' )
			res[,elem_col] <- res$element
			res[,cluster_col] <- res$cluster.serial
		}
	} else {
		res <- subset( clusters, cluster.serial %in% good.clusters.info$cluster.serial )
	}
	message( length( unique( res$element) ), ' ', elem_col, "(s) of ", length( unique( clusters$element ) ),
			" in ", length(unique(res$cluster.serial)), " stable cluster(s)" )
	res <- res[ order( res$cluster.serial, res$element ), ]
	return ( res )
}

#' Extracts chessboard biclustering with specified ID from random walk.
#' @param bimap.walk sampling walk
#' @param bimapId ID of chessboard biclustering to extract
#' @param onblock.threshold "on" blocks, are those whose "on" probe frequency is greater than this threshold
#' @returnType 
#' @return 
#' @author astukalov
BIMAP.mcmcwalk.extract_stable_biclustering <- function(
		bimap.walk, bimap.data,
		onblock.threshold = 0.6,
		min.avg.nprobes = 3, min.avg.nobjects = 3,
		min.objects.freq = 0.95, min.probes.freq = 0.95,
		objects.size.weight = 0.2, probes.size.weight = 0.2,
		greedy.algorithm = FALSE
){
	if ( nrow( bimap.walk@stable.blocks.scores ) == 0 ) {
		error( 'No stable blocks information found in the BI-MAP walk, run BIMAP.mcmcload() again with stable thresholds defined' )
	}
	proteins.clusters <- BIMAP.mcmcwalk.extract_stable_clusters(
			bimap.walk, bimap.data,
			min.avg.noccur = min.avg.nprobes, min.size = 1,
			min.included.freq = min.objects.freq,
			size.weight = objects.size.weight,
			allow.intersections = FALSE, use.object.clusters = TRUE,
			greedy.algorithm = greedy.algorithm )[,c('objects.cluster.serial','object')]
	colnames( proteins.clusters ) <- c( 'proteins.cluster', 'protein_ac' )
	samples.clusters <- BIMAP.mcmcwalk.extract_stable_clusters(
			bimap.walk, bimap.data,
			min.avg.noccur = min.avg.nobjects, min.size = 1,
			min.included.freq = min.probes.freq,
			size.weight = probes.size.weight,
			allow.intersections = FALSE, use.object.clusters = FALSE,
			greedy.algorithm = greedy.algorithm )[,c('probes.cluster.serial','probe')]
	colnames( samples.clusters ) <- c( 'samples.cluster', 'sample' )
	
	if ( is.numeric(onblock.threshold) ) {
		# build consensus blocks
		blocks <- subset( bimap.walk@stable.blocks.scores,
				enabled >= onblock.threshold * total &
						objects.cluster.serial %in% proteins.clusters$proteins.cluster
						& probes.cluster.serial %in% samples.clusters$samples.cluster,
				select=c( 'objects.cluster.serial', 'probes.cluster.serial', 'signal.mean', 'signal.var' ) )
	} else {
		# get non-empty blocks of the clustering
		blocks <- subset( bimap.walk@stable.blocks.scores,
				objects.cluster.serial %in% proteins.clusters$proteins.cluster
						& probes.cluster.serial %in% samples.clusters$samples.cluster,
				select = c("objects.cluster.serial", "probes.cluster.serial", 'signal.mean', 'signal.var' ) )
	}
	
	nBlocks <- nrow( blocks )
	if ( nBlocks == 0 ) {
		stop( "No on-blocks found in stable clustering" )
	} else {
		message( nBlocks, ' on-blocks found of ', length( unique( proteins.clusters$proteins.cluster ) ), 
				'x', length( unique( samples.clusters$samples.cluster ) ), ' possible' )
	}
	
	colnames( blocks ) <- c( 'proteins.cluster', 'samples.cluster', 'signal.mean', 'signal.sd' )
	blocks$proteins.cluster <- as.character( blocks$proteins.cluster )
	blocks$samples.cluster <- as.character( blocks$samples.cluster )
	
	proteins.clusters$proteins.cluster <- as.character( proteins.clusters$proteins.cluster )
	proteins.clusters$protein_ac <- as.character( proteins.clusters$protein_ac )
	rownames( proteins.clusters ) <- proteins.clusters$protein_ac
	
	samples.clusters$sample <- as.character( samples.clusters$sample )
	samples.clusters$samples.cluster <- as.character( samples.clusters$samples.cluster )
	rownames( samples.clusters ) <- samples.clusters$sample
	
	res <- list(
			blocks = blocks,
			proteins.clusters = proteins.clusters,
			proteins.clusters.info = bimap.walk@objects.clusters.info[ unique( proteins.clusters$proteins.cluster ),
					c( 'objects.cluster.serial', 'size',
							'nsteps', 'nsteps.included', 'avg.pairs.cooccur') ],
			samples.clusters = samples.clusters,
			samples.clusters.info = bimap.walk@probes.clusters.info[ unique( samples.clusters$samples.cluster ),
					c( 'probes.cluster.serial', 'size',
							'nsteps', 'nsteps.included', 'avg.pairs.cooccur') ],
			signals.mean = daply( blocks, c( 'proteins.cluster', 'samples.cluster' ), function ( data ) data$signal.mean ),
			signals.sd = daply( blocks, c( 'proteins.cluster', 'samples.cluster' ), function ( data ) data$signal.sd )
	)
	colnames( res$proteins.clusters.info ) <- c( 'proteins.cluster', 'size', 'nsteps', 'nsteps_included', 'avg_pairs_cooccur' )
	res$proteins.clusters.info$avg_pairs_freq <- res$proteins.clusters.info$avg_pairs_cooccur / nrow( bimap.walk@clusterings.walk )
	colnames( res$samples.clusters.info ) <- c( 'samples.cluster', 'size', 'nsteps', 'nsteps_included', 'avg_pairs_cooccur' )
	res$samples.clusters.info$avg_pairs_freq <- res$samples.clusters.info$avg_pairs_cooccur / nrow( bimap.walk@clusterings.walk )
	
	return ( res )
}
