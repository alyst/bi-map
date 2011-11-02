# Mass-Spectrometry data utilities
# 
# Author: Alexey Stukalov (astukalov@cemm.oeaw.ac.at)
###############################################################################

#' Get a list of preys ms_data, shared by given baits
#' @param ms_data mass-spec data
#' @param baits list of bait symbols
#' @returnType data.frame
#' @return data frame with records corresponding to preys shared by all baits, ms_data for each bait in column prefixed by bait's name
#' @author astukalov
#' @export
sharedPreys <- function( ms_data, baits, by = c( 'prey_ac', 'condition' ), all = FALSE )
{
    merge_res <- NULL
    for ( bait in unique( baits ) ) {
        bait_ms <- subset( ms_data, bait_ac == bait )
        colnames( bait_ms ) <- lapply( colnames( bait_ms ), 
            function( colname ) ifelse( colname %in% by, colname, paste( colname, bait, sep='.' ) )
        )
        if ( is.null( merge_res ) ) {
            merge_res <- bait_ms
        }
        else {
            merge_res <- merge( merge_res, bait_ms, by = by, all = all )
        }
    }
    # calculate count of predicted bait-prey interactions for each prey
    if ( 'iaction' %in% colnames( ms_data ) ) {
        iact_cnt_expr <- parse( text = paste( lapply( baits, 
                    function( bait ) paste( "ifelse(merge_res$iaction.", bait, ",1,0)", sep="" ) ), collapse="+" ) )
        merge_res$iactions_count <- eval( iact_cnt_expr )
    }
    
    return ( merge_res )
}

sharedPreyCounts <- function( ms_data, tuple_len = 2 )
{
    bait_prey_data <- subset( ms_data, select = c( prey_ac, bait_ac ) )
    all_baits <- unique( ms_data$bait_ac )
    cross_data <- NULL
    for ( i in 1:tuple_len ) {
        colnames( bait_prey_data ) <- c( 'prey_ac', paste( 'bait_ac', i, sep='.' ) )
        if ( is.null( cross_data ) ) {
            cross_data <- bait_prey_data
        }
        else {
            cross_data <- merge( cross_data, bait_prey_data, by = 'prey_ac' )
        }
    }
    agg_data <- aggregate( cross_data, 
        by = lapply( 1:tuple_len, function( i ) eval( parse( text = paste( 'bait_ac.', i , '=cross_data$bait_ac.', i, sep='' ) ) ) ), 
        FUN = length )
    # filter for duplicate parts
    subset_expr <- parse( text = paste("c(",
            paste( lapply(1:tuple_len, function( i ) paste( "row[['Group.", i, "']]", sep="" ) ), collapse="," ), ")") )
    agg_data <- agg_data[ apply( agg_data, 1, function( row ) anyDuplicated( eval( subset_expr ) ) == 0 ), ]
    return ( agg_data[ order( agg_data$prey_ac ), ] )
}

#' Protein abundance function
#' @param sc total spectral counts
#' @param pc unique peptides counts
#' @param seqcov sequence coverage
#' @param sl sequence length
#' @returnType number
#' @return protein abundance
#' @author astukalov
#' @export
protein.log_abundance <- function( sc, pc, seqcov, sl ) {
    return ( ifelse( sc > 0, log(sc)+4-0.3*log(sl), 0 ) )# ( sqrt(sl )), 0 ) )
}

prey_observations <- function( ms_data )
{
    prey_observations <- aggregate( data.frame( prey_ac = ms_data$prey_ac, bait_ac = ms_data$bait_ac ),
        list( prey_ac = ms_data$prey_ac ), 
        FUN = function( rows ) as.numeric( length( unique(rows) ) ) )
    prey_obs_list <- prey_observations$bait_ac
    names( prey_obs_list ) <- prey_observations$prey_ac
    return ( prey_obs_list )
}

# converts protein accession to protein label
smart_prot_id <- function( prot_acs, prot_info )
{
    return ( sapply( prot_acs, function( prot_ac ) {
                if ( length( prot_ac ) == 0 | is.null( prot_ac ) | is.na( prot_ac) ) {
                    return ( '' )
                } else if ( prot_ac == '#bait#' ) {
                    return ( prot_ac )
                }
                else if ( prot_ac %in% c( '136429', '162648' ) ) {
                    return ( unlist(strsplit( prot_info[[ as.character(prot_ac), 'description' ]], ' ' ))[[1]] )
                }
                else {
                    return ( sub( '_[A-Z]+', '', prot_info[ as.character(prot_ac), 'protein_id' ] ) )
                } } ) )
}

#' Writes table into CSV file,
#' inserting separator character at the beginning of the column names line
ms_table_write <- function( table, file = NA )
{
    table2 <- cbind( matrix( rownames( table ), ncol = 1, nrow = nrow( table ) ), table )
    write.table( table2, file = file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE )
}

ms_table_write_with_info <- function( table, col2bait, file = NA )
{
    bait_names <- matrix( col2bait[colnames(table)], nrow = 1, ncol = ncol( table ) )
    colnames(bait_names) <- colnames(table)
    rownames(bait_names) <- 'BAIT'
    write.table( bait_names, file = file, sep = '\t', quote = FALSE, row.names = TRUE, col.names = TRUE )
    write.table( table, file = file, sep = '\t', quote = FALSE, row.names = TRUE, col.names = FALSE, append = TRUE )
}

#' Converts MS data.frame into matrix (rows=proteins, cols=msruns)
#' by applying cellfunc to each DS row
ms_data.matrix <- function( ms_data, protein_info, cellfunc, msrun_column = 'msrun',
                            bait_column = 'bait_ac', prey_column = 'prey_ac',
                            sc_column = 'sc', pc_column = 'pc', seqcov_column = 'seqcov'
){
    ms_data.x <- ms_data
    ms_data.x$cell <- cellfunc( ms_data.x[,sc_column], ms_data.x[,pc_column], ms_data.x[,seqcov_column],
        protein_info[ as.character( ms_data.x[,prey_column] ), 'seqlength' ] )
    ms_data.wide <- reshape( ms_data.x[,unique(c(msrun_column,prey_column,'cell'))], 
        timevar = c( msrun_column ), idvar = c(prey_column), 
        ids = unique( as.character( ms_data[,prey_column] ) ),
        times = unique( as.character( ms_data[,msrun_column] ) ), 
        v.names = c('cell'), direction = 'wide', sep = '.' )
    ms_data.mtx <- as.matrix( ms_data.wide[ , setdiff( colnames( ms_data.wide ), c(prey_column) ) ] )
    rownames( ms_data.mtx ) <- as.character( ms_data.wide[,prey_column] )
    colnames( ms_data.mtx ) <- sub( 'cell.', '', colnames( ms_data.mtx ), fixed = TRUE )
    return ( ms_data.mtx )
}

#' Shuffle MS data
#' @param ms_data MS data.frame
#' @param shuffle.proteins whether proteins should be shuffled between samples
#' @param shuffle.samples whether samples where protein is seen should be shuffled 
#' @returnType data.frame
#' @return MS dataframe with shuffled data
#' @author astukalov
#' @export
ms_data.shuffle <- function( ms_data,
                             msrun_col = 'msrun', sample_col = 'sample',
                             bait_ac_col = 'bait_ac', prey_ac_col = 'prey_ac',
                             shuffle.proteins.rate = 1.0, shuffle.samples.rate = 1.0 )
{
    col.map <- data.frame( name = c( msrun_col,sample_col,bait_ac_col ),
                           label = c( 'msrun', 'sample', 'bait_ac' ),
                           stringsAsFactors = FALSE )
    col.map$ix <- sapply( col.map$name, function(cn) {
                           which(cn == colnames(ms_data)) } )
    rownames( col.map ) <- col.map$label
    col.map.unique <- unique( col.map[,c('name','ix')] )
    col.map.unique$label <- rownames( col.map.unique )
    msruns <- unique( ms_data[, col.map$ix ] )
    colnames( msruns ) <- col.map$label
    rownames( msruns ) <- msruns$msrun
    all.samples <- unique( msruns$sample )
    all.preys <- unique( ms_data[,prey_ac_col] )
    samples.replicates <- table( msruns$sample )
    ms_data.new <- ms_data
    if ( shuffle.samples.rate > 0.0 ) {
        ms_data.new <- ddply( ms_data.new, c(prey_ac_col), function( prey_ms_data ) {
            samples <- unique( prey_ms_data[,sample_col] )
            nsamples <- length( samples )
            nsamples.shuff <- floor( nsamples * shuffle.samples.rate )
            new.samples <- sample( samples, nsamples - nsamples.shuff )
            new.msruns <- subset( msruns, sample %in% new.samples )$msrun
            replaced.samples <- setdiff( samples, new.samples )
            old.samples <- c( new.samples, replaced.samples )
            old.msruns <- new.msruns
            avail.samples <- setdiff( all.samples, new.samples )
            for ( old.sample in replaced.samples ) {
                # select new sample that would have the same number of replicates
                nrepl <- samples.replicates[ old.sample ]
                alt.ves <- intersect( avail.samples, names( samples.replicates[ samples.replicates == nrepl ] ) )
                if ( length( alt.ves ) > 0 ) {
                    new.sample <- sample( alt.ves, 1 )
                } else {
                    new.sample <- old.sample
                }
                avail.samples <- setdiff( avail.samples, new.sample )
                new.samples <- c( new.samples, new.sample )
                new.msruns <- c( new.msruns, subset( msruns, sample == new.sample )$msrun )
                old.msruns <- c( old.msruns, subset( msruns, sample == old.sample )$msrun )
            }
            names( new.msruns ) <- old.msruns
            prey_ms_data[, col.map.unique$name ] <- msruns[ new.msruns[ prey_ms_data[,msrun_col] ],
                                                            col.map.unique$label ]
            prey_ms_data[,prey_ac_col] <- NULL
            return ( prey_ms_data )
        } )
    }
    if ( shuffle.proteins.rate > 0.0 ) {
        ms_data.new <- ddply( ms_data.new, c(sample_col), function( sample_ms_data ) {
            preys <- unique( sample_ms_data[,prey_ac_col] )
            npreys <- length( preys )
            npreys.shuff <- floor( npreys * shuffle.proteins.rate )
            new.preys <- sample( preys, npreys - npreys.shuff ) 
            new.preys <- c( new.preys, sample( setdiff( all.preys, new.preys ), npreys.shuff, replace=FALSE,) )
            names( new.preys ) <- preys
            # replace preys with new ones
            sample_ms_data[,prey_ac_col] <- new.preys[ sample_ms_data[,prey_ac_col] ]
            sample_ms_data[,sample_col] <- NULL
            return ( sample_ms_data )
        } )
    }
    return ( ms_data.new[ , colnames( ms_data ) ] )
}

