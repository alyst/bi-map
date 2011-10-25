# Mass-Spectrometry data utilities
# 
# Author: Alexey Stukalov (astukalov@cemm.oeaw.ac.at)
###############################################################################

# load RPSIbbb library
RPeptideTools.libpath <- file.path( projects_path, "pcp/peptideTools/build/release" )
#dyn.unload( file.path( RPeptideTools.libpath, paste("librpeptidetools", .Platform$dynlib.ext, sep="")) ) 
dyn.load( file.path( RPeptideTools.libpath, paste("librpeptidetools", .Platform$dynlib.ext, sep="")), 
          type = "Call" )

DigestProtein <- function( protein,
    min.peptide.length = 6,
    max.peptide.length = 100,
    max.cleavages = 0
){
    return ( .Call( "DigestProtein", protein,
                    list( min.peptide.length = min.peptide.length,
                          max.peptide.length = max.peptide.length,
                          max.cleavages = max.cleavages ) ) )
}

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

sharedPreysDubSQL <- function( bait1Symb, bait2Symb )
{
    res <- dbGetQuery( dub.dbConn, paste( "
                SELECT im.molecule_id, im.symbol, 
                avg( array_avg( md1.unique_peptides ) ) AS up1, 
                avg( array_avg( md2.unique_peptides ) ) AS up2,
                avg( md1.total_spectral_counts ) AS tsc1,
                avg( md2.total_spectral_counts ) AS tsc2,
                im.seqlength as seqlength, im.molweight, 
                EXISTS( select * from interaction i WHERE i.bait_id = b1.bait_id AND i.interactor_molecule_id = im.molecule_id ) AS exists1,
                EXISTS( select * from interaction i WHERE i.bait_id = b2.bait_id AND i.interactor_molecule_id = im.molecule_id ) AS exists2
                FROM molecule AS bm1
                JOIN bait AS b1 ON bm1.molecule_id = b1.molecule_id
                JOIN ms_data AS md1 ON md1.bait_id = b1.bait_id
                JOIN ms_data AS md2 ON md2.interactor_molecule_id = md1.interactor_molecule_id
                JOIN bait AS b2 ON md2.bait_id = b2.bait_id
                JOIN molecule AS bm2 ON bm2.molecule_id = b2.molecule_id
                JOIN molecule AS im ON im.molecule_id = md1.interactor_molecule_id
                WHERE bm1.symbol = '", bait1Symb, "' AND bm2.symbol = '", bait2Symb, "'
                AND im.entryid is NOT NULL 
                GROUP BY im.molecule_id, im.symbol, b1.bait_id, b2.bait_id, im.seqlength, im.molweight",
            sep = '' ) )
    res$bait1 <- bait1Symb
    res$bait2 <- bait2Symb
    return ( res )
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

uniprotid_to_name <- function( idlist, protein_info ) {
    lapply( idlist, function( up_id ) {
            res <- match( up_id, protein_info$primaryac )
            if ( is.na( res ) ) return ( up_id )
            protein_name <- protein_info[res, 'protein_id' ]
            return ( ifelse( !is.na(protein_name), protein_name, up_id ) )
        }
    )
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

ms_data.baits <- function( ms_data )
{
    return( unique( as.character( ms_data$bait_ac ) ) )
}

# ms data subset for bait-only ms_data, where prey_ac (being bait itself) is replaced by '#bait#' keyword
ms_data.chimeric_bait <- function( ms_data )
{
    bait_ms_data <- subset( ms_data, as.character( prey_ac ) == as.character( bait_ac ) )
    bait_ms_data$prey_ac.orig <- bait_ms_data$prey_ac # original prey symbol for chimeric bait
    bait_ms_data$prey_ac <- rep( '#bait#', nrow( bait_ms_data ) )
    return ( bait_ms_data )
}

noglob <- function( msdata ) subset( msdata, msrun != '#glob#' )

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


#' Filters isoforms of prey proteins
#' across all samples, and return ones
#' that are used for bait or isoforms detected more
#' than other isoforms of the same protein.
ms_data.isoforms.priority <- function( ms_data, specificity = c( NULL, 'sup', 'sp' ) )
{
    ms_data_noglob <- subset( ms_data, msrun != '#glob#' )
    baits_ac <- unique( as.character( ms_data_noglob$bait_ac ) )
    # map from AC to AC without isoform
    preys_ac <- sort( unique( as.character( ms_data_noglob$prey_ac ) ) )
    preys_ac_non_iso <- sapply( strsplit( preys_ac, '-', fixed = TRUE ), 
                                function( res ) res[[1]] )
    names( preys_ac_non_iso ) <- preys_ac
    #print( preys_ac_non_iso )

    # collect statistics for each isoform
    prey_stats <- ddply( ms_data_noglob, .(prey_ac), function( prey_rows ) {
            prey_ac <- prey_rows$prey_ac[1]
            return ( data.frame( prey_ac = prey_ac,
                                 prey_ac_non_iso = preys_ac_non_iso[ prey_ac ],
                                 not.bait = !any( prey_ac %in% prey_rows$bait_ac ),
                                 n.detections = nrow( prey_rows ),
                                 pc_sup.max = max( 0, prey_rows$pc_sup, na.rm = TRUE ),
                                 pc_sp.max = max( 0, prey_rows$pc_sp, na.rm = TRUE ), 
                                 pc.max = max( 0, prey_rows$pc, na.rm = TRUE,
                                 stringsAsFactors = FALSE ) ) )
        } )
    #print( prey_stats )

    # filter those isoforms, that have either:
    # 1) the same as bait of this pulldown
    # 2) most superspecific peptides
    # 3) most specific peptides
    # 4) most detections
    # 5) most peptides
    # 6) smallest isoform ID
    res <- ddply( prey_stats, .(prey_ac_non_iso), function( iso_stats ) {
            iso_stats <- iso_stats[ order( iso_stats$not.bait,
                -iso_stats$pc_sup.max,
                -iso_stats$pc_sp.max,
                -iso_stats$n.detections,
                -iso_stats$pc.max,
                iso_stats$prey_ac ), ]
            return ( data.frame(
                    prey_ac = iso_stats$prey_ac,
                    pc_sup.max = iso_stats$pc_sup.max,
                    pc_sp.max = iso_stats$pc_sp.max,
                    priority = 1:nrow(iso_stats),
                    stringsAsFactors = FALSE ) )
        } )
    if ( specificity == 'sup' ) {
        # isoform needs to have at least one superspecific peptide
        res <- subset( res, pc_sup.max > 0 )
    }
    else if ( specificity == 'sp' ) {
        # isoform needs to have at least one specific peptide
        res <- subset( res, pc_sp.max > 0 )
    }
    rownames( res ) <- res$prey_ac
    return( res )
}

ms_data.interaction_isoforms.priority <- function( ms_data, prey.isoforms.priority )
{
    # collect statistics for each isoform
    print( paste("Before isoforms prioritizing: ", nrow(ms_data) ))
    ms_data.prioritized <- ddply( ms_data, .(bait_ac, prey_ac_noiso), function( iact_rows ) {
            preys <- unique( iact_rows$prey_ac )
            preys.priority <- prey.isoforms.priority[ preys, 'priority' ]
            names( preys.priority ) <- preys
            preys.priority <- sort( preys.priority )
            preys <- names( preys.priority )
            # renumerate
            preys.priority <- 1:length( preys.priority )
            names( preys.priority ) <- preys
            iact_rows$prey_iso.priority <- preys.priority[ iact_rows$prey_ac ]
            iact_rows <- iact_rows[ order( iact_rows$prey_iso.priority, iact_rows$prey_ac,
                                    na.last = TRUE ), ]
            return ( iact_rows )
        } )
    print( paste("After isoforms prioritizing: ", nrow(ms_data.prioritized) ))
    return( ms_data.prioritized )
}

#' function appends #glob# entries to the dataset filtered by sample code
#' This is not 100% correct, since to #glob# entries, that were filtered out might also contribute 
ms_data.append_glob <- function( ms_data_filtered, ms_data )
{
    rbind( subset( ms_data_filtered, sample != '#glob#' ),
        subset( ms_data, sample == '#glob#' & paste( bait_ac, prey_ac ) %in% paste( ms_data_filtered$bait_ac, ms_data_filtered$prey_ac ) ) )
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

#' Calculate hypergeometric p-values
#' of co-occurence in the samples of preys in MS dataset
#' against the set of preselected preys
#' @param ms_data MS data (prey protein x sample)
#' @param core.preys set of ACs of prey proteins to compare with
#' @param sample_col name of MS sample column in ms_data
#' @param prey_col name of prey AC column in ms_data
#' @returnType data.frame
#' @return data.frame: prey, core.prey, p-value
#' @author Alexey Stukalov
#' @export
ms_data.preys.hyper_pvalue <- function(
        ms_data, core.preys,
        sample_col = 'sample', prey_col = 'prey_ac'
){
    core_preys.df <- ddply( ms_data[ ms_data[,prey_col] %in% core.preys, ], c( prey_col ), function( prey_data ) {
        return ( data.frame( sample = unique( prey_data[, sample_col ]) ) )
    } )
    colnames( core_preys.df ) <- c( 'prey', 'sample' )
    ncore.samples <- table( core_preys.df$prey )
    other_preys.df <- ddply( ms_data[ !(ms_data[,prey_col] %in% core.preys), ], c( prey_col ), function( prey_data ) {
        return ( data.frame( sample = unique( prey_data[, sample_col ]) ) )
    } )
    print(nrow(other_preys.df))
    colnames( other_preys.df ) <- c( 'prey', 'sample' )
    nother.samples <- table( other_preys.df$prey )
    all.samples <- unique( ms_data[,sample_col] )
    nsamples <- length( all.samples )
    prey.distances <- ddply( merge( core_preys.df, other_preys.df,
                                    by = c( 'sample' ),
                  suffixes = c( '.core', '.other' ), all = FALSE ),
           c( 'prey.core', 'prey.other' ),
           function( pair.samples ) {
               nisect <- nrow(pair.samples)
               ncore <- ncore.samples[ pair.samples$prey.core[[1]] ]
               data.frame(
                       p_value = phyper( nisect - 1, ncore, nsamples - ncore,
                                         nother.samples[ pair.samples$prey.other[[1]] ],
                                         lower.tail = FALSE, log.p = TRUE ),
                       stringsAsFactors = TRUE
               )
    } )
    return ( prey.distances )
}

#' Creates data frame of peptide modifications
#' @param peptideSeq peptide AA sequence
#' @param modificationsSeq colon-separated list of per AA modifications
#' @returnType data.frame
#' @return data.frame of non-zero modifications, including modification, offset and AA code
#' @author Alexey Stukalov
#' @export
ms_data.peptideModifications <- function( peptideSeq, modificationsSeq )
{
    modifications <- strsplit( modificationsSeq, ':', fixed=TRUE )[[1]]
    res <- data.frame()
    for ( i in 1:length(modifications) ) {
        if ( nchar( modifications[i] ) > 0 ) {
            res <- rbind( res, 
                data.frame( 
                    Peptide = peptideSeq,
                    Modification = modifications[i],
                    Offset = i - 2,
                    Residue = substr(peptideSeq, i-1, i-1 ),
                    stringsAsFactors = FALSE ) )
        }
    }
    return ( res )
}

#' Calculate the p-value of 0 hypothesis, that no phosphorylation change is observed.
#' 
#' @param modif_site_data data.frame for spectras registered at site
#' @param groups1 1st group of msruns
#' @param groups2 2nd group of msruns
#' @param modif modification to check
#' @returnType p-value of independence between 2 groups
#' @author Alexey Stukalov
#' @export
ms_data.modification.p_value <- function( modif_site_data, groups1, groups2, modif, 
                          msrun_col = 'msrun', modif_col = 'Modification',
                          quant_col = 'sc' )
{
    modif_counts <-  modif_site_data[ modif_site_data[,modif_col] == modif, ]
    if ( nrow( modif_counts ) == 0 ) return ( 1.0 ) # p_value=1 for degenerated case

    groups1_counts <- sapply( groups1, function( group ) sum( modif_counts[ modif_counts[, msrun_col] %in% group, quant_col ], na.rm = TRUE ) )
    groups2_counts <- sapply( groups2, function( group ) sum( modif_counts[ modif_counts[, msrun_col] %in% group, quant_col ], na.rm = TRUE ) )

    # p-value for hypothesis that there's no difference between groups
    m <- matrix( c( groups1_counts, groups2_counts ), nrow = 2, byrow = TRUE )
    p <- chisq.test( m, #correct = TRUE,
                     #p = rep( 1, 2 * length( groups1_counts ) ), 
                     rescale.p = TRUE )$p.value
    return ( p )
}

#' Strip down isoform info and merge down pairs of bait-prey interactions. 
#' @param ms_data 
#' @returnType data.frame
#' @return 
#' @author Alexey Stukalov
#' @export
ms_data.merge_noiso <- function( ms_data )
{
    ddply( subset( ms_data, msrun != '#glob#'), c("bait_ac_noiso", "prey_ac_noiso"),
        function( ms_data ) {
            data.frame(
                sample = paste( ms_data$sample ),
                msrun = paste( ms_data$msrun ),
                sample_origin = paste( ms_data$sample_origin ),
                msrun_origin = paste( ms_data$msrun_origin ),
                sc_max = max( ms_data$sc ),
                sc_sp_max = max( ms_data$sc_sp ),
                sc_sup_max = max( ms_data$sc_sup ),
                pc_max = max( ms_data$pc ),
                pc_sp_max = max( ms_data$pc_sp ),
                pc_sup_max = max( ms_data$pc_sup ),
                conditions = paste( sort(unique(ms_data$condition)), collapse = " "),
                stringsAsFactors = FALSE
                )
        }
    )
}

ms_data.average_msruns <- function(
    ms_data,
    cols_id = c( 'sample_origin', 'sample', 'bait_ac' ),
    col_prey = 'prey_ac',
    col_msrun = 'msrun',
    cols_msrun_related = c( 'msrun_origin' ),
    cols_to_avg = c( 'sc', 'sc_sp', 'sc_sup',
                     'pc', 'pc_sp', 'pc_sup',
                     'seqcov', 'seqcov_sp', 'seqcov_sup' )
){
    ms_noglob <- ms_data[ ms_data[ col_msrun ] != '#glob#',
                          c( cols_id, col_msrun, cols_msrun_related, col_prey, cols_to_avg ) ]
    ms_noglob.iactions <- unique( ms_noglob[,c(cols_id, col_prey)] )
    #print( nrow( ms_noglob ) )
    #print( ms_noglob.iactions[1:10,])
    #print( nrow( ms_noglob.iactions ) )
    msruns.df <- unique( ms_noglob[,c( cols_id, col_msrun )] )
    #print( nrow( msruns.df ) )
    #print( msruns.df[1:10,])
    ms_noglob.iactions.ext <- merge( msruns.df, ms_noglob.iactions, all.x = TRUE ) 
    #print( nrow( ms_noglob.iactions.ext ) )
    #print( ms_noglob.iactions.ext[1:10,])
    ms_noglob.weak <- merge( ms_noglob.iactions.ext, ms_noglob, all.x = TRUE )
    #print( nrow( ms_noglob.weak ) )
    #print( ms_noglob.weak[1:10,])
    # make a data frame with 0 sc for all replicates where interaction is unobserved (if it was observed in at least one)
    lapply ( cols_to_avg, function( col_name ) {
        ms_noglob.weak[ is.na( ms_noglob.weak[,col_name]), col_name ] <<- 0
    } )
    # average all data from tech.replicates
    return ( ddply( ms_noglob.weak, c(cols_id, col_prey), function( iact_data ) {
        res <- lapply( cols_to_avg, function( col_name ) mean( iact_data[,col_name]) )
        names( res ) <- cols_to_avg
        return ( as.data.frame( res, stringsAsFactors = FALSE ) )
    } ) )
}

#' Known direct interaction between baits or preys in matrix expansion
known_ppi.query <- function( idb.conn,
        bait_acs, prey_acs,
        internal.interaction.type = 'MI:0915', # physical association
        external.interaction.type = 'MI:0914', # association
        taxid = NULL,
        exclude.pubmedid = NULL,
        exclude.detection.methods = NULL,
        include.detection.methods = NULL
){
    # check interaction type up the vocabulary hierarchy
    check_vocab_filter <- function( tbl_prefix, vocabs ) {
        condition <- paste( "'", vocabs, "'", collapse = ',', sep = '' )
        condition <- ifelse( length( vocabs ) > 1,
                             paste( " IN (", condition, ")" ),
                             paste( " = ", condition ) )
        return ( paste( "(", tbl_prefix, ".id", condition,
                        " OR ", tbl_prefix, "_p.id ", condition,
                        " OR ", tbl_prefix, "_gp.id ", condition,
                        " OR ", tbl_prefix, "_ggp.id ", condition,
                        " OR ", tbl_prefix, "_gggp.id ", condition, ")", sep = '' ) )
    }
    # either bait-bait interaction (avoid duplicate pairs) or bait-prey(!in bait) interaction
    internal.filter <- ifelse( !is.null( internal.interaction.type ),
            paste("( ( (mol1.ac < mol2.ac AND mol2.id IN (", paste("'", intersect( bait_acs, prey_acs ), "'", sep="", collapse=","),"))
                       OR mol2.id IN (", paste("'", setdiff( prey_acs, bait_acs ), "'", sep="", collapse=","),"))
                  AND (", check_vocab_filter( "iact_type", internal.interaction.type ), "))", sep = '' ),
              "(0=1)" )
    external.filter <- ifelse( !is.null( external.interaction.type ), paste(
            "(/*(p1.experimentalrole = 'bait' OR p2.experimentalrole = 'bait') AND*/
                       ( NOT mol2.id IN (",paste("'", prey_acs, "'", sep="", collapse=","),") )
                  AND (", check_vocab_filter( "iact_type", external.interaction.type ), "))", sep = '' ),
              "(0=1)" )
    # exclude reports from specific publication
    publication.filter <- ifelse( !is.null( exclude.pubmedid ),
            paste("pub.pmid != '", exclude.pubmedid, "'", sep=''),
            '1=1' )
    method.filter.include <- ifelse( !is.null( include.detection.methods ),
            paste( "(", check_vocab_filter( "method", include.detection.methods ) ,")" ),
            "(1=1)" )
    method.filter.exclude <- ifelse( !is.null( exclude.detection.methods ),
            paste( "(NOT (", check_vocab_filter( "method", exclude.detection.methods ), "))" ),
            '(1=1)'
    )
    organism.filter <- ifelse( !is.null(taxid), paste('(mol2.taxid = \'', taxid, '\')', sep=''), '(1=1)' )

    dbGetQuery( idb.conn, paste("SELECT mol1.id as bait_ac, mol2.id as prey_ac,
            array_to_string(array_accum(DISTINCT iact.source_db), '|'::text) as source_dbs,
            array_to_string(array_accum(DISTINCT iact.interactiontype_id), '|'::text) as interactiontype_ids,
            array_to_string(array_accum(DISTINCT method.id), '|'::text) as experimentalmethod_ids,
            array_to_string(array_accum(DISTINCT coalesce(pub.first_author,'?')
                                          || '*' || coalesce(pub.journal, '?') || ' \\'' || coalesce(pub.title,'?') || '\\'' || ','
                                          || coalesce(pub.year,'?') || '(' || coalesce(pub.url, '?') || ')' ), '|'::text) as publications,
            array_to_string(array_accum(DISTINCT pub.pmid), '|'::text) as pubmed_ids,
            COUNT(DISTINCT e2p.publication_ac) as n_publications,
            COUNT(DISTINCT i2e.experiment_ac) as n_experiments
            FROM molecule mol1
            JOIN participant p1 ON p1.molecule_ac = mol1.ac
            JOIN interaction iact ON p1.interaction_ac = iact.ac
            JOIN participant p2 ON p2.interaction_ac = iact.ac 
            JOIN molecule mol2 ON p2.molecule_ac = mol2.ac
            JOIN interaction2experiment i2e ON i2e.interaction_ac = iact.ac
            JOIN experiment exp ON exp.ac = i2e.experiment_ac
            JOIN experiment2publication e2p ON e2p.experiment_ac = i2e.experiment_ac
            JOIN publication pub ON e2p.publication_ac = pub.ac
            JOIN controlled_vocab method ON exp.interactiondetectionmethod_id = method.id
            LEFT JOIN controlled_vocab method_p ON method_p.id = method.is_a
            LEFT JOIN controlled_vocab method_gp ON method_gp.id = method_p.is_a
            LEFT JOIN controlled_vocab method_ggp ON method_ggp.id = method_gp.is_a
            LEFT JOIN controlled_vocab method_gggp ON method_gggp.id = method_ggp.is_a
            JOIN controlled_vocab iact_type ON iact.interactiontype_id = iact_type.id
            LEFT JOIN controlled_vocab iact_type_p ON iact_type_p.id = iact_type.is_a
            LEFT JOIN controlled_vocab iact_type_gp ON iact_type_gp.id = iact_type_p.is_a
            LEFT JOIN controlled_vocab iact_type_ggp ON iact_type_ggp.id = iact_type_gp.is_a
            LEFT JOIN controlled_vocab iact_type_gggp ON iact_type_gggp.id = iact_type_ggp.is_a
            WHERE mol1.id in (",paste("'", bait_acs, "'", sep="", collapse=","),")
            AND (mol2.db = 'UniProt')
            AND (", internal.filter, " OR ", external.filter, ")
            AND (", method.filter.include, ")
            AND (", method.filter.exclude, ")
            AND (", publication.filter, ")
            AND (", organism.filter, ")
            /*AND ict.source_db in ('intact', 'mint', 'hprd')*/
            GROUP BY mol1.id, mol2.id"
    ) )
}

known_ppi.merge_noiso <- function( known_ppi )
{ 
    known_ppi$bait_ac_noiso <- uniprotac.removeiso( known_ppi$bait_ac )
    known_ppi$prey_ac_noiso <- uniprotac.removeiso( known_ppi$prey_ac )
    ddply( known_ppi, c("bait_ac_noiso", "prey_ac_noiso"),
    function( ki ) {
        data.frame(
            n_publications = sum( ki$n_publications ), # todo: sum is incorrect, might be duplicates
            n_experiments = sum( ki$n_experiments ), # todo: sum is incorrect, might be duplicates
            source_dbs = join_concat_lists( ki$source_dbs, '|' ),
            publications = join_concat_lists( ki$publications, '|' ),
            pubmed_ids = join_concat_lists( ki$pubmed_ids, '|' ),
            interactiontype_ids = join_concat_lists( ki$interactiontype_ids, '|' ),
            experimentalmethod_ids = join_concat_lists( ki$experimentalmethod_ids, '|' ),
            stringsAsFactors = FALSE
        )
    } )
}
