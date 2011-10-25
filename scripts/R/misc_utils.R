# Venn diagram miscellaneous utilities
# 
# Author: Alexey Stukalov (astukalov@cemm.oeaw.ac.at)
###############################################################################

require(Vennerable)

#' merge 2 dataset, keeping the order of the first one
merge.orderful <- function( ds1, ds2, ... )
{
    ds1_copy <- ds1
    ds1_copy[,'__merge_order__'] <- 1:nrow(ds1_copy)
    merge_res <- merge( ds1_copy, ds2, ... )
    # restore order and remove order column
    merge_res <- merge_res[ order(merge_res[,'__merge_order__']), setdiff( colnames( merge_res ), '__merge_order__' ) ]
    return ( merge_res )
}

#' Replace NA values with predefined value
replace.na <- function( a, na.replace = NA ) {
    return ( ifelse( is.na( a ), na.replace, a ) )
}

#' Extended lapply() -- passes the name of array element as a second
#' argument of user function.
#' @param l    list to apply to
#' @param func func( el, name_of_el, ... ) function to apply to each element 
#'             taking element and its name in the list as 2 mandatory arguments
#' @param ... extra argument for func
#' @returnType list
#' @return result of lapply
#' @author astukalov
#' @export
lapply.ex <- function( l, func, ... )
{
    res <- lapply( seq_along( l ), function( ix ) { func( l[[ix]], names(l)[[ix]], ... ) } )
    names( res ) <- names( l )
    return ( res )
}

shift <- function( v, n )
{
    n <- n %% length( v )
    if ( n == 0 ) return ( v )
    return ( v[ c( (n+1):length(v), 1:n) ] )
}

matrix.random.rowshift <- function( mtx )
{
    shifts <- sample.int( ncol(mtx), size = nrow(mtx), replace = TRUE ) - 1
    matrix( unlist( lapply( 1:nrow(mtx), function(rowix) shift( mtx[rowix,], shifts[rowix] ) ) ),
            nrow = nrow( mtx ), byrow = TRUE )
}

#' non-throwing wrapper around Venn()
plotVenn2Safe <- function( vennSets, ... )
{
    vennObj <- Venn( vennSets )
    wrn <- NA
    isets <- sapply(vennObj@IntersectionSets, length)
    if ( sum(isets) == 0 )  {
        wrn <- paste( paste(names(vennSets), collapse=' & ' ), 'empty' )
    } else if ( sum(isets[c('10','11')]) == 0 ) {
        wrn <- paste( names(vennSets)[[1]], 'empty' )
    } else if ( sum(isets[c('01','11')]) == 0 ) {
        wrn <- paste( names(vennSets)[[2]], 'empty' )
    } else if ( sum(isets[c('01','10')]) == 0 ) {
        wrn <- paste( paste(names(vennSets), collapse=' & ' ), 'equal' )
    }
    if (is.na(wrn)) {
        try( plot( vennObj, ... ) )
    } 
    else {
        grid.text( wrn, gp = gpar(cex=2) )
    }
}

#' prepare sets out of dataframe to be used with Vennerable::Venn
prepareVennSets <- function( dframe, keys, indicators )
{
    res <- lapply( indicators, function( indicator ) {
            print(indicator)
            unique( apply( dframe[ dframe[,indicator], keys ], 1, function( x ) paste( x, collapse=' ' ) ) )
        } )
    names( res ) <- indicators
    return ( res )
}

#' prepare sets out of dataframe to be used with Vennerable::Venn
prepareVennSets <- function( dframes, keys )
{
    dframes <- dframes[ sapply( dframes, nrow ) > 0 ] # skip empty
    if ( length( keys ) > 1 ) {
        res <- lapply( dframes, function( dframe ) {
                apply( unique( dframe[,keys] ), 1, function( row ) {
                        paste( row[keys], collapse=' ' ) } )
            } )
    } else {
        res <- lapply( dframes, function( dframe ) unique( dframe[,keys] ) )
    }
    #print( res )
    return ( res )
}

#' Identifies in what datasets each unique element is present and returns membership data.frame 
makeMembershipTable <- function( datasets, columns )
{
    allkeys <- NULL
    for ( dataset in datasets ) {
        thiskeys <- as.data.frame( unique( dataset[,columns] ) )
        colnames(thiskeys) <- columns
        if ( is.null(allkeys) ) {
            allkeys <- thiskeys
        }
        else {
            allkeys <- merge( allkeys, thiskeys, suffixes = c('.a','.b'), by = columns, all = TRUE )
        }
    }
    membership <- allkeys
    if ( is.null( names( datasets ) ) ) {
        datasets.names <- paste('dataset',1:length(datasets),sep='')
    }
    else {
        datasets.names <- names(datasets)
    }
    for ( i in 1:length(datasets) ) {
        dataset <- datasets[[i]]
        colname <- datasets.names[[i]]
        thiskeys <- as.data.frame( unique( dataset[,columns] ) )
        colnames(thiskeys) <- columns
        thiskeys[,colname] <- 1:nrow(thiskeys)
        thismembership <- merge( allkeys, thiskeys, by = columns, all.x = TRUE, all.y = TRUE )
        membership[,colname] <- !is.na( thismembership[,colname] )
    }
    return ( membership )
}

#' counts elements of two different sets, per group
#' @param df output of makeMembershipTable() 
#' @param group_columns column name  defining the grouping of elements
#' @param mem_column1 column name for flag of 1st set membership
#' @param mem_column2 column name for flag of 2nd set membership
#' @param order_sets if true, sets are ordered by number of elements, big first
#' @returnType data.frame
#' @return 
#' @author astukalov
#' @export
membership_counts <- function( df, group_columns, mem_column1, mem_column2, order_sets = TRUE )
{
    ddply( df, group_columns, function( sdf ) {
            A = sum( sdf[,mem_column1] & !sdf[,mem_column2 ] )
            AB = sum( sdf[,mem_column1] & sdf[,mem_column2 ] )
            B = sum( !sdf[,mem_column1] & sdf[,mem_column2 ] )
            if ( order_sets & A < B ) {
                t = A
                A = B
                B = t
            }
            S = A + AB + B
            if ( S == 0 ) return ()
            data.frame( AB = AB, A = A, B = B, S = S, rA = A/S, rAB = AB/S, rB = B/S )
        } )
}

#' calculates statistic for data generated by membership_counts()
#' 
#' @param cnts output of membership_counts() 
#' @returnType data.frame
#' @return mean ratio of A-only elements, mean ratio of A-B overlap, etc
#' @author astukalov
#' @export
membership_stats <- function( cnts )
{
    data.frame( A = mean( cnts$rA ), sd.A = sd( cnts$rA ),
        AB = mean( cnts$rAB ), sd.AB = sd( cnts$rAB ),
        B = mean( cnts$rB ), sd.B = sd( cnts$rB ) )
}

#' Calculate distance from bait to prey-prey interaction as the number direct PP interactions
#' linking prey and the bait.
#' Since among PP interactions there could be physical associations, which might include non-direct interactions,
#' direct interactions are prioritized when calculating this distance.
#' @param bait_prey1_prey2 data.frame of direct interactions
#' @param res_col name of column to add to dataset
#' @returnType integer vector
#' @return vector, each component corresponds to the row of bait_prey1_prey2
#' @author astukalov
distance_to_bait <- function( bpp ) {
    # initialize distance -- bait-bounding interaction gets 0 distance
    res <- cbind( bpp, data.frame(
        distance = ifelse( bpp$bait_ac_noiso == bpp$prey_ac_noiso.1 |
                           bpp$bait_ac_noiso == bpp$prey_ac_noiso.2, 0, NA ),
        nindirect  = 0,
        stringsAsFactors = FALSE
    ) )
    res[ !is.na( res$distance ) & !res$is_direct, 'nindirect' ] <- 1
    n_linked.last <- 0 # of proteins already linked to bait
    while ( n_linked.last < sum( !is.na( res$distance ) ) ) {
        n_linked.last <- sum( !is.na( res$distance ) )
        print( paste( "Max distance", max(res$distance, na.rm = TRUE ), "linked", n_linked.last ) )
        res.new <- res
        for ( ix in 1:nrow(res) ) {
            linked.iactions <- 1:nrow(res) != ix & !is.na( res$distance ) &
                res$bait_ac_noiso ==  res[ ix, 'bait_ac_noiso'] &
                ( res$prey_ac_noiso.1 == res[ ix, 'prey_ac_noiso.1']
                | res$prey_ac_noiso.1 == res[ ix, 'prey_ac_noiso.2']
                | res$prey_ac_noiso.2 == res[ ix, 'prey_ac_noiso.1']
                | res$prey_ac_noiso.2 == res[ ix, 'prey_ac_noiso.2'] )
            if ( any( linked.iactions ) ) {
                # prefer direct interactions over associations
                min.indirect <- min( res[ linked.iactions, 'nindirect' ] )
                linked.iactions <- linked.iactions & res$nindirect == min.indirect
                new.distance <- min( res[ linked.iactions, 'distance' ] ) + 1
                new.nindirect <- min.indirect + all( !res[ linked.iactions, 'is_direct' ] )
                # update distance if not defined before or less indirect interactions
                if ( is.na( res[ ix, 'distance' ] ) | res[ ix, 'nindirect' ] > new.nindirect ) {
                    res.new[ ix, 'distance' ] <- new.distance
                    res.new[ ix, 'nindirect' ] <- new.nindirect
                }
            }
        }
        res <- res.new
    }
    return ( res )
}

top.factor.values <- function( ds, column_name, n = 10 )
{
    freq <- table( ds[,column_name ] )
    freq <- freq[ order( freq, decreasing = TRUE )[1:n] ]
    res <- ds[ ds[,column_name] %in% names(freq), ]
    for ( colname in colnames(res) ) {
        if ( is.factor(res[,colname]) ) {
            res[,colname] <- factor(as.character(res[,colname]))
        }
    }
    return ( res )
}

# joins lists collapsed to strings
# avoiding duplicates
join_concat_lists <- function( concat_lists, sep = ' ', flanks = NULL )
{
    if ( is.null(concat_lists) || is.na(concat_lists) || length(concat_lists) == 0 ) return ( NULL )
    if ( !is.null( flanks ) ) {
        # strip flanks
        concat_lists <- lapply( concat_lists, function( s ) {
            if ( substr( s, 1, 1 ) == flanks[[1]] ) s <- substr( s, 2, nchar( s ) )
            if ( substr( s, nchar( s ), nchar( s ) ) == flanks[[2]] ) s <- substr( s, 1, nchar( s ) - 1 )
            s 
        } )
    }
    items_union  <- sort( unique( unlist( sapply( concat_lists,
                          function( concat_list ) strsplit( concat_list, sep, fixed = TRUE ) ) ) ) )
    res <- paste( items_union, collapse = sep )
    if ( !is.null( flanks ) ) {
        # add flanks
        res <- paste( flanks[[1]], res, flanks[[2]], sep = '' )
    }
    return ( res )
}

#' Order elements in the partition:
#' first -- by order of the clusters,
#' second -- alphabetiaclly within the cluster.
#' @param partition 
#' @param clusters.ordered 
#' @param cluster_col 
#' @param element_col 
#' @param decreasing 
#' @returnType 
#' @return ordered list of elements
#' @author astukalov
#' @export
partition.elements.order <- function( ordered.clusters, partition, 
        cluster_col = 'cluster', element_col = 'element',
        element_sort_cols = element_col,
        decreasing = FALSE )
{
    return ( unlist( sapply( ordered.clusters,
        function( clu ) {
            clu_mask <- partition[ , cluster_col ] == clu
            clu_elms <- partition[ clu_mask, ]
            return( clu_elms[ do.call( order, c( lapply( element_sort_cols, function( colname ) clu_elms[ , colname ] ),
                                     decreasing = decreasing ) ),
                              element_col ] )
    } ) ) )
}

elements_order.clusters_induced <- function( clustering,
        cluster_col = 'cluster', element_col = 'element',
        element_sort_cols = element_col,
        clusters.dgram = NULL, decreasing = decreasing
){
    cluster.sizes <- table( clustering[,cluster_col] )
    cluster.sizes <- cluster.sizes[ sort( names(cluster.sizes) ) ]
    res <- list()
    if ( !is.null( clusters.dgram ) ) {
        cluster.sizes <- cluster.sizes[ clusters.dgram$order ]
        res$clusters.dgram <- clusters.dgram
        res$elements.dgram <- dgram.set.n_members( as.dendrogram( clusters.dgram ), cluster.sizes )
    }
    else {
        res$clusters.dgram <- NULL
        res$elements.dgram <- NULL
    }
    # order elements by label within clusters
    res$clusters.order <- names( cluster.sizes )
    res$elements.order <- partition.elements.order( res$clusters.order,
            clustering, cluster_col = cluster_col, element_col = element_col,
            element_sort_cols = element_sort_cols,
            decreasing = decreasing )
    return ( res )
}
