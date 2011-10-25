# Plotting clusterings and MSMS trajectories.
# 
# Author: Alexey Stukalov
###############################################################################

source( file.path( ppi_scripts_path, "BIMAP.R" ) )
source( file.path( ppi_scripts_path, "misc_utils.R" ) )

require( ggplot2 )
require( lattice )
require( latticeExtra )

plot.signal_noise_sigma <- function( bimap.walk, densities = FALSE, clip.factor = 4 )
{
    walk.clusterings.frame <- bimap.walk@clusterings.walk
    walk.priors.frame <- bimap.walk@priors.walk

    if ( densities ) {
        xvals <- c( walk.clusterings.frame$noise.peak, 
                    bimap.walk@signals$signal )
        xrange <- c( min(xvals), max(xvals) ) 
        signal.grid <- seq( xrange[[1]], xrange[[2]], by = (xrange[[2]] - xrange[[1]])/100 )
        noise.density <- density( walk.clusterings.frame$noise.peak )
        noise.max <- max( noise.density$y )
        baseline.max <- max( apply( walk.priors.frame, 1, function( norm ) {
                dnorm( 0, sd = norm[['baseline.signal.sigma']] )
            }) )
        signal.density <- density( bimap.walk@signals$signal )
        signal.max <- max( signal.density$y )
        ymax <- max( min( clip.factor * noise.max, baseline.max ), min( clip.factor * baseline.max, noise.max ), signal.max )
        plot( noise.density, xlim = xrange, ylim = c( 0, ymax ), col = "blue", lwd = 2, xlab = 'signal' )
        apply( walk.priors.frame, 1, function( norm ) {
                lines( signal.grid, dnorm( signal.grid, mean = norm[['baseline.signal']], 
                                                        sd = norm[['baseline.signal.sigma']] ),
                       col = rgb(1,0.5,0,0.02), lwd = 4 )
            })
        lines( signal.density, xlim = xrange, col = "red", lwd = 2 )
    }
    else {
        yvals <- c( walk.clusterings.frame$baseline.peak, 
            walk.clusterings.frame$noise.peak,
            walk.priors.frame$baseline.signal + walk.priors.frame$baseline.signal.sigma,
            walk.priors.frame$baseline.signal - walk.priors.frame$baseline.signal.sigma )
        yranges <- c( min( yvals ), max( yvals ) )
        
        plot( walk.clusterings.frame$time, 
            walk.clusterings.frame$baseline.peak,
            col="red", type='l', lwd=2, 
            ylim = yranges
        )
        polygon( unlist( c( walk.priors.frame$time, rev( walk.priors.frame$time ) ) ),
            unlist( c( walk.priors.frame$baseline.signal + 2 * walk.priors.frame$baseline.signal.sigma,
                    rev( walk.priors.frame$baseline.signal - 2 * walk.priors.frame$baseline.signal.sigma ) ) ),
            col = rgb( 1.0, 0.2, 0.2, 0.2 ), border = NA )
        lines( walk.clusterings.frame$time, 
            walk.clusterings.frame$noise.peak, 
            col="blue", type='l', lwd = 1 )
        lines( walk.priors.frame$time, 
            walk.priors.frame$baseline.signal + walk.priors.frame$baseline.signal.sigma, 
            col="red", type='l', lty = 3 )
        lines( walk.priors.frame$time, 
            walk.priors.frame$baseline.signal - walk.priors.frame$baseline.signal.sigma, 
            col="red", type='l', lty = 3 )
        rug( side = 2, walk.priors.frame$baseline.signal - 2 * walk.priors.frame$baseline.signal.sigma,
            col = "orange" )
        rug( side = 2, walk.priors.frame$baseline.signal + 2 * walk.priors.frame$baseline.signal.sigma,
            col = "orange" )
        rug( side = 2, walk.clusterings.frame$baseline.peak, col = "red" )
        rug( side = 2, walk.clusterings.frame$noise.peak, col = "blue" )
    }
}

# group ascending integers into ranges
ranges <- function( int_list ) {
    if ( length(int_list) == 0 ) return ( list() )
    res <- list()
    start <- int_list[[1]]
    last <- int_list[[1]]
    if ( length(int_list) > 1 ) {
        for ( i in 2:length(int_list) ) {
            if ( int_list[[ i ]] > last + 1 ) {
                res <- append( res, list(c( start, last )) )
                start <- int_list[[ i ]]
            }
            last <- int_list[[ i ]]
        }
    }
    res <- append(res, list(c(start, last)) )
    return ( res )
}

plot.crosscluster <- function( proteins, samples, 
                               proteins.cluster, samples.cluster,
                               signals.matrix, signals.matrix.sd,
                               signals.frame = NULL,
                               show.abundance = TRUE,
                               bait.col = 'red',
                               border.col = NULL,
                               border.above = TRUE,
                               border.below = TRUE )
{
    #usr <- par("usr") # remember current clip
    bait_acs  <- intersect( samples$bait_ac, proteins$protein_ac )
    #clip( f[[1]], samples.cluster.pos[[2]], 
    #      proteins.cluster.pos[[1]], proteins.cluster.pos[[2]] )
    #print( which( rownames(signals.matrix) == proteins.cluster ) )
    #print( which( colnames(signals.matrix) %in% samples ) )
    rowIx <- which( rev( rownames(signals.matrix) ) == proteins.cluster )
    pushViewport( viewport( layout = grid.layout( nrow = nrow( proteins ), ncol = nrow( samples ) ) ) )
    for ( cur_bait_ac in bait_acs ) {
        bait_samples <- subset( samples, bait_ac == cur_bait_ac )
        rowIx <- which( rev( proteins$protein_ac ) == cur_bait_ac )
        for ( curSample in bait_samples$sample ) {
            colIx <- which( samples$sample == curSample )
            pushViewport( viewport( layout.pos.row = rowIx,
                                    layout.pos.col = colIx ) )
            panel.rect( 0.01, 0.01, 0.99, 0.99, 
                  border = bait.col, lwd = 3 )
            popViewport()
        }
    }
    popViewport()

    main.text = trellis.par.get('par.main.text')
    main.text$cex = 0.2
    label.color <- 'black'
    label.font <- 1
    mean.signal = signals.matrix[ as.character( proteins.cluster ), as.character( samples.cluster ) ]
    if ( is.na(mean.signal) ) {
        warning( "Cluster ", proteins.cluster, "x", samples.cluster, 
                 "have no mean signal" )
    }
    sd.signal = signals.matrix.sd[ as.character( proteins.cluster ), as.character( samples.cluster ) ]
    if ( is.na( sd.signal ) ) {
        warning( "Cluster ", proteins.cluster, "x", samples.cluster, 
                 "have no signal sd" )
    } 
    main.title = paste( format(mean.signal,digits=3), "\u0b1",
                        format(sd.signal,digits=3), sep="" )
    if ( !is.null( signals.frame ) ) {
        # embed sampling plot
        signal.frame <- subset( signals.frame, samples.cluster == samples.cluster
                  & proteins.cluster == proteins.cluster )
        if ( nrow( signal.frame ) > 0 ) {
            signal.frame$ix <- 1:nrow(signal.frame)
            xlim <- c( 1, nrow(signal.frame ) )
            ylim <- c( min( signal.frame$signal ), max( signal.frame$signal ) )
            lpoints( (signal.frame$ix - xlim[[1]]) / (xlim[[2]]-xlim[[1]]),
                     (signal.frame$signal - ylim[[1]]) / (ylim[[2]]-ylim[[1]]),
                      col = rgb( 0,0,0,0.2), pch = 3, cex = 0.1 )
#                pan( xyplot( signal~ix, signal.frame[1:5,], pch = 3, cex = 0.2, aspect = 'fill',
#                        main = NULL#main.title
#                        , xlab = NULL, ylab = NULL, scales = list(alternating=0),
#                        col = ifelse( is.bait.cluster, 'black', 'white' ),
#                        par.settings = list( main.text = main.text )
#                    ), 
#                    newpage = FALSE, more = TRUE, position = c(0.0,0.0,1.0,1.0) ) 
        }
        else {
            ltext( 0.5, 0.5, paste( main.title, "no samples" ), col = label.color, font = label.font )
        }
    }
    if ( show.abundance ) {
        panel.text( 0.5, 0.5, main.title, col = label.color, font = label.font )
    }
    #do.call("clip", as.list(usr)) # restore clipping
    if ( !is.null( border.col ) ) {
        if ( border.above && border.below ) {
            panel.rect( 0.0, 0.0, 1.0, 1.0, 
                  border = border.col, lwd = 3 )
        } else {
            panel.segments( c( 0.0, 1.0 ), c( 0.0, 0.0 ),
                            c( 0.0, 1.0 ), c( 1.0, 1.0 ),
                  col = border.col, lwd = 3 )
            if ( border.above ) {
                panel.segments( 0.0, 0.0, 1.0, 0.0,
                      col = border.col, lwd = 3 )
            }
            if ( border.below ) {
                panel.segments( 0.0, 1.0, 1.0, 1.0,
                      col = border.col, lwd = 3 )
            }
        }
    }
}

plot.crosscluster.measurements <- function(
    proteins, sample, samples,
    msrun.info,
    signals.matrix, measurements )
{
    pushViewport( viewport( layout.pos.row = 1, 
                            layout.pos.col = which( samples == sample ) ) )

    pushViewport( viewport( x = 1.0, y = 0.0,
                            width = 0.8, height = 0.3, 
                            just = c( 'right', 'bottom' ) ) )
    msruns <- msrun.info$msrun[ msrun.info$sample == sample ]
    submeas <- subset( measurements, msrun %in% msruns & prey_ac %in% proteins,
                       select = c('prey_ac', 'msrun', 'sc') )
    panel.grid( h = length(proteins)-1, v = length(msruns)-1 )
    if ( nrow( submeas ) > 0 ) {
        colpos <- sapply( submeas$msrun, function( msrun ) which( msruns == msrun )/length(msruns) ) - 0.5/length(msruns)
        rowpos <- sapply( submeas$prey_ac, function( prey_ac ) which( proteins == prey_ac )/length(proteins) ) - 0.5/length(proteins)
        ltext( colpos, rowpos, labels = submeas$sc, col = 'darkgray', cex = unit( 0.7, 'npc' ) )
    }
    popViewport()
    popViewport()
}

dgram.set.n_members <- function( dgram, sizes )
{
    calc.sizes <- function( node ) {
        if ( is.leaf(node) ) {
            return ( sizes[[ attr( node, 'label' ) ]] )
        }
        else {
            return ( do.call( sum, lapply( 1:length(node), 
                        function( nodeIx ) { calc.sizes( node[[ nodeIx ]] ) } ) ) )
        }
    }

    return ( dendrapply( dgram, function( node ) {
            attr( node, 'members' ) <- calc.sizes( node )
            return ( node )
        } ) )
}

addPositions <- function(x, order, leaf.pos)
{
    if (!is.null(attr(x, "position"))) return(x)
    else if (is.leaf(x))
    {
        attr(x, "position") <-
            list(x = leaf.pos[[which(x == order)[1]]],
                y = attr(x, "height"))
        #print(attr(x,'position'))
        return(x)
    }
    else
    {
        for (i in seq_along(x))
        {
            x[[i]] <- addPositions(x[[i]], order, leaf.pos)
        }
        attr(x, "position") <-
            list(x = mean(sapply(x, function(x) attr(x, "position")$x )),
                y = attr(x, "height"))
        return(x)
    }
}

edgeLocation <-
    function(pos.node, pos.child, type, ...)
{
    switch(type,
        rectangle = {
            data.frame(x0 = c(pos.node$x, pos.child$x),
                y0 = c(pos.node$y, pos.node$y),
                x1 = c(pos.child$x, pos.child$x),
                y1 = c(pos.node$y, pos.child$y),
                ..., stringsAsFactors = FALSE) ## 'col' can be strings
        },
        triangle = {
            data.frame(x0 = pos.node$x,  y0 = pos.node$y,
                x1 = pos.child$x, y1 = pos.child$y,
                ..., stringsAsFactors = FALSE) ## 'col' can be strings
        })
}

dendrogramGrob.fixed <- function (x, ord = order.dendrogram(x),
    side = c("right", "top"),
    add = list(), size = 5, size.add = 1,
    type = c("rectangle", "triangle"), ...) 
{
    if (size <= 0) 
        return(textGrob(label = NULL))
    type <- match.arg(type)
    native.height <- attr(x, "height")
    leaf.members <- vector(length=length(ord))
    dendrapply(x, function(x) {
            if ( is.leaf(x) ) leaf.members[which(x == ord)] <<- attr(x, 'members')
        })
    leaf.pos <- cumsum(leaf.members)
    native.xscale <- c(1, sum(leaf.members)) + c(-1, 1) * lattice.getOption("axis.padding")$factor
    #print(leaf.pos)
    xpos <- addPositions(x, ord, leaf.pos)
    nnodes <- 0
    dendrapply(xpos, function(x) {
            if (!is.leaf(x)) 
                nnodes <<- nnodes + 1
        })
    xseg <- vector(mode = "list", length = nnodes)
    i <- 0
    getSegments <- function(x, ...) {
        if (!is.leaf(x)) {
            i <<- i + 1
            pos.node <- attr(x, "position")
            xseg[[i]] <<- do.call(rbind, lapply(x, function(child) {
                        pos.child <- attr(child, "position")
                        edgeLocation(pos.node, pos.child, type = type, 
                            ...)
                    }))
        }
    }
    dendrapply(xpos, getSegments)
    all.segs <- do.call(rbind, xseg)
    nadd <- length(add)
    nleaf <- length(ord)
    native.unit <- 1/diff(native.xscale)
    #print(leaf.members)
    switch(side,
        right = {
            key.layout <-
                grid.layout(nrow = 1, ncol = 1 + nadd,
                    heights = unit(1, "null"),
                    widths =
                        unit(c(rep(size.add, length = nadd), size),
                            c(rep("lines", nadd), "lines")),
                    respect = FALSE)
            key.gf <- frameGrob(layout = key.layout)
            ## key.gf <- placeGrob(key.gf, rectGrob(gp = gpar(fill = "pink")))
            for (i in seq_len(nadd))
            {
                addi <- add[[i]]
                typei <- names(add)[i]
                switch(typei,
                    rect = {
                        key.gf <-
                            placeGrob(key.gf,
                                rectGrob(y = (leaf.pos - native.xscale[1]) * native.unit,
                                    height = native.unit,
                                    gp = do.call(gpar, addi)),
                                row = 1, col = i)
                    })
            }
            key.gf <-
                placeGrob(key.gf, 
                    with(all.segs,
                        segmentsGrob((y0 / native.height),
                            (x0 - native.xscale[1]) * native.unit,
                            (y1 / native.height),
                            (x1 - native.xscale[1]) * native.unit)),
                    row = 1, col = 1 + nadd)
            key.gf
        },
        top = {
            key.layout <-
                grid.layout(nrow = 1 + nadd, ncol = 1,
                    widths = unit(1, "null"),
                    heights =
                        unit(c(size, rep(size.add, length = nadd)),
                            c("lines", rep("lines", nadd))),
                    respect = FALSE)
            
            key.gf <- frameGrob(layout = key.layout)
            ## key.gf <- placeGrob(key.gf, rectGrob(gp = gpar(fill = "pink")))
            
            for (i in seq_len(nadd))
            {
                addi <- add[[i]]
                typei <- names(add)[i]
                switch(typei,
                    rect = {
                        key.gf <-
                            placeGrob(key.gf,
                                rectGrob(x = (leaf.pos - native.xscale[1]) * native.unit,
                                    width = native.unit,
                                    gp = do.call(gpar, addi)),
                                row = 1 + i, col = 1)
                    })
            }
            key.gf <-
                placeGrob(key.gf, 
                    with(all.segs,
                        segmentsGrob((x0 - native.xscale[1]) * native.unit,
                            (y0 / native.height),
                            (x1 - native.xscale[1]) * native.unit,
                            (y1 / native.height))),
                    row = 1, col = 1)
            key.gf
        })
}

BIMAP.signals_matrix.bihclust <- function( bimap.props, signal.na.subst = -100 )
{
    proteins.clusters <- sort( unique( bimap.props$proteins.clusters$proteins.cluster ) )
    samples.clusters <- sort( unique( bimap.props$samples.clusters$samples.cluster ) )
    signals.matrix <- matrix( NA, nrow = length( proteins.clusters ),
                                  ncol = length( samples.clusters ) )
    rownames( signals.matrix ) <- proteins.clusters
    colnames( signals.matrix ) <- samples.clusters
    signals.matrix.sd <- signals.matrix
    existing.proteins.clusters <- rownames( bimap.props$signals.mean )
    existing.samples.clusters <- colnames( bimap.props$signals.mean )
    signals.matrix[ existing.proteins.clusters, existing.samples.clusters ] <- 
        bimap.props$signals.mean[ existing.proteins.clusters, existing.samples.clusters ]
    if ( !is.null(bimap.props$signals.sd) ) {
        signals.matrix.sd[ existing.proteins.clusters, existing.samples.clusters ] <- 
            bimap.props$signals.sd[ existing.proteins.clusters, existing.samples.clusters ]
    }

    signals.clu.matrix <- signals.matrix
    signals.clu.matrix[ is.na( signals.clu.matrix ) ] <- signal.na.subst

    # hierarchical clustering of proteins
    if ( nrow( signals.clu.matrix ) > 1 ) {
        protein.clusters.dgram <- hclust( dist( signals.clu.matrix ),
                                          members = rownames( signals.clu.matrix )  )
    } else {
        protein.clusters.dgram <- NULL
    }
    proteins.ordering <- elements_order.clusters_induced( bimap.props$proteins.clusters,
            cluster_col = 'proteins.cluster', element_col = 'protein_ac',
            element_sort_cols = 'protein_label',
            clusters.dgram = protein.clusters.dgram,
            decreasing = TRUE )

    # hierarchical clustering of states
    if ( ncol( signals.clu.matrix ) > 1 ) {
        sample.clusters.dgram <- hclust( dist( t( signals.clu.matrix ) ),
                                          members = colnames( signals.clu.matrix )  )
    } else {
        sample.clusters.dgram <- NULL
    }
    samples.ordering <- elements_order.clusters_induced( bimap.props$samples.clusters,
            cluster_col = 'samples.cluster',
            element_col = ifelse( 'msrun' %in% colnames( bimap.props$samples.clusters ), 'msrun', 'sample' ),
            element_sort_cols = ifelse( 'msrun' %in% colnames( bimap.props$samples.clusters ), c( 'sample', 'msrun' ), 'sample' ),
            clusters.dgram = sample.clusters.dgram,
            decreasing = FALSE )
    return ( list( signals.matrix = signals.matrix, signals.matrix.sd = signals.matrix.sd,
                   proteins.ordering = proteins.ordering, samples.ordering = samples.ordering ) )
}

plot.bimap <- function( bimap.props, protein.info,
                                  sample.info, msrun.info,
                                  protein.info.ex = NULL, measurements = NULL,
                                  plot.samples = FALSE,
                                  show.abundance = TRUE, 
                                  show.protein.id = TRUE,
                                  show.sample.id = TRUE, show.msruns = FALSE,
                                  show.borders = FALSE,
                                  extended.result = FALSE,
                                  protein_name_col = 'short_id',
                                  protein_description_col = NULL,
                                  title = NULL,
                                  cell.func = function( sc, pc, seqcov, seqlen ) {
                                      return ( log( sc / seqlen ) ) },
                                  col = heat.colors,
                                  cells.col = heat.colors,
                                  cells.off.alpha = 0.5,
                                  aspect = 0.5,
                                  bait.col = 'green',
                                  grid.col = 'grey', grid.lty = 2, grid.lwd = 2,
                                  protein.label.width = 1.5, sample.label.width = protein.label.width )
{
    proteins.clusters <- ddply( bimap.props$proteins.clusters, c( 'proteins.cluster' ), function( cluster ) {
            data.frame( size = nrow( cluster ),
                        stringsAsFactors = FALSE ) }
    )
    colnames( proteins.clusters ) <- c( 'serial', 'size' )
    proteins.clusters$label <- proteins.clusters$serial
    rownames( proteins.clusters ) <- proteins.clusters$serial
    proteins <- data.frame(
        protein_ac = bimap.props$proteins.clusters$protein_ac,
        short_label = bimap.props$proteins.clusters$protein_ac,
        stringsAsFactors = FALSE
    )
    rownames( proteins ) <- proteins$protein_ac
    if ( !is.null(protein.info.ex) & !is.null( protein_name_col ) ) {
        proteins$short_label <- protein.info.ex[ proteins$protein_ac, protein_name_col ]
    }
    if ( !is.null(protein.info.ex) & !is.null( protein_description_col ) ) {
        proteins$description <- protein.info.ex[ proteins$protein_ac, protein_description_col ]
    }
    bimap.props$proteins.clusters <- within( bimap.props$proteins.clusters, {
        protein_label <- proteins[ bimap.props$proteins.clusters$protein_ac, 'short_label' ]
    } )

    if ( show.msruns ) {
        # join sample clustering with msruns information 
        bimap.props$samples.clusters <- merge( bimap.props$samples.clusters, msrun.info, by = c( 'sample' ),
                                            all.x = TRUE, all.y = FALSE )
    }
    samples.clusters <- ddply( bimap.props$samples.clusters, c( 'samples.cluster' ), function( cluster ) {
            data.frame( size = nrow( cluster ),
                stringsAsFactors = FALSE ) }
    )
    colnames( samples.clusters ) <- c( 'serial', 'size' )
    samples.clusters$label <- samples.clusters$serial 
    rownames( samples.clusters ) <- samples.clusters$serial
    samples <- data.frame(
        samples.cluster = bimap.props$samples.clusters$samples.cluster,
        bait_ac = sample.info[ bimap.props$samples.clusters$sample, 'bait_ac' ],
        stringsAsFactors = FALSE
    )
    if ( show.msruns ) {
        samples$sample <- bimap.props$samples.clusters$msrun
    } else {
        samples$sample <- bimap.props$samples.clusters$sample
    }
    samples$short_label = samples$sample
    rownames( samples ) <- samples$sample
    #print( samples.clusters )
    
    # compose proteins labels
    proteins$axis_label <- proteins$short_label
    if ( show.protein.id ) {
        proteins$axis_label <- paste( proteins$axis_label,
                "(", proteins$protein_ac, ")", sep='' )
    }
    proteins$samples_axis_label <- proteins$axis_label
    if ( 'description' %in% colnames(proteins) ) {
        proteins$axis_label <- paste( proteins$axis_label,
                proteins$description, sep=' | ' )
    }

    # compose sample labels
    samples$axis_label <- proteins[ samples$bait_ac, 'samples_axis_label' ]
    if ( show.sample.id ) {
        samples$axis_label <- paste( samples$axis_label, ':', samples$sample )
    }

    signals.matrix.bihclust <- BIMAP.signals_matrix.bihclust( bimap.props )
    signals.matrix <- signals.matrix.bihclust$signals.matrix
    signals.matrix.sd <- signals.matrix.bihclust$signals.matrix.sd

    # prepare binary matrix of enabled clusters
    bimap.matrix <- matrix( 0, nrow = nrow( proteins.clusters ), 
                            ncol = nrow( samples.clusters ) )
    rownames( bimap.matrix ) <- rownames( proteins.clusters )
    colnames( bimap.matrix ) <- rownames( samples.clusters )
    for ( ccIx in 1:nrow(bimap.props$cross.clusters) ) {
        cc <- bimap.props$cross.clusters[ccIx,]
        bimap.matrix[ as.character( cc[['proteins.cluster']] ),
                   as.character( cc[['samples.cluster']] ) ] <- 1
    }

    proteins <- proteins[ signals.matrix.bihclust$proteins.ordering$elements.order, ]
    proteins.clusters <- proteins.clusters[ signals.matrix.bihclust$proteins.ordering$clusters.order, ]

    samples <- samples[ signals.matrix.bihclust$samples.ordering$elements.order, ]
    samples.clusters <- samples.clusters[ signals.matrix.bihclust$samples.ordering$clusters.order, ]

    # reorder signals matrix to match dendrograms
    signals.matrix <- signals.matrix[ proteins.clusters$serial, samples.clusters$serial ]
    signals.matrix.sd <- signals.matrix.sd[ proteins.clusters$serial, samples.clusters$serial ]
    bimap.matrix <- bimap.matrix[ proteins.clusters$serial, samples.clusters$serial ]

    # plot dendrograms
    legend = list()
    if ( !is.null( signals.matrix.bihclust$samples.ordering$elements.dgram ) ) {
        legend$top <- list( fun = dendrogramGrob.fixed, args = list( signals.matrix.bihclust$samples.ordering$elements.dgram, side='top' ) )
    }
    if ( !is.null( signals.matrix.bihclust$proteins.ordering$elements.dgram ) ) {
        legend$right <- list( fun = dendrogramGrob.fixed, args = list( signals.matrix.bihclust$proteins.ordering$elements.dgram, side='right' ) )
    }

    axis.fixed <- function (side = c("top", "bottom", "left", "right"), scales, 
        components, as.table, labels = c("default", "yes", "no"), 
        ticks = c("default", "yes", "no"), ...) {
        #panel.axis( side = side, at = components[side]$ticks$at, tck = components[side]$ticks$tck, draw.labels = FALSE, ... )
        if ( side %in% c('right','top') && !is.logical(components[[side]]) ) {
            #print( components[[side]] )
            panel.axis( side = side,
                        outside = TRUE,
                        at = components[[side]]$labels$at, 
                        labels = components[[side]]$labels$labels,
                        rot = c( 90, 0 ),
                        draw.labels = TRUE, ticks = FALSE )
        }
        axis.default( side, scales, components, as.table, 'no', 
            ifelse( side == 'top', 'yes', ticks ), ... )
    }

    # prepare cells matrix
    block.matrix <- matrix( NA, nrow = nrow( proteins ), ncol = nrow( samples ) )
    colnames(block.matrix) <- samples$sample
    rownames(block.matrix) <- proteins$protein_ac
    if ( is.null( measurements ) ) {
        cells.on.matrix <- NULL
        cells.off.matrix <- NULL
    } else {
        cells.on.matrix <- ms_data.matrix( measurements, protein.info,
                cellfunc = cell.func,
                msrun_column = 'msrun', bait_column = 'bait_ac', prey_column = 'prey_ac' )
        cells.on.matrix <- cells.on.matrix[ proteins$protein_ac, samples$sample ]
        cells.off.matrix <- block.matrix
    }
    #print(signals.matrix)
    for( protsClu in proteins.clusters$serial ) {
        clu.prots <- subset( bimap.props$proteins.clusters, proteins.cluster == protsClu )$protein_ac
        for( samplesClu in samples.clusters$serial ) {
            clu.samples <- subset( samples, samples.cluster == samplesClu )$sample
            clu.val <- signals.matrix[ protsClu, samplesClu ]
            block.matrix[ clu.prots, clu.samples ] <- signals.matrix[ protsClu, samplesClu ]
            if ( !is.null( measurements ) ) {
                # move cells data to off cells
                if ( is.na( clu.val ) ) {
                    cells.off.matrix[ clu.prots, clu.samples ] <- cells.on.matrix[ clu.prots, clu.samples ]
                    cells.on.matrix[ clu.prots, clu.samples ] <- NA
                }
            }
        }
    }
    #print( block.matrix )

    layout.widths = trellis.par.get('layout.widths')
    layout.widths$left.padding = 0
    layout.widths$right.padding = 0
    layout.widths$axis.right = protein.label.width
    layout.widths$axis.left = 0.5

    layout.heights = trellis.par.get('layout.heights')
    layout.heights$bottom.padding = 0
    layout.heights$axis.top = sample.label.width
    layout.heights$axis.bottom = 0.5

    res <- levelplot( t(block.matrix),
       col.regions = col,
       colorkey = list( space='bottom' ),
       column.values = ( 1:nrow(block.matrix)-0.5 ),
       row.values = 1:ncol(block.matrix)-0.5,
       xlab = 'baits/samples', ylab = 'proteins',
       panel = function(x,y,z) {
           #print('levelplot')
           panel.levelplot( x,y,z, subscripts=1:length(x), 
                   col.regions = col, alpha.regions = ifelse( !is.null(cells.on.matrix), 0.3, 1.0 ),
                   as.table = TRUE,
#                column.values = 1:nrow(block.matrix)-0.5, 
#                row.values = 1:ncol(block.matrix )-0.5
           )
           if ( !is.null( cells.on.matrix ) ) {
               panel.levelplot( x,y,t(cells.on.matrix), subscripts=1:length(x), 
                   col.regions = cells.col,
                   as.table = TRUE,
    #               column.values = 1:nrow(block.matrix)-0.5, 
    #               row.values = 1:ncol(block.matrix )-0.5
              )
           }
           if ( !is.null( cells.off.matrix ) &&
                any( !is.na( cells.off.matrix  ) )
           ){
               panel.levelplot( x,y,t(cells.off.matrix), subscripts=1:length(x), 
                   col.regions = cells.col,
                   alpha.regions = cells.off.alpha,
                   as.table = TRUE,
    #                column.values = 1:nrow(block.matrix)-0.5, 
    #                row.values = 1:ncol(block.matrix )-0.5
               )
           }
           # borders between biclusters
           #print('refline')
           if ( !is.na( grid.col ) ) {
               panel.refline( h = cumsum(proteins.clusters[1:nrow(proteins.clusters)-1,'size']),
                              v = cumsum(samples.clusters[1:nrow(samples.clusters)-1,'size']),
                             lty = grid.lty, lwd = grid.lwd, col = grid.col )
           }
           pushViewport( viewport( layout = grid.layout( 
                         nrow = nrow( signals.matrix ), ncol = ncol( signals.matrix ),
                         widths = samples.clusters$size, 
                         heights = rev( proteins.clusters$size ) ) ) )
#           if ( !is.null( measurements ) ) {
#               #print('measurements')
#               for ( ocId in rownames( signals.matrix ) ) {
#                   cluster.proteins <- subset( proteins, protein_ac %in% subset( bimap.props$proteins.clusters, proteins.cluster == ocId )$protein_ac )$protein_ac
#                   rowIx <- which( rev( rownames( signals.matrix ) ) == ocId )
#                   for ( scId in colnames( signals.matrix ) ) {
#                       cluster.samples <- subset( samples, sample %in% subset( bimap.props$samples.clusters, samples.cluster == scId )$sample )$sample
#                       colIx <- which( colnames( signals.matrix ) == scId )
#                       pushViewport( viewport( layout.pos.col = colIx, layout.pos.row = rowIx,
#                                     layout = grid.layout( nrow = 1, ncol = length( cluster.samples ) ) ) )
#                       for ( sample in cluster.samples ) {
#                           plot.crosscluster.measurements(
#                               cluster.proteins, sample, cluster.samples,
#                               msrun.info, signals.matrix, measurements )
#                       }
#                       popViewport()
#                   }
#               }
#           }
           if ( nrow( bimap.props$cross.clusters ) > 0 ) {
               #print('cross-clu')
               cross.clusters <- bimap.props$cross.clusters
               if ( show.borders ) {
                   cross.clusters <- ddply( cross.clusters, .( samples.cluster ), function( preys.clusters ) {
                       nested.clusters <- unique( preys.clusters$nested.preys.cluster )
                       colors <- rainbow( length( nested.clusters ) )
                       names( colors ) <- as.character( nested.clusters )
                       preys.clusters$nested.color <- colors[ as.character( preys.clusters$nested.preys.cluster ) ]
                       preys.clusters$pos <- sapply( preys.clusters$proteins.cluster, function(serial) {
                           which( proteins.clusters$serial == serial ) }
                       )
                       preys.clusters$nested.cluster.above <- lapply( preys.clusters$pos, function(cur.pos) {
                           subset( preys.clusters, pos == cur.pos - 1 )$nested.preys.cluster }
                       )
                       preys.clusters$nested.cluster.below <- lapply( preys.clusters$pos, function(cur.pos) {
                           subset( preys.clusters, pos == cur.pos + 1 )$nested.preys.cluster }
                       )
                       preys.clusters$border.above <- preys.clusters$nested.preys.cluster != preys.clusters$nested.cluster.above
                       preys.clusters$border.below <- preys.clusters$nested.preys.cluster != preys.clusters$nested.cluster.below
                       preys.clusters$border.above <- ifelse( is.na(preys.clusters$border.above ), TRUE, preys.clusters$border.above )
                       preys.clusters$border.below <- ifelse( is.na(preys.clusters$border.below ), TRUE, preys.clusters$border.below )
                       print(preys.clusters)
                       return ( preys.clusters )
                   } )
               }
               if ( plot.samples ) signals.subframe <- bimap.props$signals.subframe
               else signals.subframe <- NULL
               for ( ccIx in 1:nrow(cross.clusters) ) {
                  cc <- cross.clusters[ccIx, ]
                  protClu <- cc$proteins.cluster
                  sampleClu <- cc$samples.cluster
                  cluster.proteins <- subset( proteins, protein_ac %in% subset( bimap.props$proteins.clusters, proteins.cluster == protClu )$protein_ac )
                  cluster.samples <- subset( samples, samples.cluster == sampleClu )
                  pushViewport( viewport( layout.pos.row = which( rev( rownames( signals.matrix ) ) == protClu ), 
                                          layout.pos.col = which( colnames( signals.matrix ) == sampleClu ) ) )
                  #print(cc)
                  plot.crosscluster( cluster.proteins, cluster.samples,
                      protClu, sampleClu,
                      show.abundance = show.abundance,
                      signals.matrix, signals.matrix.sd,
                      signals.subframe, bait.col = bait.col,
                      border.col = cc$nested.color,
                      border.above = is.logical(cc$border.above) & cc$border.above,
                      border.below = is.logical(cc$border.below) & cc$border.below
                  )
                  popViewport()
              }
          }
          popViewport()
       },
       legend = legend,
       #scales = list( x = list( at = seq_len( colnames(signals.matrix) ), 
       #                         labels = colnames(signals.matrix) ) ),
       main = title,
       axis = axis.fixed,
       as.table = TRUE,
       reverse.rows = TRUE,
       aspect = aspect,
       par.settings = list( layout.widths = layout.widths, layout.heights = layout.heights ),
       xscale.components = function( lim, ... )
       {
           #print('xaxis')
           res <- xscale.components.default( lim, ...)
           res$num.limit = c(0, nrow(samples) )
           res$top <- levelplot.axis.components( samples, samples.clusters,
                        label_col = 'axis_label', small.ticks = 0.5, big.ticks = 5 )
           res$bottom <- levelplot.axis.components( samples, samples.clusters,
                         label_col = NA, small.ticks = 0.5, big.ticks = 1 )
           return ( res )
       },
       yscale.components = function( lim, ... )
       {
           #print('yaxis')
           res <- yscale.components.default( lim, ...)
           res$num.limit = c(0, nrow(proteins) )
           res$left <- levelplot.axis.components( proteins, proteins.clusters,
                            label_col = NA, small.ticks = 0.5, big.ticks = 1 )
           res$right <- levelplot.axis.components( proteins, proteins.clusters,
                            label_col = 'axis_label', small.ticks = 0.5, big.ticks = 5 )
           return ( res )
       } #)
    )
    if ( extended.result ) {
        return ( c( list( plot = res ), signals.matrix.bihclust ) )
    } else {
        print( res )
    }
}

BIMAP.plot_bimap <- function( bimap.walk, bimapId,
    protein.info, sample.info, msrun.info,
    cross_cluster.threshold = 0.6, ... )
{
    bimap.props <- BIMAP.extract_clustering( bimap.walk, bimapId, TRUE,
        cross_cluster.threshold = cross_cluster.threshold )
    plot.bimap( bimap.props, protein.info, sample.info, msrun.info,
            title = paste( 'Chessboard biclustering #', bimapId, sep = '' ), ... )
}

plot_ms_data <- function( matrix )
{
    
    plot( matrix )
}

plot_cluster_samples <- function( clusters.sample.frame, cluster.serial )
{
    cluster.sample.frame <- subset( clusters.sample.frame, cluster == cluster.serial )
    samples <- unique( cluster.sample.frame$sample )
    samples.lim <- c( min( cluster.sample.frame$signal ), max( cluster.sample.frame$signal ) )
    sample.cols <- rainbow( length( samples ) )
    names( sample.cols ) <- samples
    firstState <- TRUE
    for ( sp in samples ) {
        sample.frame <- subset( cluster.sample.frame, sample == sp )
        sample.steps <- sample.frame$step
        sample.signals <- sample.frame$signal
        if ( firstState ) {
            plot( sample.steps, sample.signals, col = sample.cols[[sp]], type = 'l', ylim = samples.lim )
            title( paste('Cluster',cluster.serial,'of proteins','...' ) )
            firstState <- FALSE
        }
        else {
            lines( sample.steps, sample.signals, col = sample.cols[[sp]] )
        }
    }
    legend( "bottomright", samples, pch=rep('-', length(samples) ), col=sample.cols )
}

iteration.ranges <- function( it.subsubset, it.subset, it.set )
{
    #print( it.set )
    #print( it.subset )
    set.ranges <- ranges( which( it.set %in% it.subsubset ) )
    subset.ranges <- lapply( set.ranges, function( range ) { 
            return ( c( which( it.subset == it.set[[ range[[1]] ]] ), 
                which( it.subset == it.set[[ range[[2]] ]] ) ) ) } )
    return ( subset.ranges )
}

plot.cell.signal <- function( bimap.walk, cell.object, cell.state, all.steps = bimap.walk@clusterings.walk$step )
{
    cell.proteins.clusters <- unique( subset( bimap.walk@proteins.clusters, as.character( object ) == cell.object )$proteins.cluster )

    subframe <- subset( bimap.walk@signals, proteins.cluster %in% cell.proteins.clusters & (state == cell.state) )
    if ( nrow( subframe ) == 0 ) {
        warning( paste('Protein', cell.object,'not found in bait/PD',cell.state) )
        return()
    }
    subframe$step <- as.numeric( as.character( subframe$step ) )
    col <- rainbow( length( unique( subframe$proteins.cluster ) ) )
    names( col ) <- unique( subframe$proteins.cluster )
    all.signal.steps <- unique( subframe$step )
    ix.ranges <- lapply( names( col ), function ( objclu.serial ) {
        its <- as.numeric( unique( subset( subframe, proteins.cluster == objclu.serial )$step ) )
        return ( iteration.ranges( its, all.signal.steps, all.steps ) )
    } )
    ix.ranges <- unlist( ix.ranges[ sapply( ix.ranges, function( l ) return ( length( l ) > 0 ) ) ], 
                         recursive = FALSE )
    meanSignal <- mean(subframe$signal)
    plot( c( min(subframe$step), max(subframe$step)), 
          c( meanSignal, meanSignal ), 
        ylab = paste( 'Signal of', cell.object, cell.state ), 
        ylim = c( min(subframe$signal), max(subframe$signal) ), type='l', col = "red", lwd = 2 )
    sdcol <- rgb( 100, 100, 100, 70, maxColorValue = 255 )
    for ( ix.range in ix.ranges ) {
        subsubframe <- subframe[ ix.range[[1]]:ix.range[[2]], ]
        it.range <- subframe$step[ ix.range ]
        subMean <- mean(subsubframe$signal)
        subSd <- sd( subsubframe$signal )
        rect( it.range[[1]], subMean - subSd, it.range[[2]], subMean + subSd, col=sdcol, border=NA, lty=3 )
        points( subsubframe$step, subsubframe$signal,
                col = col[ as.character( subsubframe$proteins.cluster ) ],
                pch = 19, cex = 0.5 )
        lines( it.range, c( subMean, subMean ), col="orange" )
        rug( side = 2, subsubframe$signal, 
             col = col[ as.character( subsubframe$proteins.cluster ) ] )
    }
}

plot.object.partitions <- function( bimap.walk, n = 10 )
{
    ops.frame <- merge( bimap.walk@clusterings.walk, bimap.walk@clusterings )#[,c('clustering.serial','objects.partition.serial')]
    clusterings.frame.topn <- top.factor.values( ops.frame, 'objects.partition.serial', n = n )
    ggplot( clusterings.frame.topn, 
            aes( factor(objects.partition.serial), fill=factor(clustering.serial)) ) +
        geom_bar()
}

BIMAP.interactors.walk <- function( bimap.walk, protein.ac )
{
    object.clusters <- subset( bimap.walk@objects.clusters, object == protein.ac )$objects.cluster
    proteins.clusterings.walk <- subset( merge( bimap.walk@clusterings.walk, bimap.walk@clusterings ), select = c('step', 'objects.partition.serial') ) 
    protein.clusters.walk <- subset( merge( proteins.clusterings.walk, bimap.walk@objects.partitions, by = 'objects.partition.serial' ), 
            objects.cluster.serial %in% object.clusters, select = c( 'step', 'objects.cluster.serial' ) )
    if ( nrow(protein.clusters.walk) != nrow(bimap.walk@clusterings.walk) ) warning("walk sizes do not match")
    res <- subset( merge( bimap.walk@objects.clusters, protein.clusters.walk, by = c('objects.cluster.serial') ), 
                    select = c( 'step', 'object', 'objects.cluster.serial' ) )
    return ( res )
}

plot.BIMAP.interactors <- function( bimap.walk, protein.ac,
                                     protein.info = NULL, label_col = 'short_id', description_col = 'description',
                                     min.freq = 0.01, min.cluster.freq = 0.1,
                                     drop.infrequent = FALSE,
                                     show.proteins.clusters = TRUE, show.protein.acs = FALSE )
{
    iact.walk <- BIMAP.interactors.walk( bimap.walk, protein.ac )
    if ( nrow( iact.walk ) == 0 ) {
        warning('No interactors for ', protein.ac)
    }
    protein.name <- function( object ) {
        ifelse( object != 'other', protein.info[ object, label_col ], object )
    }
    protein.label <- function( object ) {
            print( object )
        if ( show.protein.acs ) {
            paste( object , '\n(', protein.name( object ), ')', sep = '' )
        } else {
            protein.name( object )
        }
    }
    iact.walk <- within( iact.walk, {
                #print( protein.info[ as.character( object ), label_col ] )
                object <- factor( object )
                protein <- factor( object, levels=levels(object), labels = protein.label( as.character( levels( object ) ) ) )
                protein <- reorder( protein, protein, length )
                #print( levels( object ) )
    } )
    if ( !is.null( min.freq ) ) {
        iact.freq.table <- sort( table( iact.walk$object ) )
        freq_thresh <- min.freq * length( unique( iact.walk$step ) )
        print( freq_thresh )
        objects.filter <- names( iact.freq.table )[ iact.freq.table >= freq_thresh ]
        iact.walk$protein <- factor( ifelse( iact.walk$object %in% objects.filter, as.character( iact.walk$protein ), 'other' ) )
    }
    
    # clamp objects clusters to most frequent
    if ( !is.null( min.cluster.freq ) ) {
        objects.clusters.freq <- sort( table( as.character( unique( iact.walk[,c('step','objects.cluster.serial')])$objects.cluster.serial ) ), decreasing = TRUE )
        freq_thresh <- min.cluster.freq * length( unique( iact.walk$step ) )
        print( freq_thresh )
        objects.clusters.frequent <- names( objects.clusters.freq )[ objects.clusters.freq >= freq_thresh ]
        objects.cluster <- factor( ifelse( iact.walk$objects.cluster.serial %in% objects.clusters.frequent,
                                   as.character( iact.walk$objects.cluster.serial ),
                                   "other" ) )
    } else {
        objects.cluster <- factor( iact.walk$objects.cluster.serial )
    }
    print(levels(objects.cluster))
    proteins.cluster <- objects.cluster
    proteins.cluster.labels <- sapply( levels( proteins.cluster ), 
            function( serial ) {
                print(serial)
                if ( serial == 'other' ) { 
                    return ( 'other' )
                } else {
                    return ( paste( paste( unique( sapply( sort(unique( as.character( subset( iact.walk, objects.cluster == serial )$object ) ) ),
                                                           protein.name ) ),
                                collapse = '\n' ), '\n', sep='' ) ) 
                } } )
    #print(proteins.cluster.labels)
    iact.walk$objects.cluster <- objects.cluster
    iact.walk <- within( iact.walk, {
        proteins.cluster <- factor( proteins.cluster.labels[ proteins.cluster ] )
        #proteins.cluster <- reorder( proteins.cluster, proteins.cluster, length )
    } )
    
    title <- paste( "Interactors of ", protein.label( protein.ac ), sep = '' )
    if ( drop.infrequent & !is.null(min.freq ) ) iact.walk <- subset( iact.walk, object %in% objects.filter )
    if( show.proteins.clusters ) { 
        ggplot( iact.walk ) +
                geom_bar( aes( protein, fill = proteins.cluster ) ) +
                opts( title = title, axis.text.x=theme_text(angle=-90) ) +
                scale_x_discrete('interactors')
    } else {
        ggplot( iact.walk ) +
                geom_bar( aes( protein ) )+
                opts( title = title, axis.text.x=theme_text(angle=-90) ) +
                scale_x_discrete('interactors')
    }
}

