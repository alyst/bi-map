require( plyr )

peptide.sharingFactors <- function( protein_ms_data, peptide_ms_data,
        peptide_cols = c( 'Peptide' ),
        protein_cols = 'prey_ac_noiso',
        peptide_score_col = 'mascot_score',
        exp_cols = c( 'bait_ac_noiso',
                'ap_type', 'ms_type',
                'sample', 'sample_origin',
                'msrun', 'msrun_origin' )
){
    compound_column <- function( frame, col_names, sep = '|_|' ) {
        do.call( 'paste', c( lapply( col_names, function( col_name ) as.character(frame[[col_name]]) ), sep = sep ) )
    }
    group_by_cols <- function( frame, col_names ) {
        ids <- plyr:::id(frame[c(col_names)], drop = TRUE )
        max.groups = min( max( ids ), attr( ids, "n" ) )
        message( max.groups )
        indices <- plyr:::split_indices( seq_len(nrow(frame)), ids,
                n = max.groups )
        return ( list( id_cols = col_names,
                        ids = ids, indices = indices,
                        group_size = vapply( indices, length, integer(1) ),
                        group_ixs = match(seq_len(max.groups), ids)
                ) )
    }
    # only those peptides that belong to filtered proteins
    message( 'Filtering peptide for observed proteins' )
    protein_ms_data$prot2exp_id <- compound_column( protein_ms_data, c(protein_cols, exp_cols) )
    peptide_ms_data$prot2exp_id <- compound_column( peptide_ms_data, c(protein_cols, exp_cols) )
    peptide_ms_data.filtered <- subset( peptide_ms_data, prot2exp_id %in% protein_ms_data$prot2exp_id )
    
    
#    data <- as.matrix(combined.probs[ c( "prob.xclass",
#                                         paste( 'prob', org, sep = '.' ),
#                                         "prob.VxZ" )])
    
    message( 'Building Peptide2Protein2Experiment table' )
    pep2prot2exp.groups <- group_by_cols( peptide_ms_data.filtered, c( peptide_cols, 'prot2exp_id' ) )
    pep2prot2exp <- cbind( peptide_ms_data.filtered[ pep2prot2exp.groups$group_ixs, c( peptide_cols, protein_cols, exp_cols, 'prot2exp_id' ) ],
            sc = pep2prot2exp.groups$group_size,
            pc = 1
    )
    if ( is.character( peptide_score_col ) ) {
        pep2prot2exp$max.score = vapply( pep2prot2exp.groups$indices, function( rows ) max( c( 0, peptide_ms_data.filtered[ rows, peptide_score_col ] ), na.rm = TRUE ), numeric(1) )
    }
    message( 'Building Peptide2Experiment groups' )
    pep2exp.groups <- group_by_cols( pep2prot2exp, c( peptide_cols, exp_cols ) )
    message( 'Building Protein2Experiment groups' )
    print( colnames( pep2prot2exp ))
    prot2exp.groups <- group_by_cols( pep2prot2exp, 'prot2exp_id' )
    pep2prot2exp$pep2exp_id <- as.character( pep2exp.groups$ids )
    pep2prot2exp$prot2exp_id <- as.character( prot2exp.groups$ids )
    
    message( 'Building Peptide2Experiment table' )
    pep2exp <- cbind( pep2prot2exp[ pep2exp.groups$group_ixs, c( peptide_cols, exp_cols, 'pep2exp_id' ) ],
            sc = vapply( pep2exp.groups$indices, function( rows ) sum( pep2prot2exp[ rows, 'sc'], na.rm = TRUE ), integer(1) ),
            pc = vapply( pep2exp.groups$indices, function( rows ) integer(1), integer(1) ),
            nproteins = vapply( pep2exp.groups$indices,
                    function(rows) length( unique( pep2prot2exp[ rows, 'prot2exp_id' ] ) ),
                    integer(1) )
    )
    if ( 'max.score' %in% colnames( pep2prot2exp ) ) {
        pep2exp$max.score <- vapply( pep2exp.groups$indices, function( rows ) max( pep2prot2exp[ rows, 'max.score' ], na.rm = TRUE ), numeric(1) )
    }
    rownames( pep2exp ) <- pep2exp$pep2exp_id
    
    # sc_unique is zero if this peptide is shared
    pep2prot2exp$is_proteotypic <- pep2exp[ pep2prot2exp$pep2exp_id, 'nproteins' ] == 1 
    pep2prot2exp$sc_unique <- as.integer( ifelse( pep2prot2exp$is_proteotypic, pep2prot2exp$sc, 0 ) )
    pep2prot2exp$pc_unique <- as.integer( ifelse( pep2prot2exp$is_proteotypic, pep2prot2exp$pc, 0 ) )
    
    message( 'Building Protein2Experiment table' )
    prot2exp <- cbind( pep2prot2exp[ prot2exp.groups$group_ixs, c( protein_cols, exp_cols, 'prot2exp_id' ) ],
            sc_all = vapply( prot2exp.groups$indices,
                    function(rows) sum( pep2prot2exp[ rows, 'sc' ], na.rm = TRUE ),
                    integer(1) ),
            sc_unique = vapply( prot2exp.groups$indices,
                    function(rows) sum( pep2prot2exp[ rows, 'sc_unique' ], na.rm = TRUE ),
                    integer(1) ),
            pc_unique = vapply( prot2exp.groups$indices,
                    function(rows) sum( pep2prot2exp[ rows, 'pc_unique' ], na.rm = TRUE ),
                    integer(1) )
    )
    if ( 'max.score' %in% colnames( pep2prot2exp ) ) {
        prot2exp$max.score = vapply( prot2exp.groups$indices, function( rows ) sum( pep2prot2exp[ rows, 'max.score' ], na.rm = TRUE ), numeric(1) )
    }
    rownames( prot2exp ) <- prot2exp$prot2exp_id
    
    message( 'Calculating initial protein abundance' )
    ## Calculate initial estimates for reporter protein abundances
    ## use the sum of unique spectrum counts
    ## TODO: support for other initial abundance estimates
    prot2exp$ini.abundance <- prot2exp$sc_unique
    
    # calculate sharing factor for each protein that shares given peptide
    message( 'Calculating peptide sharing factors' )
    pep2prot2exp$prot.ini.abundance <- prot2exp[ pep2prot2exp$prot2exp_id, 'ini.abundance' ]
    sum.abundance <- vapply( pep2exp.groups$indices, function( pep_rows ) sum( pep2prot2exp[ pep_rows, 'prot.ini.abundance' ] ), numeric(1) )
    names( sum.abundance ) <- pep2exp$pep2exp_id
    pep2prot2exp$share_factor <- pep2prot2exp$prot.ini.abundance / sum.abundance[ pep2prot2exp$pep2exp_id ]
    
    message( 'Calculating adjusted protein SCs' )
    # calculate adjusted spectral counts as weighted sum of share factor and spectral counts
    prot2exp$pc <- vapply( prot2exp.groups$indices,
            function(rows) sum( pep2prot2exp[ rows, 'pc' ] ),
            numeric(1) )
    prot2exp$sc_adj <- vapply( prot2exp.groups$indices,
            function(rows) pep2prot2exp[ rows, 'sc' ] %*% pep2prot2exp[ rows, 'share_factor' ],
            numeric(1) )
    prot2exp$pc_adj <- vapply( prot2exp.groups$indices,
            function(rows) pep2prot2exp[ rows, 'pc' ] %*% pep2prot2exp[ rows, 'share_factor' ],
            numeric(1) )
    
    return ( list( protein_ms_data = prot2exp,
                   peptide_ms_data = pep2prot2exp ) )
}
