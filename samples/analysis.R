# Sample AP-MS data analysis script
# 
# Author: Alexey Stukalov (astukalov@cemm.oeaw.ac.at)
###############################################################################

# base path to R scripts and BI-MAP (update to actual location)
base_path <- "/home/<user>/projects/bi-map"

# path to BI-MAP scripts
bimap_scripts_path <- file.path( base_path, "scripts/R" )
# path to data files
data_path <- file.path( base_path, "samples/tip49" )
# path to output files
results_path <- file.path( base_path, "results" )

# define path to RBIMAP library (update to actual location)
RBIMAP.libpath <- file.path( base_path, "build/release/src/rcpp" )

# load R BIMAP functions
source( file.path( bimap_scripts_path, "BIMAP.R" ) )
source( file.path( bimap_scripts_path, "BIMAP_plot.R" ) )
source( file.path( bimap_scripts_path, "BIMAP_graphml.R" ) )
source( file.path( bimap_scripts_path, "BIMAP_xlsx.R" ) )

# load AP-MS data
bimap.data <- BIMAP.msdata.load( data_path, format = 'CSV' )

# add protein description and Ensemble Gene names to proteins info
proteins_info.ensembl <- read.table( file.path( data_path, 'ensembl_info.csv' ),
                                    sep = '\t', header = TRUE, stringsAsFactors = FALSE )
rownames( proteins_info.ensembl ) <- proteins_info.ensembl$gi_number

setdiff( bimap.data$proteins$protein_ac, proteins_info.ensembl$gi_number )

bimap.data.extra <- BIMAP.msdata.extra_info( bimap.data, proteins_info.ensembl )

# load calculated BI-MAP MCMC walk
# (run samples/run_bimap_tip49.sh to get run)
bimap.walk <- BIMAP.mcmcwalk.load( file.path( results_path, 'tip49_walk.xml.bz2' ),
    0.02, 0.02 )

# or execute BI-MAP

# generate a dataframe of top n clustering models sorting according to their posterior probability
bimap.best_clusterings.df <- BIMAP.mcmcwalk.biclusterings_order( bimap.walk, n = 10 )
# extract the model with max posterior probability
bimap.best_clustering <- BIMAP.mcmcwalk.extract_biclustering( bimap.walk, bimap.best_clusterings.df[1,'clustering.serial'] )

# render the matrix view of BI-MAP model to PDF file
pdf( file= file.path( results_path, 'tip49_blocks.pdf' ),
     title="TIP49a/b", width=34, height=60 )
BIMAP.plot( bimap.best_clustering,
    bimap.data,
    show.abundance.labels = FALSE,
    show.protein_ac = TRUE, show.sample_id = TRUE, show.msruns = TRUE,
    show.measurements = TRUE,
    aspect = 3.0
)
dev.off()

# export to excel
bimap.xlsx <- BIMAP.create_xlsx( bimap.best_clustering,
    bimap.data,
    show.msruns = TRUE,
    show.measurements = TRUE
)
saveWorkbook( bimap.xlsx, file.path( results_path, 'tip49_blocks.xlsx' ) )

# construct BI-MAP network representation in GraphML format
bimap.graphml <- BIMAP.graphML( bimap.best_clustering, bimap.data )
bimap.unstyled.graphml.filename <- file.path( results_path, 'tip49_bimap.graphml' )
saveXML( bimap.graphml, file = bimap.unstyled.graphml.filename )
# apply yEd-compatible visual scheme to the graph
system( paste( "xsltproc -o", file.path( results_path, 'tip49_bimap_y.graphml' ),
        file.path( base_path, 'scripts/xslt/yfilize.xslt' ),
        bimap.unstyled.graphml.filename ) )
# the file tip49_bimap_y.graphml could now be loaded into yEd 3.8+ for layouting
