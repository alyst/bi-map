#! /bin/sh
# Run BI-MAP ON TIP49ab dataset

BIMAP_PATH=../build/release/src
DATA_PATH=tip49

$BIMAP_PATH/BIMAP-sampler --proteins_file=$DATA_PATH/proteins.txt \
                          --exp_design_file=$DATA_PATH/exp_design.txt \
                          --measurements_file=$DATA_PATH/measurements.txt \
                          --config_file=bimap_config.ini \
                          --output_file=tip49_walk.xml.bz2 \
                          --samples=10 \
                          --eesampler.burnin_iterations=10
                   

