# Serrano Remining R Docker

A plain R Docker image for running R scripts in this repository.
Example:

    docker build -t serrano-remining-rbase .

    docker run -it --rm -v $PWD:/src:ro -v $PWD/scratch:/scratch \
        serrano-remining-rbase /src/test_HDF5.R /src/Primary_146_81_features.h5
