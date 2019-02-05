#!/bin/bash

module load R/3.4.2-fasrc01
module load gcc/7.1.0-fasrc01
module load bzip2/1.0.6-fasrc01
module load hdf5/1.10.1-fasrc03

Rcode=$1

Rscript $Rcode
