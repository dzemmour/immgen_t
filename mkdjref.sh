#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-03:45:00
#SBATCH --mem 128G
#SBATCH -c 4
#SBATCH -o wrapper.log
#SBATCH -e wrapper.err
#SBATCH --mail-type=END,FAIL
# Run as sbatch /n/groups/cbdm_lab/immgen_t/immgent_wrapper_template_IGT32_40.sh /n/groups/cbdm_lab/immgen_t . IGT#
# from Broad Run directory

### Load Relevant Modules
module load bcl2fastq/2.19.1.403
module load fastqc/0.11.5
module load cellranger
module load gcc/9.2.0 R/4.1.2 hdf5/1.10.1 boost/1.62.0 openmpi/3.1.0 fftw/3.3.7 java/jdk-1.8u112 geos/3.10.2
export R_LIBS="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"
export R_LIBS_USER="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"


### Change to target directory
cd $1

fastafile=regions_IMGT.fa

cellranger mkvdjref --genome=my_vdj_ref \
                      --seq=$fastafile
