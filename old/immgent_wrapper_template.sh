#!/bin/bash
#SBATCH -c 2                    # Number of cores requested
#SBATCH -t 4-00:0:0                    # Runtime in minutes
                                # Or use HH:MM:SS or D-HH:MM:SS, instead of just number of minutes
#SBATCH -p priority                # Partition (queue) to submit to
#SBATCH --mem 32G        # 8 GB memory total
### In filenames, %j=jobid, %a=index in job array
#SBATCH -o 10X_HY2HTCCXY_wrapper.out               # Standard out goes to this file -o %j.out
#SBATCH -e 10X_HY2HTCCXY_wrapper.err               # Standard err goes to this file -e %j.err;
#SBATCH --mail-type=END         # Mail when the job ends  

##run as: sbatch -o wrapper_$(date +'%Y-%m-%d-%H-%M-%S').out -e wrapper_$(date +'%Y-%m-%d-%H-%M-%S').err 10x_wrapper.sh [runid] [layout file] [dir]
##Pipeline for Cell Rangers, multiple libraries on 1 flow cell (+- several lanes).
## Single pipeline output directory: 
###For cellranger mkfastq, the flowcell serial number is used (e.g., HAWT7ADXX)
###For cellranger count, aggr and reanalyze, the --id argument is used: outs/ contains the final pipeline output files.

##
MODIFY
runid
layoutfile
dir
fastqs

##
module load bcl2fastq/2.19.1.403
module load fastqc/0.11.5
module load cellranger
module load gcc/6.2.0 R/4.0.1 hdf5/1.10.1 boost/1.62.0 openmpi/3.1.0 fftw/3.3.7
export R_LIBS="/n/groups/cbdm_lab/dp133/R-4-0-1_libraries"
export R_LIBS_USER="/n/groups/cbdm_lab/dp133/R-4-0-1_libraries"
PATH_TO_SCRIPT=/n/groups/cbdm_lab/immgen_t


runid=BroadRun210809
layoutfile=broadrun210809_layoutfile.csv
dir=/n/groups/cbdm_lab/bv43/BroadRun210809
flowcell=HVFLNBGXJ

echo "Running cellranger mkfastq"

cellranger mkfastq --id=$runid --run=. --csv=$layoutfile


cellranger count --id=GEX \
				--sample=GEX \
				--fastqs=$dir/$runid/outs/fastq_path/$flowcell \
				--transcriptome=/n/groups/cbdm-db/bv43/M25/M25 \
				--nosecondary 

echo "starting cite-seq-count"
cd $dir/$runid/outs/fastq_path/$flowcell
mkdir combined

cat FBC_S3_L001_I1_001.fastq.gz FBC_S3_L002_I1_001.fastq.gz FBC_S3_L003_I1_001.fastq.gz FBC_S3_L004_I1_001.fastq.gz > $dir/$runid/outs/fastq_path/combined/FBC_S3_L1234_I1_001.fastq.gz
cat FBC_S3_L001_I2_001.fastq.gz FBC_S3_L002_I2_001.fastq.gz FBC_S3_L003_I2_001.fastq.gz FBC_S3_L004_I2_001.fastq.gz > $dir/$runid/outs/fastq_path/combined/FBC_S3_L1234_I2_001.fastq.gz
cat FBC_S3_L001_R1_001.fastq.gz FBC_S3_L002_R1_001.fastq.gz FBC_S3_L003_R1_001.fastq.gz FBC_S3_L004_R1_001.fastq.gz > $dir/$runid/outs/fastq_path/combined/FBC_S3_L1234_R1_001.fastq.gz
cat FBC_S3_L001_R2_001.fastq.gz FBC_S3_L002_R2_001.fastq.gz FBC_S3_L003_R2_001.fastq.gz FBC_S3_L004_R2_001.fastq.gz > $dir/$runid/outs/fastq_path/combined/FBC_S3_L1234_R2_001.fastq.gz

cd $dir

module load python/3.6.0
source /n/groups/cbdm_lab/bv43/python3/bin/activate

CITE-seq-Count -R1 $dir/$runid/outs/fastq_path/combined/FBC_S3_L1234_R2_001.fastq.gz \
		-R2 $dir/$runid/outs/fastq_path/combined/FBC_S3_L1234_R2_001.fastq.gz \
 		-t /n/groups/cbdm_lab/immgen_t/adt_hash_seq_file_April2021_panel_demux.csv -cells 100000 -cbf 1 -cbl 16 -umif 17 -umil 26 -o FBC --sliding-window

echo "Running cellranger vdj"

cellranger vdj --id=TCR \
                 --reference=/n/groups/cbdm_lab/bv43/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
                 --fastqs=$dir/$runid/outs/fastq_path/$flowcell \
                 --sample=TCR \



echo "starting seurat analysis"

PATH_TO_PROJECT=$dir

PROJECT=$PATH_TO_PROJECT/project_description.csv #adds experimental conditions and HT descriptionos

ADT_HASH_SEQ_FILE=$PATH_TO_SCRIPT/adt_hash_seq_file_April2021_panel.csv #only used in adt_qc.R !!!

ADT_HASH_TO_REMOVE_FOR_ANALYSIS_FILE=$PATH_TO_PROJECT/hash_adt_rows_to_remove_from_analysis.txt


cd $PATH_TO_PROJECT

Rscript $PATH_TO_SCRIPT/make_seurat_rna_adt_hto.R $PROJECT GEX/outs/filtered_feature_bc_matrix FBC/umi_count $ADT_HASH_TO_REMOVE_FOR_ANALYSIS_FILE > make_seurat_rna_adt_hto.log 2> make_seurat_rna_adt_hto.err

th_nFeature_RNA_lo=500
th_nFeature_RNA_hi=10000
th_percent_mito=20
Rscript $PATH_TO_SCRIPT/rna_qc.R $PATH_TO_PROJECT/seuratobject_withHTOADT_singlet.Rds $th_nFeature_RNA_lo $th_nFeature_RNA_hi $th_percent_mito > rna_qc.log 2> rna_qc.err
Rscript $PATH_TO_SCRIPT/singler.R seuratobject_withHTOADT_singlet_postRNAfiltering.Rds > singler.log 2> singler.err
Rscript $PATH_TO_SCRIPT/adt_qc.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds $ADT_HASH_SEQ_FILE 500 2 > adt_qc.log 2> adt_qc.err