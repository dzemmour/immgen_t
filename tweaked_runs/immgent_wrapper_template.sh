#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-03:00:00
#SBATCH --mem 64G
#SBATCH -c 16
#SBATCH -o wrapper.log
#SBATCH -e wrapper.err
#SBATCH --mail-type=END,FAIL
# Run as sbatch /n/groups/cbdm_lab/immgen_t/immgent_wrapper_template.sh /n/groups/cbdm_lab/immgen_t . [IGT Code]
# from Broad Run directory

### Load Relevant Modules
module load bcl2fastq/2.19.1.403
module load fastqc/0.11.5
module load cellranger
module load gcc/9.2.0 R/4.1.2 hdf5/1.10.1 boost/1.62.0 openmpi/3.1.0 fftw/3.3.7 java/jdk-1.8u112 geos/3.10.2 
export R_LIBS="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library"
export R_LIBS_USER="/n/groups/cbdm-db/ly82/scRNAseq_DB/R-4.1.2/library" 

### Change to target directory
cd $2

### Make sure to move and modify these files to the target directory
# Change Layout File 
layoutfile=layoutfile.csv
# FASTQ information
libraries=libraries.csv
# Hashtag file
feature_refs=feature_refs_2022panel.csv
# adds experimental conditions and HT descriptions; Taken from Ian's Single-Cell RNA Sequencing spreadsheet
project_description=project_description.csv 
# Hashes that are not used
hashes_to_remove=hash_adt_rows_to_remove_from_analysis.txt
# Information for various antibody-derived tags, only used in adt_qc.R !!!
adt_hash_seqs=adt_hash_seq_file_2022_panel.csv 
#sample names per Hashtag
samples=samples.csv
#Location of BCL run  folder
bcl_location=$(<bcl_location.csv)

# CellRanger Pipeline
echo "Running cellranger mkfastq"
cellranger mkfastq --id=fastqs --run=..$bcl_location --csv=$layoutfile 
# --use-bases-mask=Y26n*,I10,N10,Y90n* --barcode-mismatches=2

echo "Running cellranger count"
cellranger count  --id=combined_sample \
                   --libraries=$libraries \
                   --feature-ref=$feature_refs \
                   --transcriptome=/n/groups/cbdm-db/bv43/M25/M25 \
                   --expect-cells=100000

### Script to split output combined_sample into GEX and FBC should go here
mkdir -p GEX
mkdir -p FBC

### CellRanger Pipeline (TCR)
echo "Running cellranger vdj"
cellranger vdj --id=TCR \
                 --reference=/n/groups/cbdm_lab/bv43/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
                 --fastqs=fastqs/outs/fastq_path/ \
                 --sample=TCR

# Seurat Analysis
echo "starting seurat analysis"

Rscript $1/filter.R combined_sample/outs/raw_feature_bc_matrix combined_sample/outs/filtered_feature_bc_matrix $2

Rscript $1/make_seurat_rna_adt_hto.R $project_description GEX/filtered_feature_bc_matrix FBC/filtered_feature_bc_matrix $hashes_to_remove > make_seurat_rna_adt_hto.log 2> make_seurat_rna_adt_hto.err

th_nFeature_RNA_lo=500
th_nFeature_RNA_hi=10000
th_percent_mito=20

Rscript $1/rna_qc.R seuratobject_withHTOADT_singlet.Rds $th_nFeature_RNA_lo $th_nFeature_RNA_hi $th_percent_mito > rna_qc.log 2> rna_qc.err

Rscript $1/singler.R seuratobject_withHTOADT_singlet_postRNAfiltering.Rds > singler.log 2> singler.err

Rscript $1/naming_samples.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR $2

Rscript $1/adt_qc.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds $adt_hash_seqs 500 2 > adt_qc.log 2> adt_qc.err

# Rscript $1/adt_qc_distribution.R FBC/outs GEX/outs > adt_qc_distribution.log 2> adt_qc_distribution.err

## TCR Analysis
Rscript $1/tcr_info.R $3 > tcr_info.log 2> tcr_info.err

# ## MiXCR Analysis
# mkdir -p fastqs/outs/fastq_path/combined_tcr
# mkdir -p MiXCR

# cwd=$(pwd)
# flowcell=${cwd: -9}
# cat fastqs/outs/fastq_path/$flowcell/TCR*R1* > fastqs/outs/fastq_path/combined_tcr/TCR_R1.fastq.gz
# cat fastqs/outs/fastq_path/$flowcell/TCR*I1* > fastqs/outs/fastq_path/combined_tcr/TCR_I1.fastq.gz
# cat fastqs/outs/fastq_path/$flowcell/TCR*R2* > fastqs/outs/fastq_path/combined_tcr/TCR_R2.fastq.gz

# Rscript $1/MiXCR_prep.R > mixcr.log 2> mixcr.err

# ulimit -n 65536
# ulimit -u 65536

# module load java/jdk-11.0.11
# module load python/3.7.4 
# source /n/groups/cbdm_lab/nip330/python3_7/bin/activate

# cd MiXCR
# for ht in */
# do
# 	echo $ht
#     cd $ht
#     mkdir -p umicollapse
#     mkdir -p vdjca
#     mkdir -p clna
#     mkdir -p clones

#     barcode_splitter --bcfile barcodes.txt ../../fastqs/outs/fastq_path/combined_tcr/TCR_R1.fastq.gz ../../fastqs/outs/fastq_path/combined_tcr/TCR_R2.fastq.gz ../../fastqs/outs/fastq_path/combined_tcr/TCR_I1.fastq.gz --idxread 1 --suffix .fastq
#     rm unmatched-read-*
#     rm multimatched-read-*

#     for fq in *2.fastq
#     do
#     	file_name=${fq%%-*}
#         java -jar /n/groups/cbdm_lab/nip330/UMICollapse/umicollapse.jar fastq -i $fq -o umicollapse/${file_name}_umicollapse.fastq --merge avgqual
#     	/n/groups/cbdm_lab/nip330/mixcr-3.0.13/mixcr align --library imgt --species mmu umicollapse/${file_name}_umicollapse.fastq vdjca/${file_name}.vdjca
#         /n/groups/cbdm_lab/nip330/mixcr-3.0.13/mixcr assemble -r clna/${file_name}_report.txt -a vdjca/${file_name}.vdjca clna/${file_name}.clna
#         /n/groups/cbdm_lab/nip330/mixcr-3.0.13/mixcr exportClones --preset full clna/${file_name}.clna clones/${file_name}_clones.txt
#     done


#     cd ..
# done

# cd ..
# module load java/jdk-1.8u112

# Rscript $1/MiXCR_Formatting.R > mixcr_format.log 2> mixcr_format.err
