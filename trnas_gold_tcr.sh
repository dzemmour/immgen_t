#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-03:00:00
#SBATCH --mem 10G
#SBATCH -c 4
#SBATCH -o wrapper.log
#SBATCH -e wrapper.err
#SBATCH --mail-type=END,FAIL
# Run as sbatch /n/groups/cbdm_lab/immgen_t/immgent_wrapper_template_IGT32_40.sh .
# from Broad Run directory

### Load Relevant Modules



### Change to target directory
cd $1

source activate Blast

makeblastdb -in all_contig.fasta -out all_contig -parse_seqids -dbtype nucl
#magicblast -query all_contig.fasta -db trans_goldrath -out blast_0320_allcontigfasta.tsv -outfmt tabular -infmt fasta -no_unaligned -limit_lookup F
magicblast -query trans_goldrath.fasta -db all_contig -out blast_0320_allcontigfasta.tsv -outfmt tabular -infmt fasta -no_unaligned -limit_lookup F



#makeblastdb -in possorted_genome.fasta -out possorted_genome -parse_seqids -dbtype nucl
#magicblast -query trans_goldrath.fasta -db possorted_genome -out blast_output_tcr_bam.tsv -outfmt tabular -infmt fasta -no_unaligned -limit_lookup F

# Seurat Analysis
#echo "starting seurat analysis"

#Rscript $1/filter2.R combined_sample/outs/raw_feature_bc_matrix combined_sample/outs/filtered_feature_bc_matrix $2 $IGT Automatic $RNA_Cutoff $ADT_Cutoff > filter.log 2> filter.err



#Rscript $1/make_seurat_rna_adt_hto_noProj.R $project_description GEX/filtered_feature_bc_matrix FBC/filtered_feature_bc_matrix $hashes_to_remove $spleen  > make_seurat_rna_adt_hto.log 2> make_seurat_rna_adt_hto.err 

#th_nFeature_RNA_lo=300
#th_nFeature_RNA_hi=10000
#th_percent_mito=20


#Rscript $1/rna_qc_DM.R seuratobject_withHTOADT_singlet.Rds $th_nFeature_RNA_lo $th_nFeature_RNA_hi $th_percent_mito > rna_qc.log 2> rna_qc.err

#Rscript $1/singler.R seuratobject_withHTOADT_singlet_postRNAfiltering.Rds > singler.log 2> singler.err

#Rscript $1/adt_qc.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds $adt_hash_seqs 500 2 > adt_qc.log 2> adt_qc.err

#Rscript $1/Tcell_filt_part1.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds $genelist

#Rscript $1/naming_samples.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds $2 


## TCR Analysis
#Rscript $1/tcr_info.R $3 > tcr_info.log 2> tcr_info.err


#Rscript $1/Tcell_filt_part2_IGT32.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds $T_cutoffs $IGT $spleen  > T2.log 2> T2.err

#Rscript $1/dietSeurat.R seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering_postTfiltering.Rds seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds $2 $IGT

#Rscript $1/naming_samples2.R dataset.Rds dataset_clean.Rds $2 $IGT
#Rscript $1/Generating_QC_table.R seuratobject_withHTOADT_singlet_preRNAfiltering_QC_stats.txt seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_preADTfiltering.txt seuratobject_withHTOADT_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv seuratobject_IGT_singlet_postRNAfi
