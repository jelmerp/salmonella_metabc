# DIRS, FILES, AND SETTINGS ----------------------------------------------------
## Dirs
dir_fq_org=data/fastq_original/211210_Girish_Dhanashree_RunAll-313451170/
dir_fq_all=../data/fastq/all
dir_fq=data/fastq/
dir_fastqc=results/fastqc/
dir_multiqc=results/multiqc/
dir_trim=results/cutadapt/
dir_asv=results/dada
dir_ps=results/phyloseq
dir_tax=results/taxonomy

## Output files
seqtab="$dir_asv"/seqtab.rds
tree="$dir_asv"/tree.rds
tax_dada="$dir_tax"/taxa_dada.rds
ps_raw="$dir_ps"/ps_dadatax_raw.rds
ps_filt="$dir_ps"/ps_dadatax_filt.rds

## Input files and settings
meta=metadata/meta.tsv
prim_f="GAGTGCCAGCMGCCGCGGTAA"
prim_r="ACGGACTACHVGGGTWTCTAAT"

config_dada=workflow/config/config_dada.R
config_ps=workflow/config/config_ps-filter.R

## Scripts
bin=mcic-scripts/metabarcoding


# WORKFLOW ---------------------------------------------------------------------
## Process metadata
sbatch scripts/meta.R metadata/meta_raw.tsv "$meta"

## Copy all FASTQ files into a single directory,
## then take subset belonging to the Salmonella project 
mkdir -p "$dir_fq_all" "$dir_fq"
find "$dir_fq_org" -name "*fastq.gz" -exec cp {} "$dir_fq_all" \;
cp "$dir_fq_all"/*_S24[1-9]*.gz \
    "$dir_fq_all"/*_S25[0-9]*.gz  \
    "$dir_fq_all"/*_S26[0-7]*.gz  \
    "$dir_fq_all"/*_S29[1-9]*gz \
    "$dir_fq_all"/*_S30[0-3]*gz \
    "$dir_fq"

## Run FastQC
for fq in "$dir_fq"/*fastq.gz; do
    sbatch mcic-scripts/qc/fastqc.sh -i "$fq" -o "$dir_fastqc"
done

## Run cutadapt
for fq in "$dir_fq"/*_R1_*fastq.gz; do
    sbatch "$bin"/cutadapt.sh -i "$fq" -o "$dir_trim" -f "$prim_f" -r "$prim_r"
done

## Run MultiQC
sbatch mcic-scripts/qc/multiqc.sh -i results -o "$dir_multiqc"

## ASV inference w/ dada
sbatch --mem=32G -c 8 -t 60 "$bin"/dada.R -i "$dir_trim" -o "$dir_asv" -c "$config_dada"
sbatch "$bin"/dada_qc.R -i "$dir_asv"/qc/nseq_summary.txt -o "$dir_asv"/qc

## Assign taxonomy, construct tree, and make phyloseq object
sbatch "$bin"/tax_assign.sh -i "$seqtab" -o "$tax_dada" -a dada
sbatch "$bin"/tree_build.R -i "$seqtab" -o "$tree"
sbatch "$bin"/ps_make.R -s "$seqtab" -x "$tax_dada" -t "$tree" -m "$meta" -o "$ps_raw"
sbatch "$bin"/ps_filter.R -i "$ps_raw" -o "$ps_filt" -q "$dir_ps"/qc -c "$config_ps"

sbatch "$bin"/ps_agglomtaxa.R -i "$ps_filt" -o "$dir_ps"/glom
