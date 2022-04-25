#!/bin/bash
## run as bsub -n 1 -M 4 -W 16:00 ./runGistic.sh
#BSUB -n 1 -M 12 -W 359
#<usage>
[[ $# -gt 0 ]] || {
    echo "Description:"
    echo "This script to runs GISTIC for 20kb with parameters optimized for bulkDNA log-ratio cutoffs"
    echo ""
    echo "Usage:"    
    echo "bsub -n 1 -M 12 -W 359 ./src/05_runGISTIC_SC.sh <confidence> path/to/file.cbs.seg <file.extension> <5k|20k|50k bin_level> <out_dir>"
    echo "bsub -n 1 -M 12 -W 359 ./src/05_runGISTIC_SC.sh 80 path/to/file.cbs.seg 20k .cbs.seg gistic_out/"
    echo ""
    exit 1;
}
#</usage>

set -e -x -o pipefail -u

## path to GISTIC hg19.mat file is located
HG19MAT=${HOME}/genomes/homo_sapiens/Ensembl/GRCh37.p13/Sequence/GISTIC2/hg19.mat


## confidence intervals, seg file, and file extension
CONF=$1
INFILE=$2
EXTENSION=$3
BIN_LEVEL=$4
## VBINS50K=res/gistic/grch37.bin.boundaries.50k.bowtie.k50.makerFile.txt
## VBINS20K=res/gistic/grch37.bin.boundaries.20k.bowtie.k50.makerFile.txt
## VBINS5K=res/gistic/grch37.bin.boundaries.5k.bowtie.k50.makerFile.txt

## path to marker files based on 5k, 20k, and 50k 
if [ $BIN_LEVEL == "5k" ]
then
    VBINS=res/gistic/grch37.bin.boundaries.5k.bowtie.k50.makerFile.txt
elif [ $BIN_LEVEL == "20k" ]
then
    VBINS=res/gistic/grch37.bin.boundaries.20k.bowtie.k50.makerFile.txt
elif [ $BIN_LEVEL == "50k" ]
then
    VBINS=res/gistic/grch37.bin.boundaries.50k.bowtie.k50.makerFile.txt
fi

## creating sample output diretory :: removes extension
[[ ! -z "$EXTENSION" ]] || EXTENSION=.seg
OUT_STEM=$( basename $INFILE $EXTENSION )
# output directory to collect all samples, for this directory use "./"
OUT_DIR=$5

## create the final path/to/output
BD=${OUT_DIR}/out_DNA_${OUT_STEM}_${BIN_LEVEL}_c${CONF}  ## out_[DNA|SC]_hg19_1M_c${CONF}
[[ -d $BD ]] || mkdir $BD

## DNA: Default thresholds are kept for amplifications and deletions
##  SC: Set thresholds for copy number to AMP = log2(2.5)/2, and DEL = log2(1.2)/2
~/bin/gistic2 -refgene ${HG19MAT} -b ${BD} -seg ${INFILE} -mk $VBINS \
	      -conf 0.${CONF} -broad 1 -twoside 1 \
	      -res 0.05 -genegistic 1 -rx 0 -js 8 \
	      -savegene 1 -gcm median -v 10 -armpeel 1


## __END__
## Usage: gp_gistic2_from_seg -b base_dir -seg segmentation_file
## -refgene ref_gene_file  [-mk markers_file] [-alf array_list_file(def:empty)]
## [-cnv cnv_file] [-ta amplifications_threshold(def=.1)] [-td deletions_threshold(def=.1)]
## [-js join_segment_size(def=8)] [-ext extension] [-qvt qv_thresh(def=0.25)]
## [-rx remove_x(def=1)] [-v verbosity_level(def=0)] [-cap cap_val(def=1.5]]
## [-broad run_broad_analysis(def=0)] [-brlen broad_length_cutoff(def=0.98)]
## [-maxseg max_sample_segs(def=2500)] [-res res(def=0.05)] [-conf conf_level(def=0.75)]
## [-genegistic do_gene_gistic(def=0)] [-smalldisk save_disk_space(def=0)]
## [-smallmem use_segarray(def=1)] [-savegene write_gene_files(def=0)]
## [-arb do_arbitration(def=1)] [-twosides use_two_sided(def=0)] [-peaktype peak_types(def=robust)]
## [-saveseg save_seg_data(def=1)] [-savedata write_data_files(def=1)]
## [-armpeel armpeel(def=1)] [-gcm gene_collapse_method(def=mean)]
## [-scent sample_center(def=median)] [-maxspace max_marker_spacing]
## [-logdat islog(def=auto-detect)]
