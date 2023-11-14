#!/bin/bash

# Simulate cross-sample contamination by creating fastq with x% reads from sample A and y% from sample by
# Strategy: 
# If simulating 5% contamination
# 1. Create a new "Sample C" with 95% reads from A, and 5% from B.
# 2. If total reads 80M
# Total reads:  80M reads (P23993 has 92.5M reads)
# Reads from A: 76M
# Reads from B: 4M

# Get input
usage() { echo "Usage: $0 [-p <percentage contamination to simulate>] [-r <total M reads in output fastq>] [config file]" 1>&2; exit 1; }

while getopts ":p:r:" o; do
    case "${o}" in
        p)
            perc=${OPTARG}
            ;;
        r)
            totreads=$(echo "${OPTARG} * 1000000" | bc )
            ;;
        *)
            usage
            ;;
    esac
done
shift $((OPTIND-1))

if [ -z "${perc}" ] || [ -z "${totreads}" ]; then
    echo "Error: Please enter -p and -r "
    usage
    exit
fi

if (( $(echo "${perc} < 0 " |bc -l) )) || (( $(echo "${totreads} < 0" | bc -l) )); then
    usage
    exit
fi


set_params() {
    outdir=${PWD}/cont${perc}_reads${totreads}
    mkdir -p $outdir

    # Get percentage of A and B in final file
    # A is "true" sample
    # B is "contaminating" sample
    percA=$(echo "100 - $perc" | bc)
    percB=$perc

    # Get number of reads that will make those percentages
    readsA=$(echo "$totreads * $percA / 100" | bc)
    readsB=$(echo "$totreads * $percB / 100" | bc)
}

print_params() {
    echo "Total reads:      $totreads"
    echo "% contamination:  $perc %"
    echo "-------"
    echo "% A reads:        $percA %"
    echo "nr reads A:       $readsA"
    echo "-------"
    echo "% B reads:        $percB %"
    echo "nr reads B:       $readsB"
}

timeCMD() {
    start=$(date +%s)
    $1
    end=$(date +%s)
    t=$(($end - $start))
    printf 'Time taken: %d minutes, %d seconds\n' "$(( t/60 ))" "$(( t - 60*(t/60) ))"
}

process_sample() {
    echo ""
    echo "--------------------------------------------"
    echo "=> processing $1 Sample A: "
    echo "- Downsampling $A1 to $readsA reads.. "
    timeCMD "seqkit sample -s11 -2 -n $readsA -o $outdir/A.$1.R1.tmp.fq.gz $A1"
    echo "- Downsampling $A2 to $readsA reads.. "
    timeCMD "seqkit sample -s11 -2 -n $readsA -o $outdir/A.$1.R2.tmp.fq.gz $A2"

    echo ""
    echo "--------------------------------------------"
    echo "=> processing $1 Sample B: "
    echo "- Downsampling $B1 to $readsB reads.. "
    timeCMD "seqkit sample -s11 -2 -n $readsB -o $outdir/B.$1.R1.tmp.fq.gz $B1"
    echo "- Downsampling $B2 to $readsB reads.. "
    timeCMD "seqkit sample -s11 -2 -n $readsB -o $outdir/B.$1.R2.tmp.fq.gz $B2"
    
    echo ""
    echo "--------------------------------------------"
    CR1="$outdir/C.$1.${perc}_percent_contamination.R1.fq.gz"
    CR2="$outdir/C.$1.${perc}_percent_contamination.R2.fq.gz"
    echo "- Concatenate A+B samples to create C.$1 fq"
    echo "- creating $CR1 .."
    start=$(date +%s)
    cat $outdir/A.$1.R1.tmp.fq.gz $outdir/B.$1.R1.tmp.fq.gz > $CR1
    echo "- creating $CR2 .."
    cat $outdir/A.$1.R2.tmp.fq.gz $outdir/B.$1.R2.tmp.fq.gz > $CR2
    end=$(date +%s)
    t=$(($end - $start))
    printf 'Time taken: %d minutes, %d seconds\n' "$(( t/60 ))" "$(( t - 60*(t/60) ))"
    echo ""
    echo "--------------------------------------------"
    echo "- Final contaminated $1 sample written to "
    echo "R1: $CR1"
    echo "R2: $CR2"
    echo ""
    echo "--------------------------------------------"
    # echo "- Gzip C files"
    # echo "- R1"
    # if [ -f ${CR1}.gz ]; then
    #    rm -f ${CR1}.gz
    # fi
    # timeCMD "gzip $CR1"

    # echo "- R2"
    # if [ -f ${CR2}.gz ]; then
    #    rm -f ${CR2}.gz
    # fi
    # timeCMD "gzip $CR2"
    echo ""
    echo "--------------------------------------------"
    echo "-- Remove tmp files"
    rm -f $outdir/B.$1.R1.tmp.fq.gz
    rm -f $outdir/B.$1.R2.tmp.fq.gz
    rm -f $outdir/A.$1.R1.tmp.fq.gz
    rm -f $outdir/A.$1.R2.tmp.fq.gz

}

create_manifest() {
    echo "- Create Manifest yaml"
    mani="$outdir/manifest.yaml"
    echo "manifestVersion: v1alpha" > $mani
    echo "name: cont${perc}_${totreads}" >> $mani
    echo "bed: /samples/P23993/S31285117_grch38.bed" >> $mani
    echo "client: NOI" >> $mani
    echo "cancer: Lung Cancer" >> $mani
    echo "normalReads:" >> $mani
    echo "  - C.normal.${perc}_percent_contamination.R1.fq.gz " >> $mani
    echo "  - C.normal.${perc}_percent_contamination.R2.fq.gz " >> $mani
    echo "tumorReads: " >> $mani
    echo "  - C.tumor.${perc}_percent_contamination.R1.fq.gz " >> $mani
    echo "  - C.tumor.${perc}_percent_contamination.R2.fq.gz " >> $mani
    echo "rnaReads: " >> $mani
    echo "  - C.RNA.${perc}_percent_contamination.R1.fq.gz " >> $mani
    echo "  - C.RNA.${perc}_percent_contamination.R2.fq.gz " >> $mani
    
}

run() {

    set_params
    print_params

    # Read config file
    . $config 
    
    # Install seqkit via conda, if config says so
    if [[ $install == true ]]
    then
        echo "# - conda install seqkit .."
        conda install -c bioconda seqkit
    fi

    # Tumor DNA
    echo "==============================="
    echo "-> Processing Tumor DNA samples "
    echo "==============================="
    # P23993 tumor DNA
    A1="$tA1"
    A2="$tA2"
    # S2957 tumor DNA - contaminating
    B1="$tB1"
    B2="$tB2"
    process_sample "tumor"

    ## NORMAL DNA
    echo "==============================="
    echo "-> Processing Normal DNA samples "
    echo "==============================="
    # P23993 normal DNA
    A1="$nA1"
    A2="$nA2"
    # S2957 normal DNA
    B1="$nB1"
    B2="$nB2"
    process_sample "normal"

    echo " ========================= "
    echo " -> Processing tumor RNA sample "
    echo " ========================= "
    # P23993 tumor RNA
    A1="$rA1"
    A2="$rA2"
    # S2957 tumor RNA
    B1="$rB1"
    B2="$rB2"
    process_sample "RNA"

    create_manifest
}

config=$1
run config

