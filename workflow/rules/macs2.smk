################################################################################
# Performing peak detection with macs2
# WARNING: This code is made to work with another workflow
################################################################################

# This rule retrieve read size and number of reads
rule retrieve_singleinfo:
  input:
    "../results/qc/fastqc/raw_fastq/single/{singleEndName}_fastqc.zip"
  output:
    "../results/qc/info_macs2/single/{singleEndName}.txt"
  shell:
    """
    unzip {input}
    foldname=`test=`ls {input} | sed -e "p;s/\.zip//"`
    cd $foldname
    grep -e "Total Sequences" -e "Sequence length" fastqc_data.txt > {ouput}
    cd ..
    rm -r $foldname
    """

rule macs2_narrow_single:
  input:
    chipexp = "../results/bam/single/bowtie2_results/{genome}/{singleexpsync}.bam",
    controlexp = "../results/bam/single/bowtie2_results/{genome}/{singleinputsync}.bam",
    macs2info = "../results/qc/info_macs2/single/{samplenames}.txt"
  output:
    "../results/peak_detection/single/macs2/{genome}/{qvalthres}/{modeltype}/{singleexpsync}_control_{singleinputsync}_info_{samplenames}_peaks.narrowPeak"
  params:
    genome_size = config["genome"]["size"]
  shell:
    """
    nbseq=`grep "Total Sequences" {input.macs2info} | sed -e "s/Total\sSequences\s//""` 
    sqlength=`grep "Sequence length" {input.macs2info} | sed -e "s/Sequence\slength\s//"`
    thres=`echo "scale=0 ; $nbseq / 7000000" | bc`
    outfold=`dirname {output}`
    echo "---- Creating nomodel wihtout broad\n"
    macs2 callpeak -t {input.chipexp} -c {input.controlexp} -n {wildcards.singleexpsync} --outdir $outfold -f 'BAM' -g {params.genome_size} -s $sqlength -q {wildcards.qvalthres} --nomodel --extsize 150 --keep-dup $thres
    """

# rule macs2_broad_single:
#   input:
#     rules.remove_duplicates_bowtie_single.output.bamFile
#   output:
#     "../results/peak_detection/single/macs2/{genome}/{qvalthres}/no_model_broad/{singlebestmulti}_peaks.broadPeak"

