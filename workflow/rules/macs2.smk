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

def retrieve_elongation(wildcards):
    FilePath = "../results/qc/bowtie2_saturation_percentmultireads/wildcards.genome/single/wildcards.samplenames.txt"
    SampleName = wildcards.singleexpsync
    with open(FilePath, 'r') as f:
      mat = [[element.strip() for element in line.split('\t')] for line in f]
    
    # Retrieve alignment threshold from samplenames
    bowtie2thres = SampleName.split("_trimmed_")[1].split("_sorted")[0]



rule macs2_narrow_single:
  input:
    chipexp = "../results/bam/single/bowtie2_results/{genome}/{singleexpsync}.bam",
    controlexp = "../results/bam/single/bowtie2_results/{genome}/{singleinputsync}.bam",
    macs2info = "../results/qc/info_macs2/single/{samplenames}.txt",
    elongation = "../results/qc/bowtie2_saturation_percentmultireads/{genome}/single/{samplenames}.txt"
    !!!!!!!!!!!!!!!!!!!!
!!! ADD THE RETRIEVAL OF THE ELONGATION SIZE FOR MACS. CHECK IF singleexpsync IS EQUIVALENT TO singlebestmulti FOR sorted_nodups
!!!!!!!!!!!!!!!!!!!!

  output:
    "../results/peak_detection/single/macs2/{genome}/{qvalthres}/{modeltype}/{singleexpsync}_control_{singleinputsync}_info_{samplenames}_peaks.narrowPeak"
  threads: 1
  conda: "../envs/macs2.yaml"
  benchmark: "benchmark/macs2_narrow_single/{genome}/single/{samplenames}.tsv"
  params:
    genome_size = config["genome"]["size"]
  shell:
    """
    ## Retrieve in the file obtained with the above rule 'retrieve_singleinfo' the number and length of reads
    nbseq=`grep "Total Sequences" {input.macs2info} | sed -e "s/Total\sSequences\s//""` 
    sqlength=`grep "Sequence length" {input.macs2info} | sed -e "s/Sequence\slength\s//"`

    ## Compute the threshold for duplicates which increases of 1 every 7 million reads
    thres=`echo "scale=0 ; $nbseq / 7000000" | bc`

    ## Define the output folder and the threshold used for bowtie2 alignment
    outfold=`dirname {output}`
    thresalign=`ls *.txt | awk -F'_trimmed_' '{print $2}' | awk -F'_sorted' '{print $1}'`

    echo "---- Creating nomodel wihtout broad\n"
    macs2 callpeak -t {input.chipexp} -c {input.controlexp} -n {wildcards.singleexpsync} --outdir $outfold -f 'BAM' -g {params.genome_size} -s $sqlength -q {wildcards.qvalthres} --nomodel --extsize 150 --keep-dup $thres
    """




rule macs2_broad_single:
  input:
    chipexp = "../results/bam/single/bowtie2_results/{genome}/{singleexpsync}.bam",
    controlexp = "../results/bam/single/bowtie2_results/{genome}/{singleinputsync}.bam",
    macs2info = "../results/qc/info_macs2/single/{samplenames}.txt",
    !!!!!!!!!!!!!!!!!!!!
!!! ADD THE RETRIEVAL OF THE ELONGATION SIZE FOR MACS. CHECK IF singleexpsync IS EQUIVALENT TO singlebestmulti FOR sorted_nodups
!!!!!!!!!!!!!!!!!!!!

  output:
    "../results/peak_detection/single/macs2/{genome}/{qvalthres}/no_model_broad/{singleexpsync}_control_{singleinputsync}_info_{samplenames}_peaks.broadPeak"
  threads: 1
  conda: "../envs/macs2.yaml"
  !! benchmark:
  params:
    genome_size = config["genome"]["size"]
  shell:
    """
    nbseq=`grep "Total Sequences" {input.macs2info} | sed -e "s/Total\sSequences\s//""` 
    sqlength=`grep "Sequence length" {input.macs2info} | sed -e "s/Sequence\slength\s//"`
    thres=`echo "scale=0 ; $nbseq / 7000000" | bc`
    outfold=`dirname {output}`
    macs2 callpeak -t {input.chipexp} -c {input.controlexp} -n {wildcards.singleexpsync} --outdir $outfold -f 'BAM' -g {params.genome_size} -s $sqlength --nomodel --extsize 150 --keep-dup $thres --broad --broad-cutoff {wildcards.qvalthres}
    """
