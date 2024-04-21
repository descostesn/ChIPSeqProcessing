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
    bamfile = "../results/bam/single/bowtie2_results/{genome}/{singlebestmulti}_sorted_nodups.bam",
    macs2info = "../results/qc/info_macs2/single/{singleEndName}.txt"
  output:
    "../results/peak_detection/single/macs2/{genome}/{qvalthres}/{modeltype}/{singlebestmulti}_peaks.narrowPeak"
  shell:
    """
    nbseq=`grep "Total Sequences" 
    """
!!
echo "---- Creating nomodel wihtout broad\n"
macs2 callpeak -t {input} -c $input_file_vector -n $experiment_name --outdir $output_folder_nomodel -f $format -g $genome_size -s $tag_size -q $qvalue --nomodel --extsize $elongation_size --keep-dup $artefact_threshold";
!!




rule macs2_broad_single:
  input:
    rules.remove_duplicates_bowtie_single.output.bamFile
  output:
    "../results/peak_detection/single/macs2/{genome}/{qvalthres}/no_model_broad/{singlebestmulti}_peaks.broadPeak"

