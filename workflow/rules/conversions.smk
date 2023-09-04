################################################################################
# Convert sra files to fastq after downloading them followed by compression.
# WARNING: This code is made to work with TEbench
################################################################################

rule fasterq_dump_single:
  input:
    rules.download_sra_single.output.singleSRA
  output:
    singleFastq = "../results/data/single/{singleEndName}.fastq"
  threads: 20
  conda: "../../envs/parallelfastqdump.yaml"
  benchmark: "benchmark/fasterq_dump_single/{singleEndName}.tsv"
  shell:
    """
    echo "Converting to fastq"
    cd ../results/data/sra/single
    mkdir fastqdump-{wildcards.singleEndName}
    parallel-fastq-dump -s ./{wildcards.singleEndName} -t {threads} -O ./ --tmpdir ./fastqdump-{wildcards.singleEndName}
    sleep 10s
    mv ./{wildcards.singleEndName}.fastq ../../single
    rm -r fastqdump-{wildcards.singleEndName}
    cd ../../..
    sleep 10s
    """

rule fasterq_dump_paired:
  input:
    rules.download_sra_paired.output.pairedSRA
  output:
    pairedFastq_1 = "../results/data/paired/{pairedEndName}_1.fastq",
    pairedFastq_2 = "../results/data/paired/{pairedEndName}_2.fastq"
  threads: 20
  conda: "../../envs/parallelfastqdump.yaml"
  benchmark: "benchmark/fasterq_dump_paired/{pairedEndName}.tsv"
  shell:
    """
    echo "Converting to fastq"
    cd ../results/data/sra/paired 
    mkdir fastqdump-{wildcards.pairedEndName}
    parallel-fastq-dump -s ./{wildcards.pairedEndName} -t {threads} -O ./ --tmpdir ./fastqdump-{wildcards.pairedEndName} --split-files
    sleep 10s
    mv ./{wildcards.pairedEndName}_1.fastq ../../paired
    mv ./{wildcards.pairedEndName}_2.fastq ../../paired
    rm -r fastqdump-{wildcards.pairedEndName}
    cd ../../..
    sleep 10s
    """

rule gzip_fastq_single:
  input:
    rules.fasterq_dump_single.output.singleFastq
  output:
    singleFastq = "../results/data/single/{singleEndName}.fastq.gz"
  threads: 1
  benchmark: "benchmark/gzip_fastq_single/{singleEndName}.tsv"
  shell:
    """
    echo "Compressing files"
    gzip -c {input} > {output.singleFastq}
    """

rule gzip_fastq_paired1:
  input:
    pairedFastq_1 = "../results/data/paired/{pairedEndName}_1.fastq"
  output:
    pairedFastq_1 = "../results/data/paired/{pairedEndName}_1.fastq.gz"
  threads: 1
  benchmark: "benchmark/gzip_fastq_paired1/{pairedEndName}.tsv"
  shell:
    """
    echo "Compressing files"
    gzip -c {input.pairedFastq_1} > {output.pairedFastq_1}
    """
    
rule gzip_fastq_paired2:
  input:
    pairedFastq_2 = "../results/data/paired/{pairedEndName}_2.fastq"
  output:
    pairedFastq_2 = "../results/data/paired/{pairedEndName}_2.fastq.gz"
  threads: 1
  benchmark: "benchmark/gzip_fastq_paired2/{pairedEndName}.tsv"
  shell:
    """
    echo "Compressing files" 
    gzip -c {input.pairedFastq_2} > {output.pairedFastq_2}
    """
    