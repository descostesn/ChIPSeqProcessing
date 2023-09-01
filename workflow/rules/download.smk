################################################################################
# Download sra files
#
# from the csv of TEbench. The file can be found in config/data.csv
# WARNING: This code is made to work with TEbench
################################################################################


rule download_sra_single:
  output:
    singleSRA = "../results/data/sra/single/{singleEndName}"
  params:
    linksingle = lambda wildcards: samples_single_df.loc[wildcards.singleEndName,
         "link"],
    srasingle = lambda wildcards: samples_single_df.loc[wildcards.singleEndName, 
        "sraname"]
  threads: 1
  benchmark: "benchmark/download_sra_single/{singleEndName}.tsv" 
  shell:
    """
    echo "Downloading sra {params.linksingle}"
    wget --directory-prefix=../results/data/sra/single/ {params.linksingle}
    sleep 10s
    mv ../results/data/sra/single/{params.srasingle} {output.singleSRA}  
    """

rule download_sra_paired:
  output:
    pairedSRA = "../results/data/sra/paired/{pairedEndName}"
  params:
    linkpaired = lambda wildcards: samples_paired_df.loc[wildcards.pairedEndName, 
        "link"],
    srapaired = lambda wildcards: samples_paired_df.loc[wildcards.pairedEndName, 
        "sraname"]
  threads: 1
  benchmark: "benchmark/download_sra_paired/{pairedEndName}.tsv"
  shell:
    """
    echo "Downloading sra {params.linkpaired}"
    wget --directory-prefix=../results/data/sra/paired/ {params.linkpaired}
    sleep 10s
    mv ../results/data/sra/paired/{params.srapaired} {output.pairedSRA}
    """

################################################################################
# Download genome related files
#
# The links to the files are in the configuration file in ../config/config.yaml
# WARNING: This code is made to work with TEbench
################################################################################

rule download_genome_fasta:
  output:
    "../results/data/fasta/{genome}/{prefix}.dna.chromosome.fa",
  threads: 1
  benchmark: "benchmark/download_genome_fasta/{genome}/{prefix}.tsv"
  params:
    genomeLinks = lambda wildcards: GENOMELINKS
  shell:
    """
    echo "Downloading {wildcards.genome} fasta"
    for link in {params.genomeLinks}
    do
        wget $link
    done
    
    echo "Merging files"
    cat {wildcards.prefix}.dna.chromosome.*.fa.gz > final-{wildcards.prefix}
    rm {wildcards.prefix}.dna.chromosome.*.fa.gz
    mv final-{wildcards.prefix} {wildcards.prefix}.dna.chromosome.fa.gz
    gunzip {wildcards.prefix}.dna.chromosome.fa.gz

    echo "Organizing files to destination folders"
    mkdir -p ../results/data/fasta/{wildcards.genome}
    mv {wildcards.prefix}.dna.chromosome.fa ../results/data/fasta/{wildcards.genome}
    """

rule download_genome_annotations:
  output:
    "../results/data/gtf/{genome}/{prefix}.{version}.chr.gtf",
    "../results/data/gff/{genome}/{prefix}.{version}.chr.gff3"
  threads: 1
  benchmark: "benchmark/download_genome_annotations/{genome}/{prefix}.{version}.tsv"
  params:
    pathGTF = lambda wildcards: config["genome"]["gtfLink"],
    pathGFF = lambda wildcards: config["genome"]["gffLink"]
  shell:
    """
    echo "Downloading GTF"
    wget {params.pathGTF}
    
    echo "Filtering MT from GTF"
    gunzip {wildcards.prefix}.{wildcards.version}.chr.gtf.gz
    grep -v MT {wildcards.prefix}.{wildcards.version}.chr.gtf > tmp-{wildcards.prefix}
    rm {wildcards.prefix}.{wildcards.version}.chr.gtf
    mv tmp-{wildcards.prefix} {wildcards.prefix}.{wildcards.version}.chr.gtf
    
    echo "Downloading GFF"
    wget {params.pathGFF}
    
    echo "Filtering MT from GFF"
    gunzip {wildcards.prefix}.{wildcards.version}.chr.gff3.gz
    grep -v MT {wildcards.prefix}.{wildcards.version}.chr.gff3 > tmp
    rm {wildcards.prefix}.{wildcards.version}.chr.gff3
    mv tmp {wildcards.prefix}.{wildcards.version}.chr.gff3
    
    echo "Organizing files to destination folders"
    mkdir -p ../results/data/gtf/{wildcards.genome}
    mkdir -p ../results/data/gff/{wildcards.genome}
    mv {wildcards.prefix}.{wildcards.version}.chr.gtf ../results/data/gtf/{wildcards.genome}
    mv {wildcards.prefix}.{wildcards.version}.chr.gff3 ../results/data/gff/{wildcards.genome}
    """
