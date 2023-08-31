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
