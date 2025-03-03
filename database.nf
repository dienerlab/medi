#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.threads = 20
params.out = "${launchDir}/data"

workflow {
    def foodb = "https://foodb.ca/public/system/downloads/foodb_2020_4_7_csv.tar.gz"
    def genbank_summary = "https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt"
    def taxdump = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"

    download(foodb, genbank_summary) | get_taxids
    download_taxa_dbs(taxdump)
    get_lineage(get_taxids.out.combine(download_taxa_dbs.out))
        | match_taxids
        | download_sequences

    download_sequences.out.map{it[0]}.flatten().set{seqs}

    seqs | sketch
    ANI(sketch.out.collect())

    food_mappings(download.out, match_taxids.out, curate_taxids.out)
}



process download {
    cpus 1
    publishDir "${params.out}/dbs"

    input:
    val foodb
    val genbank_summary

    output:
    tuple path("foodb"), path("genbank_summary.tsv")

    script:
    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${params.foodb} -O foodb.tgz && \
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 ${params.genbank_summary} -O genbank_summary.tsv && \
    tar -xf foodb.tgz && \
    mv foodb_*_csv foodb
    """
}

process curate_taxids {
    cpus 1
    memory "4 GB"
    time "30 m"
    publishDir "${params.data_dir}"

    input:
    tuple path(foodb), path(gb_summary)
    path(curated)

    output:
    path("curated_foods.csv")

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    foods <- fread("foodb/Food.csv")
    foods[, "revised_taxonomy_id" := ncbi_taxonomy_id]
    setkey(foods, id)
    curated <- fread("curated")
    foods[curated[["id"]], revised_taxonomy_id := curated[["taxid"]]]
    fwrite(foods, "curated_foods.csv")
    """
}

process get_taxids {
    cpus 1
    memory "4 GB"
    time "1 h"

    input:
    tuple path(foodb), path(gb_summary)
    path(curated)

    output:
    tuple path("foodb"), path("taxids.tsv"), path("${gb_summary}")

    script:
    """
    #!/usr/bin/env Rscript

    library(data.table)

    dt <- fread("${gb_summary}", sep="\t")[
        grepl("ftp.ncbi.nlm.nih.gov", ftp_path, fixed = TRUE)
    ]
    genbank <- dt[!is.na(taxid), .(taxid = as.character(unique(taxid)))]
    genbank[, "source" := "genbank"]
    dt <- fread("${curated}")
    foodb <- dt[!is.na(ncbi_taxonomy_id), .(taxid = ncbi_taxonomy_id)]
    foodb[, "source" := "foodb"]
    fwrite(rbind(genbank, foodb), "taxids.tsv", col.names=F, sep="\\t")
    """
}

process download_taxa_dbs {
    cpus 1
    memory "500 MB"
    time "2 h"

    input:
    val(taxdump)

    output:
    path("taxdump")

    script:
    """
    wget --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 4 \
        ${taxdump} && \
        mkdir taxdump && tar -xf taxdump.tar.gz --directory taxdump
    """
}

process get_lineage {
    cpus 1
    memory "8 GB"
    time "1 h"

    input:
    tuple path(foodb), path(taxids), path(gb_summary), path(taxadb)

    output:
    tuple path("$foodb"), path("lineage.txt"), path("lineage_ids.txt"), path("${gb_summary}")

    script:
    """
    taxonkit lineage --data-dir $taxadb -i 1 $taxids > raw.txt && \
    taxonkit reformat --data-dir $taxadb -i 3 raw.txt > lineage.txt && \
    taxonkit reformat --data-dir $taxadb -t -i 3 raw.txt > lineage_ids.txt
    """
}

process match_taxids {
    cpus 1
    publishDir params.out
    memory "8 GB"
    time "1 h"

    input:
    tuple path(foodb), path(lineage), path(lineage_ids), path(gb_summary)

    output:
    path("matches.csv")

    script:
    """
    match.R $lineage_ids $gb_summary
    """
}

process download_sequences {
    cpus 8
    memory "64 GB"

    publishDir params.out
    time "48h"

    input:
    path(matches)

    output:
    tuple path("sequences/*.fna.gz"), path("manifest.csv")

    script:
    """
    download.R $matches $task.cpus sequences
    """
}


process food_mappings {
    cpus 1
    memory "64 GB"
    publishDir "${params.out}/dbs"
    time "30 m"

    input:
    tuple path(foodb), path(gb_summary)
    path(matches)
    path(curated)

    output:
    tuple path("food_matches.csv"), path("food_contents.csv.gz")

    script:
    """
    food_mapping.R ${foodb} ${matches} ${curated}
    """
}

process sketch {
    cpus 2
    memory "4 GB"
    publishDir "${params.out}/sketches"
    time "8 h"

    input:
    path(seq)

    output:
    path("*.sig")

    script:
    """
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000 ${seq}
    """
}

process ANI {
    cpus params.threads
    memory "64 GB"
    publishDir "${params.out}", mode: "copy", overwite: true
    time "8 h"

    input:
    path(sigs)

    output:
    path("mash_ani.csv")

    script:
    """
    sourmash compare -k 21 --ani -p ${task.cpus} --csv mash_ani.csv ${sigs}
    """
}
