#!/usr/bin/env nextflow
params.out_dir = "${launchDir}/data"
params.data_dir = "${launchDir}/data"
params.db = "${launchDir}/data/medi_db"
params.foods = "${params.db}/food_matches.csv"
params.food_contents = "${params.db}/food_contents.csv.gz"
params.single_end = false
params.trim_front = 5
params.min_length = 50
params.quality_threshold = 20
params.read_length = 150
params.threshold = 10
params.consistency = 0.95
params.entropy = 0.1
params.multiplicity = 4
params.confidence = 0.3
params.mapping = false
params.batchsize = 400
params.maxcpus = 24
params.dbmem = null
params.help = false

/* Helper functions */

// Helper to calculate the required RAM for the Kraken2 database
def estimate_db_size(hash) {
    def db_size = null

    // Calculate db memory requirement
    if (params.dbmem) {
        db_size = MemoryUnit.of("${params.dbmem} GB")
    } else {
        db_size = MemoryUnit.of(file(hash).size()) + 6.GB
        log.info("Based on the hash size I am reserving ${db_size.toGiga()}GB of memory for Kraken2.")
    }

    return db_size
}

def helpMessage() {
    log.info """
    ~~~ MEDI - Metagenomic Estimation of Dietary Intake ~~~

    Usage:
        nextflow run quant.nf -resume -with-conda /path/to/conda/env --db /path/to/medi_db

    Required arguments:
        --db <db_dir>               Directory containing the MEDI database.

    Options:
      General:
        --single_end                Set if the data is single-end. Default: false.
        --trim_front <int>          Number of bases to trim from the 5' end of reads. Default: 5.
        --min_length <int>          Minimum length of reads to keep. Default: 50.
        --quality_threshold <int>   Minimum quality score to keep a base. Default: 20.
        --batchsize <int>           Number of samples to process in a batch. Default: 400.
        --maxcpus <int>             Maximum number of CPUs to use. Default: 24.
        --dbmem <int>               Memory to reserve for the Kraken2 database. Default: calculate automatically from DB.
      Bracken:
        --read_length <int>         Read length for Bracken. Default: 150.
        --threshold <int>           Bracken read threshold. Default: 10.
        --mapping                   Generate mapping summaries. Default: false.
      Architeuthis (decoy filtering):
        --consistency <float>       Minimum consistency score for a classification to be kept. Default: 0.95.
        --entropy <float>           Maximum entropy score for a classification to be kept. Default: 0.1.
        --multiplicity <int>        Maximum multiplicity of a classification to be kept. Default: 4.
""".stripIndent()
}

/********************/

/* Workflow definition starts here */

workflow {
        // Show help message
    if (params.help) {
        helpMessage()
        exit 0
    }

    Channel
        .fromList(["D", "G", "S"])
        .set{levels}

    // find files
    if (params.single_end) {
        Channel
            .fromPath("${params.data_dir}/raw/*.fastq.gz")
            .map{row -> tuple(row.baseName.split("\\.fastq")[0], tuple(row))}
            .set{raw}
    } else {
        Channel
            .fromFilePairs([
                "${params.data_dir}/raw/*_R{1,2}_001.fastq.gz",
                "${params.data_dir}/raw/*_{1,2}.fastq.gz",
                "${params.data_dir}/raw/*_R{1,2}.fastq.gz"
            ])
            .ifEmpty { error "Cannot find any read files in ${params.data_dir}/raw!" }
            .set{raw}
    }

    // quality filtering
    preprocess(raw)

    // quantify taxa abundances
    batched = preprocess.out    // buffer the samples into batches
        .collate(params.batchsize)
        .map{it -> it.collect{a -> a[1]}.flatten().sort()}

    // run Kraken2 per batch
    kraken(batched)

    // filter Kraken2 results and generate reports
    kraken.out
        .flatMap{it[0]}
        .map{tuple it.baseName.split(".k2")[0], it}
        | architeuthis_filter | kraken_report
    count_taxa(kraken_report.out.combine(levels))
    count_taxa.out.map{s -> tuple(s[1], s[2])}
        .groupTuple()
        .set{merge_groups}
    merge_taxonomy(merge_groups)

    if (params.mapping) {
        // Get individual mappings
        summarize_mappings(architeuthis_filter.out)
        summarize_mappings.out.collect() | merge_mappings
    }

    // Add taxon lineages
    add_lineage(merge_taxonomy.out)

    // Quantify foods
    add_lineage.out.collect() | quantify

    // quality overview
    multiqc(merge_taxonomy.out.map{it[1]}.collect())
}

/* Process definitions */

process preprocess {
    cpus 4
    memory "6 GB"
    publishDir "${params.out_dir}/preprocessed"
    time "1h"

    input:
    tuple val(id), path(reads)

    output:
    tuple val(id),
        path("${id}_filtered_R*.fastq.gz"),
        path("${id}_fastp.json"),
        path("${id}.html")

    script:
    if (params.single_end)
        """
        fastp -i ${reads[0]} -o ${id}_filtered_R1.fastq.gz \
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus}
        """

    else
        """
        fastp -i ${reads[0]} -I ${reads[1]} \
            -o ${id}_filtered_R1.fastq.gz -O ${id}_filtered_R2.fastq.gz\
            --json ${id}_fastp.json --html ${id}.html \
            --trim_front1 ${params.trim_front} -l ${params.min_length} \
            -3 -M ${params.quality_threshold} -r -w ${task.cpus}
        """
}

process kraken {
    cpus params.maxcpus
    memory { estimate_db_size("${params.db}/hash.k2d") }
    time { 2.h + reads.size() * 0.5.h }
    scratch false
    publishDir "${params.data_dir}/kraken2"

    input:
    path(reads)

    output:
    tuple path("*.k2"), path("*.tsv")

    script:
    """
    #!/usr/bin/env python

    import sys
    import os
    from subprocess import run

    base_args = [
        "kraken2", "--db", "${params.db}",
        "--confidence", "${params.confidence}",
        "--threads", "${task.cpus}", "--gzip-compressed"
    ]

    reads = "${reads}".split()

    se = ${params.single_end ? "True" : "False"}
    if se:
        fwd = reads
    else:
        fwd = reads[0::2]
        rev = reads[1::2]
        base_args += ["--paired"]
        assert len(fwd) == len(rev)
        assert (
            [f.split("_filtered_R")[0] for f in fwd] ==
            [r.split("_filtered_R")[0] for r in rev]
        )

    for i, read in enumerate(fwd):
        idx = read.split("_filtered_R")[0]
        args = base_args + [
            "--output", f"{idx}.k2",
            "--report", f"{idx}.tsv",
            "--memory-mapping", read
        ]
        if not se:
            args.append(rev[i])
        res = run(args)
        if res.returncode != 0:
            if os.path.exists(f"{idx}.k2"):
                os.remove(f"{idx}.k2")
            sys.exit(res.returncode)
    """
}

process architeuthis_filter {
    cpus 1
    publishDir "${params.out_dir}/kraken2", overwrite: true
    time 1.h
    memory "2 GB"

    input:
    tuple val(id), path(k2)

    output:
    tuple val(id), path("${id}_filtered.k2")

    script:
    """
    architeuthis mapping filter ${k2} \
        --data-dir ${params.db}/taxonomy \
        --min-consistency ${params.consistency} --max-entropy ${params.entropy} \
        --max-multiplicity ${params.multiplicity} \
        --out ${id}_filtered.k2
    """
}

process kraken_report {
    cpus 1
    memory "200 MB"
    publishDir "${params.out_dir}/kraken2", overwrite: true
    time 30.m

    input:
    tuple val(id), path(k2)

    output:
    tuple val(id), path("*.tsv")

    script:
    """
    kraken2-report ${params.db}/taxo.k2d ${k2} ${id}.tsv
    """
}

process summarize_mappings {
    cpus 1
    publishDir "${params.out_dir}/architeuthis"
    time 1.h

    input:
    tuple val(id), path(k2), path(report)

    output:
    path("${id}_mapping.csv")

    script:
    """
    architeuthis mapping summary ${k2} --data-dir ${params.db}/taxonomy --out ${id}_mapping.csv
    """
}

process merge_mappings {
    cpus 1
    publishDir "${params.out_dir}", mode: "copy", overwrite: true
    time 1.h

    input:
    path(mappings)

    output:
    path("mappings.csv")

    script:
    """
    architeuthis merge ${mappings} --out mappings.csv
    """
}

process count_taxa {
    cpus 4
    memory "640 MB"
    publishDir "${params.out_dir}/bracken", overwrite: true
    time 1.h

    input:
    tuple val(id), path(report), val(lev)

    output:
    tuple val(id), val(lev), path("${lev}/${lev}_${id}.b2")

    script:
    """
    mkdir ${lev} && \
        fixk2report.R ${report} ${lev}/${report} && \
        bracken -d ${params.db} -i ${lev}/${report} \
        -l ${lev} -o ${lev}/${lev}_${id}.b2 -r ${params.read_length} \
        -t ${params.threshold} -w ${lev}/${id}_bracken.tsv
    """
}

process quantify {
    cpus 1
    memory "2 GB"
    publishDir "${params.out_dir}", mode: "copy", overwrite: true
    time 2.h

    input:
    path(files)

    output:
    tuple path("food_abundance.csv"), path("food_content.csv")

    script:
    """
    quantify.R ${params.foods} ${params.food_contents} ${files}
    """
}

process merge_taxonomy {
    cpus 1
    memory "1 GB"
    time 2.h

    input:
    tuple val(lev), path(reports)

    output:
    tuple val(lev), path("${lev}_merged.csv")

    script:
    """
    architeuthis merge ${reports} --out ${lev}_merged.csv
    """
}

process add_lineage {
    cpus 1
    memory "4 GB"
    publishDir "${params.out_dir}", mode: "copy", overwrite: true
    time 2.h

    input:
    tuple val(lev), path(merged)

    output:
    path("${lev}_counts.csv")

    script:
    """
    architeuthis lineage ${merged} --data-dir ${params.db}/taxonomy --out ${lev}_counts.csv
    """
}


process multiqc {
    cpus 1
    memory "2 GB"
    publishDir "${params.out_dir}", mode: "copy", overwrite: true
    time 2.h

    input:
    path(report)

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc ${params.out_dir}/preprocessed ${params.out_dir}/kraken2
    """
}
