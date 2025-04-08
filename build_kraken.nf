#!/usr/bin/env nextflow
include { food_mappings } from './database.nf'

nextflow.enable.dsl = 2

params.additionalDbs = ["bacteria", "archaea", "human", "viral", "plasmid", "UniVec_Core"]
params.maxDbSize = "750 GB"
params.confidence = 0.3
params.threads = 20
params.rebuild = false
params.out = "${launchDir}/data"
params.db = "${params.out}/medi_db"
params.additionalDecoys ="${params.out}/decoys"

/* Helper functions */

// Helper to calculate the required RAM for the Kraken2 database
def estimate_db_size(hash, extra) {
    def db_size = null

    // Calculate db memory requirement
    db_size = MemoryUnit.of(file(hash).size()) + extra
    log.info("Based on the hash size I am reserving ${db_size.toGiga()}GB of memory for Kraken2.")

    return db_size
}

workflow {
    if (!params.rebuild) {
        Channel.fromPath("${params.out}/sequences/*.fna.gz").set{food_sequences}

        Channel.fromPath("${params.additionalDecoys}/*.fna.gz")
        .set{decoy_sequences}

        sequences = decoy_sequences.concat(food_sequences)

        taxonomy = setup_kraken_db()
        add_existing(taxonomy, params.additionalDbs)
        add_sequences(sequences, taxonomy)
        lib = assemble_library(
            taxonomy,
            add_existing.out.collect(),
            add_sequences.out.collect()
        )
    } else {
        taxonomy = Channel.of(file("${params.db}/taxonomy"))
        lib = Channel.of(file("${params.db}/library"))
    }

    build_kraken_db(taxonomy, lib)
    self_classify(build_kraken_db.out)
    build_bracken(self_classify.out) | add_info
}


process setup_kraken_db {
    cpus 1
    memory "32 GB"
    time "12 h"

    output:
    path("medi_db/taxonomy")

    script:
    """
    kraken2-build --download-taxonomy --db medi_db
    """
}

process add_sequences {
    cpus 4
    memory "32 GB"
    time "48 h"

    input:
    path(fasta)
    path(taxonomy)

    output:
    path("medi_db/library/added/*.*")

    script:
    """
    mkdir medi_db && mv ${taxonomy} medi_db
    gunzip -c $fasta > ${fasta.baseName} && \
    kraken2-build --add-to-library ${fasta.baseName} --db medi_db --threads ${task.cpus} && \
    rm ${fasta.baseName}
    """
}

process add_existing {
    cpus 4
    memory "32 GB"
    time "48 h"

    input:
    path(taxonomy)
    each group

    output:
    path("medi_db/library/$group")

    script:
    if (group == "human")
        """
        mkdir medi_db && mv ${taxonomy} medi_db
        kraken2-build --download-library $group --db medi_db --no-mask --threads ${task.cpus}
        """
    else
        """
        mkdir medi_db && mv taxonomy medi_db
        kraken2-build --download-library $group --db medi_db --threads ${task.cpus}
        """
}

process assemble_library {
    cpus 1
    memory "8 GB"
    time "12h"
    publishDir params.db

    input:
    path(taxonomy)
    path(existing)
    path(sequences)

    output:
    path("medi_db/library")

    script:
    """
    mkdir medi_db && mkdir medi_db/library && mkdir medi_db/library/added
    mv ${taxonomy} medi_db
    mv ${existing} medi_db/library
    mv ${sequences} medi_db/library/added
    """
}

process build_kraken_db {
    cpus params.threads
    memory { MemoryUnit.of(params.maxDbSize) + 16.GB }
    cpus params.max_threads
    time "72 h"
    publishDir params.db

    input:
    path(taxonomy)
    path(library)

    output:
    tuple path("medi_db/taxonomy"), path("medi_db/*.k2d")

    script:
    """
    mkdir medi_db && mv ${taxonomy} ${library} medi_db
    kraken2-build --build --db medi_db \
        --threads ${task.cpus} \
        --max-db-size ${task.memory.toGiga()}
    """
}


process self_classify {
    cpus 1
    memory { estimate_db_size("${db}/hash.k2d", 64.GB) }
    time "48 h"

    input:
    path(db)

    output:
    path(db)

    script:
    """
    kraken2 --db ${db} --threads ${task.cpus} \
        --confidence ${params.confidence} \
        --threads ${task.cpus} \
        --memory-mapping ${db}/library/*/*.f*a > ${db}/database.kraken
    """
}

process build_bracken {
    cpus 20
    memory "64 GB"
    publishDir params.out
    time "12 h"

    input:
    path(db)

    output:
    path("$db")

    script:
    """
    bracken-build -d $db -t ${task.cpus} -k 35 -l 100 && \
    bracken-build -d $db -t ${task.cpus} -k 35 -l 150
    """
}

process library {
    cpus 1
    memory "4 GB"
    time "1 h"

    input:
    path(db)

    output:
    path("$db/library/*/*.f*a")

    script:
    """
    ls ${db}/library/*/*.f*a | wc -l
    """
}

process add_info {
    cpus 1
    memory "1 GB"
    publishDir params.out

    input:
    path(db)

    output:
    path("$db")

    script:
    """
    cp ${params.downloads}/dbs/{food_matches.csv,food_contents.csv.gz} ${db}
    cp ${params.downloads}/manifest.csv ${db}
    """
}