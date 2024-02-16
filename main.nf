params.bam = null
params.genome = null
params.cytobands = true
params.results = "results"

log.info """\
    W G S  C O V E R A G E
    ======================
        bam: ${params.bam}
     genome: ${params.genome}
  cytobands: ${params.cytobands}
    results: ${params.results}
    """

process mosdepth {
    container 'quay.io/biocontainers/mosdepth:0.3.6--hd299d5a_0'
    publishDir "${params.results}/${sample}"

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    path '*.per-base.d4', emit: 'per_base_d4'
    path '*.mosdepth.global.dist.txt', emit: 'global_dist'
    path '*.mosdepth.summary.txt', emit: 'summary'

    script:
    """
    mosdepth -t ${task.cpus} --d4 ${sample} $bam
    """
}

process samtools_index {
    container 'quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0'
    
    input:
    path bam

    output:
    path "*.bai", emit: "bai"

    script:
    """
    samtools index $bam
    """
}

process plot_coverage {
    conda "environments/plot_coverage.yaml"
    publishDir "${params.results}/${sample}"

    input:
    val sample
    path coverage
    val genome
    path cytobands

    output:
    path '*.total_coverage.png'

    script:
    cytoband_arg = ""
    if (params.cytobands) {
        cytoband_arg = "--cytobands $cytobands"
    }

    if (genome) {
        genome = "--reference-genome $genome"
    }

    """
    plot_coverage.py $cytoband_arg $genome -o ${sample}.total_coverage.png $coverage
    """
}

process guess_genome {
    container "quay.io/biocontainers/pysam:0.22.0--py39hcada746_0"

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    stdout

    script:
    """
    guess_genome.py $bam
    """
}

def is_newer(a, b) {
    return file(a).lastModified() > file(b).lastModified()
}

workflow {
    if (params.bam == null) {
        error("No bam file provided")
    }

    bam = file(params.bam, checkIfExists: true)
    bai = file("${params.bam}.bai")
    sample = bam.getBaseName().split("_").first()

    if (!bai.exists() || !is_newer(bai, bam)) {
        bai_ch = samtools_index(Channel.fromPath(bam))
    } else {
        bai_ch = Channel.fromPath(bai)
    }

    bam_ch = Channel.of([sample, bam]).combine(bai_ch)

    if (!params.genome) {
        genome_ch = guess_genome(bam_ch)
        genome_ch.view { log.info "guessing genome build is $it" }
    } else {
        genome_ch = Channel.value(params.genome)
    }

    cytoband_ch = genome_ch.map { file("data/cytoBand.${it}.txt") }

    coverage = mosdepth(bam_ch)
    plot_coverage(sample, coverage.per_base_d4, genome_ch, cytoband_ch)
}
