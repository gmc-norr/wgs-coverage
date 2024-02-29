params.bam = null
params.bamlist = null
params.genome = null
params.cytobands = true
params.results = "results"
params.regions = null

log.info """\
    W G S  C O V E R A G E
    ======================
        bam: ${params.bam}
    bamlist: ${params.bamlist}
     genome: ${params.genome}
  cytobands: ${params.cytobands}
    results: ${params.results}
    regions: ${params.regions}
    """

process mosdepth {
    container 'quay.io/biocontainers/mosdepth:0.3.6--hd299d5a_0'
    publishDir "${params.results}/${sample}"

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), path('*.per-base.d4'), emit: per_base_d4
    tuple val(sample), path('*.mosdepth.global.dist.txt'), emit: global_dist
    tuple val(sample), path('*.mosdepth.summary.txt'), emit: summary

    script:
    """
    mosdepth -t ${task.cpus} --d4 ${sample} $bam
    """
}

process samtools_index {
    container 'quay.io/biocontainers/samtools:1.19.1--h50ea8bc_0'

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path("*.bai"), emit: bai

    script:
    """
    samtools index $bam
    """
}

process plot_coverage {
    conda "${projectDir}/environments/plot_coverage.yaml"
    publishDir "${params.results}/${sample}/plots"

    input:
    tuple val(sample), path(coverage)
    path cytobands
    path regions

    output:
    tuple val(sample), path('*.png'), emit: plots

    script:
    cytoband_arg = ""
    if (params.cytobands) {
        cytoband_arg = "--cytobands $cytobands"
    }

    region_arg = ""
    if (regions.name != "NO_FILE") {
        region_arg = "--regions $regions"
    }

    """
    plot_coverage.py $region_arg $cytoband_arg -o ${sample} $coverage
    """
}

process plot_gene_coverage_distributions {
    conda "${projectDir}/environments/plot_coverage.yaml"
    publishDir "${params.results}/summary"

    input:
    path coverage
    path regions

    output:
    path "*.png", emit: plots

    """
    plot_gene_coverage_distributions.py $regions $coverage
    """
}

process split_bed {
    container 'quay.io/biocontainers/csvtk:0.29.0--h9ee0642_0'

    input:
    path bed

    output:
    path "*.${extension}", emit: split_bed

    script:
    extension = bed.getExtension()
    """
    csvtk split --no-header-row --tabs -f 4 $bed
    """
}

process guess_genome {
    container "quay.io/biocontainers/pysam:0.22.0--py39hcada746_0"

    input:
    tuple val(sample), path(bam), path(bai)

    output:
    tuple val(sample), stdout

    script:
    """
    guess_genome.py $bam
    """
}

def is_newer(a, b) {
    return file(a).lastModified() > file(b).lastModified()
}

workflow {
    if (!params.bam && !params.bamlist) {
        error("No bam file(s) provided")
    }

    if (params.bam && params.bamlist) {
        error("Cannot provide both a bam file and a bam list")
    }

    if (params.bamlist) {
        bams = file(params.bamlist).readLines().collect { file(it, checkIfExists: true) }
    } else {
        bams = [file(params.bam, checkIfExists: true)]
    }

    bam_ch = Channel
        .fromList(bams)
        .map { bam ->
            sample = bam.getBaseName().split("_").first()
            [sample, bam]
        }

    bam_wo_index = bam_ch.filter { it ->
        bai = file("${it[1]}.bai")
        !bai.exists() || !is_newer(bai, it[1])
    }

    bai_ch = samtools_index(bam_wo_index)

    bam_ch = bam_ch
        .join(bai_ch, remainder: true)
        .map({ [it[0], it[1], it[2] ? it[2] : "${it[1]}.bai"] })

    if (!params.genome) {
        genome_ch = guess_genome(bam_ch)
        genome_ch.view { log.info "guessing genome build is ${it[1]} for ${it[0]}" }
        cytoband_ch = genome_ch.map { file("${projectDir}/data/cytoBand.${it[1]}.txt") }
    } else {
        genome_ch = Channel.value(params.genome)
        cytoband_ch = Channel.value(file("${projectDir}/data/cytoBand.${params.genome}.txt"))
    }

    regions_ch = Channel.value(file("${projectDir}/assets/NO_FILE"))
    if (params.regions) {
        regions_ch = Channel.value(file(params.regions))
    }

    coverage = mosdepth(bam_ch)
    plot_coverage(coverage.per_base_d4, cytoband_ch, regions_ch)
    plot_gene_coverage_distributions(coverage.per_base_d4.collect { it[1] }, regions_ch)
}
