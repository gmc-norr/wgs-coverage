# WGS clinical coverage

A Nextflow pipeline that calculates and visualises read coverage in WGS data.

## Usage

```
nextflow run gmc-norr/wgs-coverage \
    --bam <input.bam> | --bamlist <input.txt> \
    [--genome hg19|hg38] \
    [--cytobands true|false] \
    [--results <output_directory>]
    [--regions <regions.bed>]
```

The main input is either a single BAM file, or a text file contining a list of paths to BAM files, one per line.

If `--genome` is not supplied, the version will be guessed based on the BAM header. Cytoband definitions are stored under `data`, and the version is based on either the supplied or the guessed genome version.

Coverage for custom regions can plotted by supplying a `--regions` bed file. This requires that the name column of the bed file has an entry, and one plot per unique entry in this column will be produced.

## Output

By default, the following output will be created:

```
results/
    <sample>/
        <sample>.mosdepth.global.dist.txt
        <sample>.mosdepth.summary.txt
        <sample>.per-base.d4
        plots/
            <sample>.total_coverage.png
            <sample>.<region>.png (one per region, if --regions supplied)
```
