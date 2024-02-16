# WGS clinical coverage

A Nextflow pipeline that calculates and visualises read coverage in WGS data.

## Usage

```
nextflow run gmc-norr/wgs-coverage \
    --bam <input.bam> \
    [--genome hg19|hg38] \
    [--cytobands true|false] \
    [--results <output_directory>]
```

If `--genome` is not supplied, the version will be guessed based on the BAM header. Cytoband definitions are stored under `data`, and the version is based on either the supplied or the guessed genome version.

## Input

- Alignments in BAM format

## Output

- PNG visualisation of the coverage across all chromosomes
