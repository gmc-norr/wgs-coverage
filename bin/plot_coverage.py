#!/usr/bin/env python

import pathlib
from collections import defaultdict

import click
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import style
from pyd4 import D4File

INCLUDE_CHROMS = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY",
]


CYTOBAND_COLORS = {
    "acen": "#AA0000",
    "gneg": "#FFFFFF",
    "gpos100": "#000000",
    "gpos75": "#404040",
    "gpos50": "#7F7F7F",
    "gpos25": "#BFBFBF",
    "gvar": "#00008B",
    "stalk": "#ADD8E6",
}


def parse_cytobands(cytobands_file) -> dict:
    cytobands = defaultdict(list)
    with open(cytobands_file) as f:
        for line in f:
            chrom, *rest = line.strip().split()
            if chrom not in INCLUDE_CHROMS:
                continue
            cytobands[chrom].append(
                {
                    "start": int(rest[0]),
                    "end": int(rest[1]),
                    "name": rest[2],
                    "stain": rest[3],
                }
            )
    return cytobands


def plot_coverage(
    coverage: D4File,
    cytobands: dict,
    bin_size: int = 1_000_000,
    n_col: int = 10,
    max_coverage: int = 100,
):
    panel_factor = 1
    plot_cytobands = len(cytobands) > 0
    n_row = 1 + len(INCLUDE_CHROMS) // n_col

    if plot_cytobands:
        n_col *= 2
        panel_factor = 2

    if plot_cytobands:
        chrom_layout = [
            x
            for pair in zip([f"{c}_cb" for c in INCLUDE_CHROMS], INCLUDE_CHROMS)
            for x in pair
        ]
    else:
        chrom_layout = INCLUDE_CHROMS[:]

    chrom_layout.extend(
        ["." for _ in range((n_row * n_col) - len(INCLUDE_CHROMS) * panel_factor)]
    )
    chrom_layout = [
        chrom_layout[i : i + n_col] for i in range(0, len(chrom_layout), n_col)
    ]

    width_ratios = np.repeat(1, n_col)
    if plot_cytobands:
        width_ratios = [x for _ in range(n_col // 2) for x in (0.15, 1)]

    cov = {}
    for c in INCLUDE_CHROMS:
        resampled_chrom_cov = coverage.resample(c, bin_size=bin_size)[0]
        bin_start = range(1, bin_size * len(resampled_chrom_cov), bin_size)
        cov[c] = {
            "coverage": resampled_chrom_cov,
            "position": bin_start,
        }

    max_cov = min(max([max(cov[c]["coverage"]) for c in cov]), max_coverage)
    max_pos = max([max(cov[c]["position"]) for c in cov])

    # style.use("ggplot")

    fig, axs = plt.subplot_mosaic(
        chrom_layout, layout="constrained", width_ratios=width_ratios
    )
    for chrom in INCLUDE_CHROMS:
        ax = axs[chrom]
        x = cov[chrom]["coverage"]
        y = cov[chrom]["position"]

        ax.axvline(x=0, color="gray", linewidth=0.5, linestyle="solid")
        ax.axvline(x=30, color="firebrick", linewidth=0.5, linestyle=(0, (5, 5)))

        ax.plot(x, y)
        ax.set_title(chrom)
        ax.set_ylim(1.05 * max_pos, -0.05 * max_pos)
        ax.set_xlim(-0.05 * max_cov, max_cov * 1.05)

        ax.yaxis.set_visible(False)
        ax.spines.left.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.spines.right.set_visible(False)
        # ax.tick_params(axis="y", which="both", length=0, labelsize=0)

        if plot_cytobands:
            ax = axs[f"{chrom}_cb"]
            ax.set_ylim(1.05 * max_pos, -0.05 * max_pos)
            ax.set_xlim(0, 10)

            ax.xaxis.set_visible(False)
            ax.spines.top.set_visible(False)
            ax.spines.right.set_visible(False)
            ax.spines.bottom.set_visible(False)

            for band in cytobands[chrom]:
                if chrom == "chr1":
                    print(band)
                rect = mpatches.Rectangle(
                    (0, band["start"]),
                    10,
                    band["end"] - band["start"],
                    facecolor=CYTOBAND_COLORS[band["stain"]],
                )
                ax.add_patch(rect)

    if plot_cytobands:
        fig.set_size_inches(20, 10)
    else:
        fig.set_size_inches(15, 10)

    fig.supxlabel("Coverage")
    fig.supylabel("Position")

    return fig


@click.command()
@click.argument("coverage")
@click.option(
    "-c",
    "--cytobands",
    "cytobands_file",
    type=click.Path(path_type=pathlib.Path),
    help="Cytoband file",
)
@click.option(
    "-o",
    "--output",
    "plot_file",
    type=click.Path(path_type=pathlib.Path),
    help="Output image file",
)
@click.option(
    "-r",
    "--regions",
    "regions",
    type=click.Path(path_type=pathlib.Path),
    help="BED file(s) with regions to cover",
)
@click.option("--dpi", default=100, show_default=True, help="Plot resolution")
def main(coverage, cytobands_file, regions, dpi, plot_file):
    f = D4File(coverage)

    cytobands = {}
    if cytobands_file:
        cytobands = parse_cytobands(cytobands_file)

    p = plot_coverage(f, cytobands)
    p.savefig(plot_file, dpi=dpi)


if __name__ == "__main__":
    main()
