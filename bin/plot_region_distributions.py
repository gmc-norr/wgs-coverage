#!/usr/bin/env python

import logging
from pathlib import Path
from typing import Dict, List, Union

import click
import matplotlib as mpl
import numpy as np
import polars as pl
import scipy
from matplotlib import pyplot as plt
from plot_coverage import parse_bed
from pyd4 import D4File


def plot_region_distributions(
    coverage: List[Path], regions: List[Dict[str, Union[str, int]]], gene: str
):
    chromosome = regions[0]["chrom"]
    min_pos = min([x["start"] for x in regions])
    max_pos = max([x["end"] for x in regions])

    x = np.arange(min_pos, max_pos)

    gs = mpl.gridspec.GridSpec(len(coverage), 1)
    fig = plt.figure(
        figsize=(8, len(coverage) // 3),
    )

    max_y = 0

    axs = []

    data_df = pl.DataFrame()

    for i, p in enumerate(sorted(coverage)):
        logging.debug("getting data from %s", p)
        cov = D4File(str(p))

        ax = fig.add_subplot(gs[i, :])
        axs.append(ax)

        x = np.arange(1, 100)
        y = cov.resample((chromosome, min_pos, max_pos), bin_size=1)[0]
        try:
            cov_kde = scipy.stats.gaussian_kde(y)
            dy = cov_kde.evaluate(x)
        except np.linalg.LinAlgError:
            cov_kde = None
            dy = np.zeros(len(x))

        data_df = data_df.vstack(
            pl.DataFrame(
                {
                    "sample": p.name.split(".")[0],
                    "gene": gene,
                    "x": x,
                    "y": dy,
                }
            )
        )

        ax.plot(x, dy, lw=0.5, color="black")
        ax.fill_between(x, dy, color="steelblue")

        max_y = max(max_y, max(dy))

        ax.text(-1, 0, p.name.split(".")[0], ha="right")

        for s in ["top", "right", "bottom", "left"]:
            ax.spines[s].set_visible(False)
        ax.yaxis.set_visible(False)
        ax.xaxis.set_visible(False)
        ax.patch.set_alpha(0)

    for ax in axs:
        ax.set_ylim(0, max_y)
        ax.set_xlim(1, 60)

    axs[-1].xaxis.set_visible(True)
    axs[-1].set_xlabel("Coverage")
    gs.update(hspace=-0.75)

    data_df.rechunk()

    return data_df, fig


@click.command()
@click.argument(
    "regions",
    type=click.Path(
        path_type=Path,
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
)
@click.argument(
    "coverage",
    nargs=-1,
    type=click.Path(
        path_type=Path,
        exists=True,
        file_okay=True,
        dir_okay=False,
    ),
)
@click.option(
    "-g",
    "--gene",
    "gene",
    help="Limit plotting to this gene",
)
@click.option(
    "--dpi",
    "dpi",
    type=int,
    default=300,
    show_default=True,
    help="DPI of output files",
)
@click.option(
    "-v",
    "--verbose",
    "verbose",
    is_flag=True,
    help="Enable debug logging",
)
def main(regions, coverage, gene, dpi, verbose):
    """
    When plotting the coverage distributions, zeroes are ignored.
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=log_level)

    regions = parse_bed(regions)
    for g, d in regions.items():
        if gene and gene != g:
            continue
        logging.debug("plotting %s", g)
        d, p = plot_region_distributions(coverage, d, g)
        p.savefig(
            f"{g}.distribution.png",
            bbox_inches="tight",
            dpi=dpi,
        )
        d.write_csv(f"{g}.distribution.tsv", separator="\t")
        if gene:
            break
    else:
        logging.error("could not find gene %s", gene)


if __name__ == "__main__":
    main()
