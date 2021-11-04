from argparse import ArgumentParser
from pathlib import Path
from matplotlib import pyplot


def plot_umi_per_cell(filename, ax=None):
    umis = read_umi_per_cell(filename)

    if ax is None:
        f = pyplot.figure()
        ax = f.add_subplot(1, 1, 1)

    ax.plot(umis, range(len(umis)), linewidth=2)
    ax.set_xlabel('Cumulative number of barcodes')
    ax.set_ylabel('count of reported UMIs')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid(color="dimgrey", linestyle="-", linewidth=0.5, which="both", alpha=0.3)
    return ax


def read_umi_per_cell(filename):
    values = []
    with open(filename, "rt") as instream:
        for line in instream:
            values.append(float(line.rstrip()))

    return values


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument("--umi-per-cell")
    parser.add_argument("-o", "--output", required=True)

    args = parser.parse_args(cmdline)
    ax = None
    if args.umi_per_cell is not None:
        ax = plot_umi_per_cell(args.umi_per_cell)

    if ax is None:
        parser.error("Please select a plot to generate")

    ax.figure.savefig(args.output)


if __name__ == "__main__":
    main()
