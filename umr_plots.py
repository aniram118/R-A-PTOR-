#!/usr/bin/env python
from __future__ import print_function

disclaimer = """
This script produces comprehensive analysis plots of UMR regions from 3' FastQ file.
"""

import sys
import os
import matplotlib

matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import datetime
import seaborn as sns
import pandas as pd
import numpy as np
import collections

plt.style.use("seaborn-colorblind")

matplotlib.rcParams['axes.titlesize'] = 13
matplotlib.rcParams['axes.labelsize'] = 12
matplotlib.rcParams['xtick.labelsize'] = 11
matplotlib.rcParams['ytick.labelsize'] = 11
matplotlib.rcParams['figure.titlesize'] = 20


def fastq_reader(filename):
    i = 0
    with open(filename, 'r') as infile:
        name = infile.readline().rstrip()
        while True:
            i += 1
            seq = ""
            for s in infile:
                if s[0] == '+':
                    commentp = s.rstrip()
                    break
                else:
                    seq += s.rstrip()
            qual = ""
            for q in infile:
                if len(qual) > 0 and q[0] == '@':
                    yield name, seq, qual
                    name = q.rstrip()
                    break
                else:
                    qual += q.rstrip()
            else:
                yield name, seq, qual
                return


def get_fastq_info(filename):
    seq_lengths = []
    mean_qualities = []
    kmers_start = []
    kmers_end = []
    nts_A = []
    nts_G = []
    nts_T = []
    nts_C = []
    nts_U = []
    # channels = []
    # start_times = []

    for head, seq, qual in fastq_reader(filename):
        name, comment = head.split(" ", 1)

        seq_lengths.append(len(seq))
        seq=seq.replace('T','U')
        mean_qualities.append(round(np.mean(bytearray(qual, "ascii")) - 33, 2))

        kmers_start.append(seq[0:4])
        kmers_end.append(seq[-4:])

        ntc = collections.Counter()
        ntc.update(seq)
        nts_A.append(ntc['A'])
        nts_G.append(ntc['G'])
        nts_T.append(ntc['T'])
        nts_C.append(ntc['C'])
        nts_U.append(ntc['U'])

        # PromethION/MinION
        # info = dict(part.split("=") for part in comment.split(" "))

        # if "ch" in info:
        #     ch = int(info["ch"])
        #     channels.append(ch)

        # if "start_time" in info:
        #     start_time = info["start_time"]
        #     t = pd.to_datetime(start_time, infer_datetime_format=False).floor('S')
        #     start_times.append(t)

    df = pd.DataFrame({"seq_length": seq_lengths,
                       "mean_quality": mean_qualities,
                       "kmers_start": kmers_start,
                       "kmers_end": kmers_end,
                       "nt_A": nts_A,
                       "nt_G": nts_G,
                       "nt_T": nts_T,
                       "nt_C": nts_C,
                       "nt_U": nts_U,
                       })
    return df


def xlabel_pos_right(ax):
    label = ax.xaxis.get_label()
    x_lab_pos, y_lab_pos = label.get_position()
    label.set_position([1.0, y_lab_pos])
    label.set_horizontalalignment('right')
    ax.xaxis.set_label(label)


## Plots

def plot_nt_content(ax, df):
    ax.set_title("Nucleotide Content", loc="left")
    ax.set_xlabel("NT")
    ax.set_ylabel("Percent")
    nt = {}
    nt['A'] = np.sum(df.nt_A)
    nt['T'] = np.sum(df.nt_T)
    nt['G'] = np.sum(df.nt_G)
    nt['C'] = np.sum(df.nt_C)
    nt['U'] = np.sum(df.nt_U)
    s = np.sum((nt['A'], nt['T'], nt['G'], nt['C'], nt['U']))
    for k, v in nt.items(): nt[k] = (v / s) * 100
    if nt['U'] > 0:
        sns.barplot(['G', 'C', 'A', 'U'], [nt['G'], nt['C'], nt['A'], nt['U']], ax=ax)
    else:
        sns.barplot(['G', 'C', 'A', 'T'], [nt['G'], nt['C'], nt['A'], nt['T']], ax=ax)


def plot_length_binned(ax, df, num_bins=1000):
    ax.set_title("Histogram of sequence lengths 1000 bins", loc="left")
    ax.set_xlabel("Lengths")
    ax.set_ylabel("Count (Log)")
    ax.hist(df.seq_length, num_bins)
    ax.set_yscale('log', nonposy='clip')
    xlabel_pos_right(ax)
    # sns.distplot(x, ax=ax, bins=num_bins, kde=False)


def plot_length_percentile(ax, df, percentile=90):
    cutoff = int(round(np.percentile(df.seq_length, percentile)))
    num_bins = int(round(cutoff))
    ax.set_title("Histogram of Sequence Lengths from 0 to " + str(cutoff) + " (" + str(percentile) + "%ile) in " + str(
        num_bins) + " bins", loc="left")
    ax.set_xlabel("Length")
    ax.set_ylabel("Count (Log)")
    ax.set_yscale('log', nonposy='clip')
    xlabel_pos_right(ax)
    ax.hist(df.seq_length, num_bins, range=(0, cutoff))


def plot_nucleotides_per_length(ax, df):
    blens = {}
    for l in df.seq_length:
        bin_l = int(l / 100) * 100
        if bin_l in blens:
            blens[bin_l] += l
        else:
            blens[bin_l] = l

    len_sum = df["seq_length"].sum()
    mean_sum = len_sum / 2.0
    n50 = len_sum
    l = 0
    for k in sorted(df["seq_length"], reverse=True):
        l += k
        if l > mean_sum:
            break
        n50 = k

    print("N50\t" + str(n50))

    ax.set_title("Nucleotides per Sequence Length in 100nt bins (N50:" + str(n50) + ")", loc="left")
    ax.set_xlabel("Length")
    ax.set_ylabel("GBases")
    lens = np.asfarray(list(blens.values()))
    ax.bar(list(blens.keys()), lens / float(1000 ** 3), width=100)  # , s=point_size) #, range=(0,20000))
    ax.axvline(x=n50, color="red")
    xlabel_pos_right(ax)


def plot_qualities(ax, df):
    ax.set_title("Histogram of Mean Qualities", loc="left")
    ax.hist(df.mean_quality, 160, range=(0, 40))  # , s=point_size) #, range=(0,20000))
    ax.set_xlim(left=0, right=40)
    ax.set_xlabel("Quality")
    xlabel_pos_right(ax)
    # df.hist(ax=ax, column=['mean_quality'], by=['seq_length'])


# Matplotlib hist2d isn't so nice with log scales, so use pcolormesh directly
def plot_hist2d_helper(ax, xbins, ybins, x, y):
    counts, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
    ax.pcolormesh(xbins, ybins, counts.T, cmap="Blues", norm=matplotlib.colors.LogNorm())


def plot_quality_seqlength(ax, df):
    ax.set_title("Read Length vs Quality", loc="left")
    max_len = np.max(df.seq_length)
    xbins = np.arange(0, 40, 0.1)
    ybins = np.arange(1, np.log(max_len), 0.01) ** 10
    plot_hist2d_helper(ax, xbins, ybins, df.mean_quality, df.seq_length)
    ax.set_xlim(left=0, right=40)
    ax.set_xlabel("Mean Quality")
    ax.set_ylabel("Read Length Log10")
    xlabel_pos_right(ax)
    ax.set_yscale('log', nonposy='clip')


# def plot_time_seqlength(ax, df):
#     ax.set_title("Read Length by Time", loc="left")
#     max_len = np.max(df.seq_length)
#     deltas = df.start_times - df.start_times.min()
#     hours = deltas / pd.Timedelta(hours=1)
#     xbins = np.arange(0, np.max(hours), 1)
#     ybins = np.arange(1, 6, 0.01) ** 10
#     plot_hist2d_helper(ax, xbins, ybins, hours, df.seq_length)
#     ax.set_xlabel("Time")
#     ax.set_ylabel("Read Length Log10")
#     xlabel_pos_right(ax)
#     ax.set_yscale('log', nonposy='clip')
#     ax.set_xlim(left=0, right=np.max(hours))
#     xlabel_pos_right(ax)


# def plot_nts_by_time(ax, df):
#     ax.set_title("Cumulative Nucleotides since start time", loc="left")
#     deltas = df.start_times - df.start_times.min()
#     hours = deltas / pd.Timedelta(hours=1)
#     ax.plot(hours, df.seq_length.cumsum() / float(1000 ** 3), linewidth=6)
#     ax.set_xlabel("Time (hours)")
#     ax.set_ylabel("GBases")
#     ax.set_xlim(left=0, right=np.max(hours))
#     xlabel_pos_right(ax)



# adapted from https://github.com/mattloose/flowcellvis
# def get_coords(channel, flowcellsize):
#     if flowcellsize == 3000:
#         # find which block of 12 we are in:
#         block = (channel - 1) // 250
#         remainder = (channel - 1) % 250
#         row = remainder // 10
#         column = remainder % 10 + block * 10
#         return (column, row)
#     if flowcellsize == 128:
#         return (channel // 8, channel % 8)
#     else:
#         return minION_flowcell_layout(channel)


# def plot_nts_by_layout(ax, df):
#     title = "Channel activity"
#     ax.set_title(title, loc="left")
#     flowcellsize = 512
#     if df.channels.max() > 512:
#         flowcellsize = 3000
#     coords = [get_coords(channel, flowcellsize) for channel in df.channels]
#     ax.scatter([c[0] for c in coords], [c[1] for c in coords], alpha=0.1, s=3)
#     ax.set_xticks([])
#     ax.set_yticks([])
#     xlabel_pos_right(ax)


def plot_kmer_start(ax, df):
    df.kmers_start.value_counts().head(n=40).plot.bar(ax=ax)
    ax.set_title("Start of Sequence k-mer Content (k=4) top 40", loc="left")
    ax.set_xlabel("K-mer")
    ax.set_ylabel("Count")
    xlabel_pos_right(ax)


def plot_kmer_end(ax, df):
    df.kmers_end.value_counts().head(n=40).plot.bar(ax=ax)
    ax.set_title("End of Sequence k-mer Content (k=4) top 40", loc="left")
    ax.set_xlabel("K-mer")
    ax.set_ylabel("Count")
    xlabel_pos_right(ax)


# Assemble Plots
def plot_fastq_info(df):
    (fig, axsa) = plt.subplots(ncols=2, nrows=2)
    axs = axsa.flatten().tolist()
    axs.reverse()
    point_size = 8
    marker_size = 3

    plot_length_binned(axs.pop(), df, 1000)

    plot_qualities(axs.pop(), df)

    # plot_nucleotides_per_length(axs.pop(), df)

    plot_quality_seqlength(axs.pop(), df)

    #  plot_length_percentile(axs.pop(), df, 90)

    #  plot_length_percentile(axs.pop(), df, 70)

    # plot_nts_by_time(axs.pop(), df)

    # plot_nts_by_channel(axs.pop(), df)
    plot_nt_content(axs.pop(), df)

    # plot_time_seqlength(axs.pop(), df)

    # plot_nts_by_layout(axs.pop(), df)

    #  plot_kmer_start(axs.pop(), df)

    #  plot_kmer_end(axs.pop(), df)

    plt.subplots_adjust(left=0.2, wspace=1.0, top=2.8)
    basename = os.path.basename(filename)
    plt.suptitle(
        basename + " (" + str(len(df)) + " reads, " + str(round(df.seq_length.sum() / (1000.0 ** 3), 2)) + "GBases)")
    fig.set_size_inches(21, 12)
    fig.tight_layout()
    plt.subplots_adjust(top=0.9)
    return plt


## MAIN ##

print(disclaimer, file=sys.stderr)

if len(sys.argv) < 2:
    print("USAGE: " + sys.argv[0] + " FASTQ_File(s)\n")

for filename in sys.argv[1:]:

    basename = os.path.basename(filename)

    df_filename = basename + ".df.tsv"
    if os.path.isfile(df_filename):  # read alreads computed dataframe
        df = pd.read_csv(df_filename, sep="\t", encoding='utf-8')
    else:
        df = get_fastq_info(filename)
        df.to_csv(df_filename, sep='\t', encoding='utf-8')

    print(df.describe())
    plt = plot_fastq_info(df)
    plt.savefig(basename + ".png")